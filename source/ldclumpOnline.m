function [R, C] = ldclumpOnline(rsid, pval, p)
% does LD-clumping using IEU OpenGWAS API as documented here:
% http://gwas-api.mrcieu.ac.uk/docs/ Since there may be some variants
% missing from the LD reference panel used by the API, LD clumping will be
% performed for those variants using LD structure of unrelated white
% British in UKBB (so pop must be set to EUR and is only valid for UKBB).
% 
% Caution: this clumping approach may be wrong specially when index
%          variants are mistakenly determined (due to presence of missing
%          index variants from the LD ref panel). 
% Output:
%   R: a string array containing top physically independent signals
%   C: a string array of available ids in the LD ref panel.
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, Sep 2021.

arguments
    rsid {mustBeVector, mustBeText}
    pval {mustBeVector, mustBeNonempty, mustBeInRange(pval, 0, 1)}
    p.pthresh (1,1) double = 5e-8
    p.r2 (1,1) double = 0.01
    p.kb (1,1) double = 1000
    p.pop (1,1) {mustBeMember(p.pop, ["EUR", "SAS", "EAS", "AFR", "AMR", "legacy"])} = "EUR"
    p.reflookup (1,1) logical = true % check if input ids are present in the LD ref panel. If some ids are missing, then index signals will be further compared with these missing ids internally with UKBB data.
    p.pos {mustBeVector, mustBeNumeric} = [];
    p.chr {mustBeVector} = [];
    p.minorAllele {mustBeVector} = "";
    p.parallel (1,1) logical = false % for internal LD calculations
end

headers = {'Content-Type' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', 50000);

% keep only variants with p <= pthresh
idx = pval > p.pthresh;
rsid(idx) = []; pval(idx) = [];
if all(~isempty(p.pos))
    p.pos(idx) = [];
    p.chr(idx) = [];
end
if all(p.minorAllele ~= "")
    p.minorAllele(idx) = [];
end

if p.reflookup
    fprintf('finding missing variants from LD reference panel...')
    cdata = struct('rsid', rsid, 'pop', p.pop);
    C = string(webwrite("http://gwas-api.mrcieu.ac.uk/ld/reflookup", cdata, options));
    fprintf('\b\b Done.\n')
else
    C = [];
end

data = struct('rsid', rsid, 'pval', pval, 'pthresh', p.pthresh, ...
    'r2', p.r2, 'kb', p.kb, 'pop', p.pop);

url = "http://gwas-api.mrcieu.ac.uk/ld/clump";
fprintf('performing online clumping...')
R = webwrite(url, data, options);
R = string(R);
fprintf('\b\b Done.\n\n')

% check if some variants are missing from input variants, and further use
% individuals level data from UKBB to infer their clump. In this case both
% pos and chr options must be provided. For multi-allelic variants,
% minorAllele should also be provided to avoid index errors. 
if p.reflookup && ~isempty(C)
    fprintf('clumping the missing variants from LD reference panel...\n')

    qceid = getQCEID(3); % white unrelated Brits
    
    % build two tables for current index variants and those missing from LD
    % ref panel
    idx = ismember(rsid, R); 
    cX = table(rsid(idx), pval(idx), p.chr(idx), p.pos(idx), ...
        'VariableNames', {'snp', 'p', 'chr', 'pos'});
    if height(cX) ~= numel(R) % multi-allelic index variants
        fprintf('WARNING: some of index variants are multi-allelic!\n')
    end
    [~, idx2] = setdiff(rsid, C);
    mX =  table(rsid(idx2), pval(idx2), p.chr(idx2), p.pos(idx2), ...
        'VariableNames', {'snp', 'p', 'chr', 'pos'});
    
    if all(p.minorAllele ~= "")
        cX.minorAllele = p.minorAllele(idx);
        mX.minorAllele = p.minorAllele(idx2);
    end
    clear idx idx2
    
    % compare each index variant with missing variants from LD ref panel to see
    % if any of them lie within -/+ kb margin from the index signal with a
    % smaller p-value
    for i = 1:size(cX, 1)
        idx = mX.pos >= (cX.pos(i) - p.kb*1e3) & ...
            mX.pos <= (cX.pos(i) + p.kb*1e3) & ...
            mX.chr == cX.chr(i) & ...
            mX.p < cX.p(i);

        if any(idx)
            idx = find(idx);
            snp = [cX.snp(i); mX.snp(idx)];
            chr = string([cX.chr(i); mX.chr(idx)]);
            if all(p.minorAllele ~= "")
                geno = getbulkgeno(snp, chr, 'minorAllele', ...
                    [cX.minorAllele(i); mX.minorAllele(idx)],...
                    'verbose', false, 'parallel', p.parallel);
            else
                geno = getbulkgeno(snp, chr, 'verbose', false, 'parallel', p.parallel);
            end
            
            fi = fieldnames(geno);
            qcidx = ismember(geno.(fi{1}).eid, qceid);
            [~, snpIdx] = ismember(snp, geno.(fi{1}).snp);
            R2 = corr(geno.(fi{1}).bed(qcidx, snpIdx)).^2;
            R2(2:end, :) = []; % only keep the LD r2 with current index variant
            idx = idx(R2(2:end) >= p.r2); % find variants in LD with current index variant
            if any(idx)
                newindex = mX(idx, :);
                [~, idx] = min(newindex.p); % strongest association in this locus
                fprintf('%s replaced current index (%s)\n', newindex.snp(idx), cX.snp(i))
                cX(i, :) = newindex(idx, :); % new index variant in this locus
            end

        end
    end
    
    R = cX.snp;
    fprintf('\b\b Done.\n\n')
end

end %END