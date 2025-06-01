function clumped = ldclumpPeaks(gss, opts)
% uses MATLAB findpeaks to identify approximately independent loci.

arguments
    gss
    opts.sig (1,1) double = 5e-8 % significance threshold (can be as -log10 P)
    opts.log10p (1,1) logical = false % p values are in -log10 scale so as 'sig'
    opts.light (1,1) logical = true % light version of Manhattan plot by excluding variants with P < 0.05
    opts.parallel (1,1) logical = false % use tall/gather for reading 'infile' summary stat file
    opts.verbose (1,1) logical = true
    opts.win (1,1) double = 1e4 % window for clumping in bp
    opts.n (1,1) double = 2e4 % number of randomly selected individuals for pairwise LD
    opts.r2 (1,1) double = 0.01 % same as PLINK clump r2
    opts.bgenhome {mustBeFolder} 
    opts.bgenSexSubfolder {mustBeTextScalar} = "sex.chrom"% subfolder within 'bgenhome' where sex chromosomes are present, if left empty, looks for files withi 'bgenhome' directory.
    opts.refGenome {mustBeTextScalar, mustBeMember(opts.refGenome, ["GRCh37", "GRCh38"])} =  "GRCh37" % for nearest gene
    opts.checkld (1,1) logical = true
end

% gss can be a table or a file
if ~istable(gss)
    gss = readGWASfile(gss, 'parallel', opts.parallel, 'full', false, ...
        'light', opts.light, 'p', 0.05);
else
    gss.pos = double(gss.pos);
end

cols = colnames(gss);
gss.Properties.VariableNames = cols.lower;
    
% reconcile zero p-values passing machine precision
if opts.log10p
    f_zero = isinf(gss.p);
    if any(f_zero)
        gss.p(f_zero) = -log10(realmin*eps);
    end
else
    f_zero = gss.p == 0;
    if any(f_zero)
        gss.p(f_zero) = realmin*eps;
    end
end

if ~opts.log10p
    if any(gss.p <= 0) || any(gss.p > 1)
        error('P values must be in [0, 1)!')
    end
    
    gss.p = -log10(gss.p); % log10 transofrmation of p-values
    opts.sig = -log10(opts.sig);
end

if opts.verbose
    fprintf('finding peaks (loci)...')
    labtime = tic;
end

gss = sortrows(gss, "pos", "ascend"); % findpeaks needs strictly increasing points
if isstring(gss.chr)
    [~, idx] = natsort(gss.chr);
    gss = gss(idx, :);
else
    gss = sortrows(gss, "chr", "ascend");
end

g = groupsummary(gss, "chr", @max, "pos");

if isstring(g.chr)
    [~, idx] = natsort(g.chr);
    g = g(idx, :);
end

g.sum = cumsum(g.fun1_pos);
gss.posx = gss.pos;
chrU = g.chr;
for i = 2:numel(chrU)
    idx = gss.chr == chrU(i);
    gss.posx(idx) = gss.posx(idx) + g.sum(i-1) + 1e9;
end
gss.id = string(gss.chr) + ":" + gss.pos;
[~, idx] = unique(gss.id, "stable");
gss = gss(idx, :);

% find sex chromosomes
sexidx = ismember(string(gss.chr).lower, ["x", "y", "23", "24"]);

clumped = peakClump(gss(~sexidx, :), opts);

if any(sexidx) % fetch sex chr LD
    if isfield(opts, "bgenSexSubfolder")
        opts.bgenhome2 = fullfile(opts.bgenhome, opts.bgenSexSubfolder);
    end
    clumped = [clumped; peakClump(gss(sexidx, :), opts)];
end

if opts.verbose
    fprintf('\b\b done (%.2f sec)\n', toc(labtime))
end

end % END

%% subfunctions ===========================================================
function gss = peakClump(gss, opts)

[~, posx] = findpeaks(gss.p, gss.posx, MinPeakDistance=opts.win, MinPeakHeight=opts.sig);
idx = find(ismember(gss.posx, posx));

% findpeaks needs two neighbouring points for each peak, meaning at least 3
% snps per each chr should be present.
chrn = groupsummary(gss, "chr");
chrn(chrn.GroupCount > 2, :) = [];
if ~isempty(chrn)
    tmp = cell(height(chrn), 1);
    for k = 1:height(chrn)
        tmp{k, 1} = gss.posx(gss.chr == chrn.chr(k) & gss.p >= opts.sig);
    end
    tmp = vertcat(tmp{:});
    idx2 = find(ismember(gss.posx, tmp));
    idx = union(idx, idx2, "stable");
end

% find variants in LD
if opts.checkld
    checkld = gss(idx, :);
    checkld.keep = true(height(checkld), 1);
    eid = datasample(getQCEID(3, false), opts.n); % n randomly selected white British
    
    snpcol = find(ismember(colnames(checkld).lower, ["annotate", "snp"]), 1);
    
    ld = getInsampleLD(checkld.(snpcol), checkld.chr, ...
            "ldsample", eid, "parallel", opts.parallel, "bgenhome", opts.bgenhome);
    
    % first sort p, so that filtering be done based on strongest hit
    % first, and then checking weaker associations.
    [checkld, sortedIdx] = sortrows(checkld, "p", "descend");
    idx = idx(sortedIdx);
    
    checkld = checkKeepLD(checkld, ld, snpcol, opts);
    
    idx(~checkld.keep) = [];
else
    snpcol = find(ismember(colnames(gss).lower, ["annotate", "snp"]), 1);
end

gss.locus = repmat("", height(gss), 1);
gss.locus(idx) = gss.(snpcol)(idx);
gss(:, ["id", "posx"]) = [];
emptyLocus = gss.locus == "";

vep = mapVariant2gene(gss(~emptyLocus, :), 'verbose', true, ...
    'bgenhome', opts.bgenhome, 'refGenome', opts.refGenome);
[f1, f2] = ismember(gss.locus, vep.id); f2(f2<1) = [];
gss.locus(f1) = vep.nearestGene(f2);

% gss(gss.locus == "", :) = []; % return full summary stats for manplotter

end

%% ------------------------------------------------------------------------
function checkld = checkKeepLD(checkld, ld, snpcol, opts)

    [f1, f2] = ismember(checkld.(snpcol), ld.snp); f2(f2<1) = [];
    checkld(~f1, :) = [];
    ld.ld = ld.ld.^2;
    ld.ld = ld.ld(f2, :); ld.ld = ld.ld(:, f2);
    ld.chr = ld.chr(f2);
    uchr = unique(ld.chr);
    ld.r2 = false(size(ld.ld));
    for j = 1:numel(uchr)
        chridx = ld.chr == uchr(j);
        ld.r2(chridx, chridx) = ld.ld(chridx, chridx) >= opts.r2;
    end
    ld = ld.r2;

    for j = 1:size(ld, 2)
        if any(ld(:, j))
            if ~checkld.keep(j)
                continue % this hit is not a lead within the clump, so we don't care to check for it's LD friends
            end

            checkld.keep(ld(:, j) & checkld.p < max(checkld.p(ld(:, j)))) = false;
            checkld.p(~checkld.keep) = 0; % this avoids ignoring snp1 (p = 1e-8), in a situation where snp2 (p = 1e-9) is in ld with snp3 (p = 1e-10), and snp1 is in ld with snp2.
            ld(:, j) = false; ld(j, :) = false; % don't include this hit anymore after this loop
        end
    end

end