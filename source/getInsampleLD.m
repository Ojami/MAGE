function ld = getInsampleLD(variants, chr, opts)
% gets in-sample LD structure (Pearson R) for input variants.
% This function doesn't filter bi-allelic variants (should be considered
% for future developments) at the moment. So, it's important (depends) that
% 'variants' are unique (function doesn't check).
% 
% Oveis Jamialahmadi, University of Gothenburg, Feb 2022.
% 
% @21APR2023: 'maf' filter was added. variants below this cut-off (default
%            is 0.01) will be removed. 
% @09SEP2024: 'chunk' was included as user defined option (by default is
%            determined automatically).

arguments
    variants {mustBeText, mustBeVector} % either variant names or in format of pos-pos
    chr {mustBeVector}
    opts.ldsample double {mustBeVector}
    opts.bgenhome {mustBeFolder}
    opts.parallel (1,1) logical = false
    opts.a2freq (1,1) logical = true % calculate freq of a2 allele
    opts.matchid {mustBeText, mustBeVector} % list of variant ids (chr:pos:a1:a2, a2 is effect allele) for matching. must have the same length as variants
    opts.removenan (1,1) logical = true % remove missing genotypes (no need for UKBB bgen files)
    opts.infoscore (1,1) double {mustBeInRange(opts.infoscore, 1e-6, 1)} = 1.1e-6 % INFO score cut-off: variants < INFO score are removed.
    opts.maf (1,1) double {mustBeInRange(opts.maf, 1e-30, 0.5)} = 0.01
    opts.chunk (1,1) double

    %@18NOV2024
    opts.range (1,1) logical = false % input variants are range of pos1-pos2?
end

if isfield(opts, 'matchid')
    opts.matchid = string(opts.matchid);
    ucolon = unique(count(opts.matchid, ":"));
    if numel(ucolon) > 1 || ~all(ucolon == 3) || numel(opts.matchid) ~= numel(variants)
        error('getInsampleLD:wrongArgument', 'matchid must be in format of chr:pos:a1:a2 where a2 is effect allele')
    end
else
    opts.matchid = [];
end

if ~isempty(duplicates(variants)) && ~isfield(opts, 'matchid')
    fprintf("\ngetInsampleLD::warning: duplicates was found in variants\n")
    fprintf("\tThis results in multiple reading of the same variant\n")
    fprintf("\tas duplicates can be assigned to different chunks on parallel pool\n")
    fprintf("\tremove duplicates since in the end they'll be read twice anyway!\n")
end

if ~isfield(opts, 'ldsample')
    opts.ldsample = getQCEID(3, false); % unrelated white-British from UK Biobank
end

if isfield(opts, "chunk")
    chunksz = round(opts.chunk);
else
    if numel(variants) < 1e3
        chunksz = inf; % use internal parfor for bgenreader
    else
        chunksz = round(numel(variants)/70);
    end
end

if ~isfield(opts, 'bgenhome')
    ld = getbulkgeno(variants, string(chr), 'parallel', opts.parallel, ...
        'datatype', "single", 'verbose', false, 'chunk', chunksz,...
        'eid', opts.ldsample, "merge", true, "infoscore", opts.infoscore);
else
    ld = getbulkgeno(variants, string(chr), 'parallel', opts.parallel, ...
        'datatype', "single", 'home', opts.bgenhome, 'verbose', false, ...
        'chunk', chunksz, 'eid', opts.ldsample, "merge", true, "infoscore", opts.infoscore);
end

% @18NOV2024
if opts.range
    variants = ld.snp;
end

if opts.removenan
    fnan = any(isnan(ld.bed), 2);
    ld.bed(fnan, :) = [];
    ld.eid(fnan) = [];
end

if isempty(opts.matchid)
    % change to original format since bgenreader by default flips a1/a2 so that
    % a2 is always minor allele (should be changed). This affects r (for
    % fine-mapping and colocalization) and not r2 (for regional plots).
    ld.bed(:, ld.swap) = 2 - ld.bed(:, ld.swap); 
    newa2 = ld.a1(ld.swap);
    ld.a1(ld.swap) = ld.a2(ld.swap);
    ld.a2(ld.swap) = newa2;

else
    % match to custom 'matchid' variant ids
    ld.id = join([string(ld.chr), ld.pos, ld.a1, ld.a2], ":");
    ld.idf = join([string(ld.chr), ld.pos, ld.a2, ld.a1], ":"); % flip ids

    idxf = ismember(ld.idf, opts.matchid); 
    if any(idxf) % flip a1/a2
        ld.id(idxf) = ld.idf(idxf);
        ld.bed(:, idxf) = 2 - ld.bed(:, idxf); 
        newa2 = ld.a1(idxf);
        ld.a1(idxf) = ld.a2(idxf);
        ld.a2(idxf) = newa2;
    end
    ld = rmfield(ld, "idf");

    % remove non-matching ids from calls
    rmidx = ~ismember(ld.id, opts.matchid);
    if any(rmidx)
        ld.bed(:, rmidx) = [];
        fis = setdiff(fieldnames(ld), ["eid", "bed"]);
        for i = 1:numel(fis), ld.(fis{i})(rmidx) = []; end
    end

    % set orders same as 'matchid' field
    [idx1, idx2] = ismember(opts.matchid, ld.id);
    ld.bed = ld.bed(:, idx2(idx1));
    fis = setdiff(fieldnames(ld), ["eid", "bed"]);
    for i = 1:numel(fis), ld.(fis{i}) = ld.(fis{i})(idx2(idx1)); end
end

idx_missing = ~(ismember(ld.snp, variants) | ismember(ld.pos + "-" + ld.pos, variants));

% apply MAF cut-off
ld.a2freq = (mean(ld.bed, 1)./2).';
maf = ld.a2freq;
idx = maf > 0.5;
maf(idx) = 1 - maf(idx);
idx = (maf < opts.maf) | idx_missing;

if any(idx)
    fi = setdiff(fieldnames(ld), ["eid", "bed"]);
    ld.bed(:, idx) = [];
    for i = 1:numel(fi), ld.(fi(i))(idx) = []; end
end

if ~opts.a2freq % this flag should be removed
    ld = rmfield(opts, "a2freq");
end
ld.ld = corr(ld.bed);

fi = ["eid", "swap", "bed", "eid"];
ld = rmfield(ld, fi);
 
end % END