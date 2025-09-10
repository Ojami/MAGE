function makeFinemaptable(opts)
% uses finemapping results from with final clumped loci and creates a table
% based on highest PIPs per each clumped locus

arguments
    opts.fm_path {mustBeFolder} % path to susie sumstat files
    opts.fm_suffix {mustBeText, mustBeVector} = ".susie.txt" % default suffix of call_susie function. Can be output of susie-inf pipeline too
    opts.out_file {mustBeTextScalar} = "finemap" % path/to/output excel file
    opts.index {mustBeA(opts.index, "table")} % lead variants within each region 
    opts.PIPmax (1,1) double = 0.05 % max PIP to keep (to avoid CS with large set of variants)
    opts.refGenome {mustBeMember(opts.refGenome, ["GRCh37", "GRCh38"])} = "GRCh37" % for variant annotation
    opts.win (1,1) double % windows size used for fine-mapping. Default is empty, automatical detection (needs genome-wide regions)
    opts.tag {mustBeText, mustBeVector, mustBeNonempty} % optional tags for different finemappers (by default inferred autmatically)
    opts.log10BF (1,1) double = 2 % remove any CS with log10BF <= this value (see https://www.nature.com/articles/s41586-022-05473-8#Sec11)
    opts.pthresh (1,1) double = 5e-8 % same as above from FINNGEN papaer: an independent CS must meet both conditions above
    opts.keep_functional (1,1) logical = false % keep functional (pLoF + missense) variants irrespective of other filters? (under development)

    % FLAMES/CALDERA options (to create credible set files for FLAMES pipeline)
    opts.flames (1,1) logical = false
    opts.caldera (1,1) logical = false
    opts.cs_path {mustBeTextScalar} = fullfile(pwd, "cs_files") % output path for CS files of FLAMES/CALDERA
end

% check FLAMES/CALDERA args
if opts.flames || opts.caldera
    assert(isscalar(opts.fm_path), "fm_path must be scalar when flames/caldera flags is on!")
    
    if ~isfolder(opts.cs_path), mkdir(opts.cs_path); end
    
    % run once for FLAMES and another time for CALDERA
    is_caldera = opts.caldera;

    if opts.flames
        ifile = fullfile(fileparts(opts.cs_path), "indexfile.txt");
        if isfile(ifile)
            fprintf("\tFLAMES indexfile already exists: %s\n", ifile)
        else
            opts.caldera = false;
            MakeSingleFMtable(opts);
        end
    end
    
    opts.caldera = is_caldera;
    opts.flames = false;
    if opts.caldera

        ifile = fullfile(opts.cs_path, "cs_file.txt");
        if isfile(ifile)
            fprintf("\tCALDERA cs_file already exists: %s\n", ifile)
        else
            MakeSingleFMtable(opts);
        end
    end
    
    return
end

% if fm_suffix contains different fine-mapping outputs (SuSiE and FINEMAP),
% the final table should be merged
if numel(opts.fm_path) > 1
    assert(numel(opts.fm_path) == numel(opts.fm_suffix), "mismatch size of fine-mapping path/suffix")
end

df = cell(numel(opts.fm_suffix), 1);

if isfield(opts, "tag")
    colTags = string(opts.tag);
    assert(numel(opts.tag) == numel(opts.fm_suffix))
else
    opts.tag = strings(size(df));
    colTags = [];
end

% functional variants to be kept?
opts.cats = ["stop_gained"; "start_lost"; "splice_donor_variant"; ...
            "frameshift_variant"; "splice_acceptor_variant"; "stop_lost";...
            "transcript_ablation"; "missense_variant"; ...
            "protein_altering_variant"; "inframe_insertion"; "inframe_deletion"];
if opts.keep_functional
    fprintf("\t[INFO]: pLoF/missense variants with PIP >= %.2f will be kept\n", opts.PIPmax)
end

for k = 1:numel(opts.fm_suffix)
    inopts = opts;
    inopts.fm_suffix = inopts.fm_suffix(k);

    if numel(opts.fm_path) > 1
        inopts.fm_path = opts.fm_path(k);
    end
    
    fprintf("[%d of %d] Parsing fine-mapping results: %s\n", ...
        k, numel(opts.fm_suffix), inopts.fm_path)
    [df{k}, opts.tag(k)] = MakeSingleFMtable(inopts);
    fprintf("===========================================================\n")
end

if isempty(colTags)
    colTags = opts.tag;
end

fprintf("[INFO]: writing fine-mapping results to %s", opts.out_file)
if numel(df) > 1
    res = mergeFinemapTabs(df, colTags);
else
    res = df{1};
end

% write to output xlsx table
opts.out_file = regexprep(opts.out_file, [".xlsx$", ".xls$"], "");
[pth, opts.out_file, ext] = fileparts(opts.out_file);
opts.out_file = opts.out_file + ext;
if isempty(pth) || pth == ""
    pth = pwd();
elseif ~isfolder(pth)
   mkdir(pth)
end
opts.out_file = fullfile(pth, opts.out_file + ".xlsx");
writetable(res, opts.out_file, Sheet="overall");

if numel(df) > 1
    colTags = matlab.lang.makeValidName(colTags);
    for k = 1:numel(df)
        writetable(df{k}, opts.out_file, Sheet=colTags(k));
    end
end

disp(". Done.")

end % END

%% subfunctions ===========================================================
function [res, method] = MakeSingleFMtable(opts)

% get susie summary file names
sfiles = getfilenames(opts.fm_path, fullpath=false).x_;
sfiles(~sfiles.endsWith(opts.fm_suffix)) = [];
assert(~isempty(sfiles), "cannot find any fine-mapping summary stats!")

% check if they come from susie-inf (bgz/gz)
gz_idx = sfiles.endsWith([".bgz", ".gz"]);
opts.delete = [];
if any(gz_idx)
    opts.fm_suffix = opts.fm_suffix.erase([".bgz", ".gz"] + textBoundary("end"));
    gz_sfiles = sfiles(gz_idx);
    arrayfun(@(x)gunzip(x, opts.fm_path), fullfile(opts.fm_path, gz_sfiles))
    opts.delete = gz_sfiles.erase([".bgz", ".gz"] + textBoundary("end"));
    sfiles(gz_idx) = opts.delete;
    opts.delete = fullfile(opts.fm_path, opts.delete);
end

% unify colnames of index table
opts.index.Properties.VariableNames = createGWASheader(colnames(opts.index));

% check regions <-> index variants mapping --------------------------------
regions = sfiles.extractBefore(opts.fm_suffix);
patt = "chr" + digitsPattern + "_" + digitsPattern + "_" + digitsPattern;
assert(all(regions.contains(patt)), "regions should be named as chrN_pos1_pos2!")

% clean up regions
regions = arrayfun(@(x)strsplit(x, "_"), ...
    regions.erase(textBoundary("start") + "chr"), uni=false);
regions = array2table(vertcat(regions{:}), VariableNames=["chr", "start", "end"]);
regions = convertvars(regions, "chr", underlyingType(opts.index.CHR));
regions = convertvars(regions, ["start", "end"], @double);

if isfield(opts, "win")
else
    % what was the window size for fine-maping?
    regions.win = (regions.end - regions.start)./2; 
    opts.win = mode(regions.win); % note: some regions may start from 1 so 'win' may be not unique.
    fprintf("\t[INFO]: discoved a default window of %.2f mb\n", opts.win/1e6)
end

% center point: end - win
regions.center = regions.end - opts.win;

% first check: are they default MCH/similar regions or overallping regions
% were merged for fine-mapping?
[out, index] = mergeRegions(opts.index, opts.win);
merged_regions = find(cellfun(@(x)numel(x), out.idx) > 1);
if ~isempty(merged_regions) % regions were merged
    fprintf("\t[INFO] %d regions were merged during fine-mapping.\n", numel(merged_regions))
    
    % create a new map for regions
    regions.map = regions.chr + ":" + regions.start + ":" + regions.end;
    
    % modify index
    append_index = cell(numel(merged_regions), 1);
    for k = 1:numel(merged_regions)
        df1 = out(merged_regions(k), :);
        midx = ismember(index.idx, df1.idx{1});
        df2 = index(midx, :);
        index(midx, :) = [];
        
        % double check
        assert(min(df2.pos1) == df1.pos1 & max(df2.pos2) == df1.pos2, ...
            "merged regions don't match! check mergeRegions subfunction!")

        % check if this region maps to any of fine-mapped regions, if not
        % it is missing from fine-mapping probably due to MHC exclusion. We
        % reconcile it by removing those regions
        tmp_map = df1.CHR + ":" + df1.pos1 + ":" + df1.pos2;
        if ~any(regions.map == tmp_map)
            % what region do we have in fine-mapping that matches df1?
            idx1 = find(ismember(regions.start, df2.pos1) & regions.chr == df2.CHR(1));
            idx2 = find(ismember(regions.end, df2.pos2) & regions.chr == df2.CHR(1));
            if isscalar(idx1) && isscalar(idx2) && (idx1 == idx2)
                r_map = regions(idx1, ["start", "end"]);

                % find existing loci in index table mapping to this region
                idx1 = df2.pos1 == r_map.start;
                idx2 = df2.pos2 == r_map.end;
                keep_rows_idx = idx1 | idx2;
                df2 = df2(keep_rows_idx, :);
            end              
        end

        df2 = sortrows(df2, "P", "ascend"); % top locus/SNP comes first
        df2.SNP(:) = join(df2.SNP, ",");
        df2.locus(:) = join(df2.locus, ",");
        df2.pos1(:) = min(df2.pos1); df2.pos2(:) = max(df2.pos2(:));
        df2.idx = [];
        append_index{k} = df2(1, :);
        append_index{k}.idx = df1.idx;
    end
    index.idx = num2cell(index.idx);
    index = [index; vertcat(append_index{:})];
    
    index.map = index.CHR + ":" + index.pos1 + ":" + index.pos2;
    opts.index = index; clear index
    assert(all(ismember(regions.map, opts.index.map)), "index/fine-mapping merged regions mismatch!")
    idx_miss = ~ismember(opts.index.map, regions.map);

else
    regions.map = regions.chr + ":" + regions.center; 

    % map to index variants
    opts.index.map = opts.index.CHR + ":" + opts.index.BP;
    
    assert(all(ismember(regions.map, opts.index.map)), "index/fine-mapping regions mismatch!")
    
    % check what's not covered in fine-mapping
    idx_miss = ~ismember(opts.index.map, regions.map);
end


if any(idx_miss)
    fprintf("\t[missed in fine-mapping] %d region(s):\n", nnz(idx_miss))
    tmp = opts.index(idx_miss, ["CHR", "BP", "SNP", "locus"]);
    tmp = convertvars(tmp, 1:width(tmp), @string);

    wCHR = max([tmp.CHR.strlength; 3]);
    wPOS = max([tmp.BP.strlength; 3]);
    wSNP = max(tmp.SNP.strlength);
    wLOC = max(tmp.locus.strlength);

    fmtRow = "%-" + wCHR + "s  %-" + wPOS + "s  %-" + wSNP + "s  %-" + wLOC + "s\n";

    % Print header
    fprintf(fmtRow, "CHR", "POS", "SNP", "LOCUS"); %#ok

    % Print rows (vectorized with compose for clarity)
    lines = compose(fmtRow, tmp.CHR, tmp.BP, tmp.SNP, tmp.locus);
    fprintf("%s", lines);
    opts.index(idx_miss, :) = [];
end

% merge regions/fine-mapping files into index table
[~, keep_idx] = ismember(regions.map, opts.index.map);
opts.index = opts.index(keep_idx, :);

index = opts.index;
index.files = fullfile(opts.fm_path, sfiles);
index.regions = regions.chr + ":" + regions.start + "-" + regions.end;

% read/prune fine-mapping results -----------------------------------------
% read each region and remove:
% -CS -1 unless 
%   - 'index' lead variant is there
%   - Keep the max PIP -> diagnositcs due to purity filter
%   - keep_functional flag is true (as long as PIP >= PIPmax)
% -PIP < PIPmax -> avoid overcrowding CSs
% 

% decide on method/model/tool: currently only susieR (call_susie.R
% function part of MAGE), susie-inf and finemap-inf
opts.method = "SuSiE"; % default

if opts.flames || opts.caldera
    opts.PIPmax = 0;
    fprintf("\t[INFO]:FLAMES|CALDERA flag is on. Set PIPmax to 0 (keeping all variants in CS)\n")
end

opts.csfilter = 0;

res = cell(numel(index.files), 1);
for k = 1:height(index.files)
    df_opts = detectImportOptions(index.files(k), ...
        FileType="text", ...
        TextType="string", ...
        NumHeaderLines=0, ...
        VariableNamingRule="preserve", ...
        TreatAsMissing="NA");
    
    keep_cols = string(df_opts.SelectedVariableNames);
    
    if k == 1 
        if any(keep_cols.startsWith("lbf_variable"))
            opts.method = "SuSiE-INF";
        elseif any(keep_cols == "post_mean_cond")
            opts.method = "FINEMAP-INF";
        end
        fprintf("\t[fine-mapping method] %s\n", opts.method)
        opts.method = lower(opts.method);
    end

    if opts.method == "susie"
        null_cols = ["omega", "lbf_variable", "tausq", "sigmasq"];
        
        if ~(opts.flames || opts.caldera)
            % keep alpha for susie and use them. Note that sum(PIPs) per CS
            % can be > 1 and causes warning in FLAMES. This happens because
            % FLAMES originally was developed using FINEMAP (assumption of
            % 1 causal variant). So, using alphaX where X is the CS is more
            % precise/relevant in this case
            null_cols = union("alpha", null_cols);
        end
        
        keep_cols(keep_cols.startsWith(null_cols)) = [];
        keep_cols(keep_cols.startsWith("mu" + digitsPattern)) = [];
    else
        keep_cols(keep_cols.startsWith(["mu" + digitsPattern, ...
            "omega", "tausq", "sigmasq"])) = [];
    end
    
    keep_cols = setdiff(keep_cols, ["Z", "cs_specific_prob", ...
            "BETA", "SE", "beta", "se", "post_mean_cond", "post_sd_cond"], "stable");
    df_opts.SelectedVariableNames = keep_cols;
    df = readtable(index.files(k), df_opts);

    % Log10 bayes factor of comparing the solution of this model (cs
    % independent credible sets) to cs -1 credible sets
    % https://github.com/FINNGEN/finemapping-pipeline
    df = susieinf_cs_log10bf(df);
    if ~isnumeric(df.CS_LOG10BF), df.CS_LOG10BF = double(df.CS_LOG10BF); end

    % remove alpha/lbf_variable cols
    keep_cols = colnames(df);

    final_null_cols = ["alpha", "omega", "lbf_variable", "tausq", "sigmasq", "SNPVAR", "idx"];

    if opts.flames || opts.caldera
        % see above note: alphaX over PIP
        final_null_cols(final_null_cols == "alpha") = [];
    end

    keep_cols(keep_cols.startsWith(final_null_cols)) = [];
    keep_cols(keep_cols.startsWith("mu" + digitsPattern)) = [];
    
    df = df(:, keep_cols);
    if ~any(colnames(df) == "lead_r2")
        df.lead_r2(:) = nan;
    else
        df.lead_r2 = double(df.lead_r2);
    end

    if opts.method ~= "susie"
        old_cols = ["chromosome", "position", "rsid", "allele1", ...
            "allele2", "p", "prob", "post_mean", "cs"];
        new_cols = ["CHR", "POS", ...
            "SNP", "A1", "A2", "P", "PIP", "BETA_MEAN", "CS"];

        if opts.method == "finemap-inf"
            old_cols(end) = []; new_cols(end) = [];
        end
        df = renamevars(df, old_cols, new_cols);
    end

    % if any of valid CS (>0) have PIPs < PIPmax, we keep only 1 row, just
    % not to exclude that CS from final table
    [G, ~] = findgroups(df.CS);
    rowN = splitapply(@(x) {x}, (1:height(df))', G);
    max_pip_idx2 = cellfun(@(ix, x) ix(x), rowN, ...
                  num2cell(splitapply(@(x) find(x == max(x), 1), df.PIP, G)));
    max_pip_idx = false(height(df), 1);
    max_pip_idx(max_pip_idx2) = true;

    
    if opts.method ~= "finemap-inf"
        %@31AUG2025: finemap-inf doesn't generate CS
        cs_size = groupsummary(df, "CS");
        cs_size(cs_size.CS < opts.csfilter, :) = [];

        % remove negative CS (keep lead variant at this region)
        % max_pip_idx = df.PIP == max(df.PIP);
        keep_idx = (df.CS > opts.csfilter & df.PIP >= opts.PIPmax) | ...
            ismember(df.SNP, index.SNP(k).split(",")) | max_pip_idx;
        df(~keep_idx, :) = [];
    
        % Add CS size
        [idx1, idx2] = ismember(df.CS, cs_size.CS);
        df.CS_size(:) = nan;
        df.CS_size(idx1) = cs_size.GroupCount(idx2(idx1));

    else
        % remove negative CS (keep lead variant at this region)
        % max_pip_idx = df.PIP == max(df.PIP);
        keep_idx = df.PIP >= opts.PIPmax | df.SNP == index.SNP(k) | max_pip_idx;
        df(~keep_idx, :) = [];
        df.CS(:) = nan; df.CS_size(:) = nan;
    end

    df.Index(:) = index.SNP(k);
    df.Locus(:) = index.locus(k);
    if opts.method == "finemap-inf"
        df = sortrows(df, "PIP", "descend");
    else
        df = sortrows(df, ["CS", "PIP"], ["ascend", "descend"]);
    end

    % double check index SNPs for multiallelic variants and exclude those
    % not GW significant (based on 'pthresh') and throw a warning. This may
    % happen if index variant at a specific locus is duplicated
    % (multiallelic)
    lead_variants = index.SNP(k).split(",");
    lead_idx = ismember(df.SNP, lead_variants) & df.P > opts.pthresh & ...
        ~(df.CS > opts.csfilter & df.PIP >= opts.PIPmax);
    if any(lead_idx)
        fprintf("\t[WARNING-multiallelic index] %d variant was removed (P > %.2E)\n", nnz(lead_idx), opts.pthresh)
        df(lead_idx, :) = [];
    end

    % log10 Bayes Factor filtering
    bf_filter = df.CS_LOG10BF <= opts.log10BF & df.CS > opts.csfilter;
    if any(bf_filter)
        tmp = df(bf_filter, :);
        n_cs_out = numel(unique(tmp.CS));
        n_cs_all = numel(unique(df.CS(df.CS > opts.csfilter)));
        df(bf_filter, :) = [];
        fprintf("\t[LBF filtering]:\n" + ...
            "\t\tLocus:%s\n" + ...
            "\t\tremoved %d (out of %d) CS with log10 BF <= %.1f\n", ...
            tmp.Locus(1), n_cs_out, n_cs_all, opts.log10BF);
    end

    % GW-P filtering (remove CSs that don't contain any P > 'pthresh'
    cs_p = groupsummary(df, "CS", @min, "P");
    cs_p(cs_p.CS <= opts.csfilter | cs_p.fun1_P <= opts.pthresh, :) = [];
    if ~isempty(cs_p)
        p_filter = ismember(df.CS, cs_p.CS);
        tmp = df(p_filter, :);
        df(p_filter, :) = [];
        n_cs_all = numel(unique(df.CS(df.CS > opts.csfilter)));
        fprintf("\t[P filtering]:\n" + ...
            "\t\tLocus:%s\n" + ...
            "\t\tremoved %d (out of %d) CS without any GW P <= %.2E\n", ...
            tmp.Locus(1), height(cs_p), n_cs_all, opts.pthresh);
    end

    res{k} = df;
    clear df
end

res(cellfun(@isempty, res)) = [];
res = vertcat(res{:});

if opts.flames || opts.caldera
    writeFlamesCSfiles(res, opts);
    return
end

% annotate variants
t1 = tic;
fprintf("\t[INFO]: annotating variants...")
vep = mapVariant2gene(res, refGenome=opts.refGenome, verbose=false);
t1 = seconds(toc(t1));
fprintf('\b\b done in %.2f min.\n', minutes(t1))

[~, idx] = ismember(res.SNP, vep.id);
res.Consequence = vep.consequence(idx);
res.variant = vep.snp(idx);
res.nearestGene = vep.nearestGene(idx);
res = movevars(res, ["Index", "Locus"], Before=1);
res = movevars(res, ["variant", "Consequence", "nearestGene"], After="SNP");
res = renamevars(res, "BETA_MEAN", "Beta");
if any(colnames(res) == "BETA_SD")
    res = renamevars(res, "BETA_SD", "SD");
else
    res.SD(:) = nan;
end

if ~isempty(opts.delete)
    arrayfun(@delete, opts.delete);
end
method = opts.method;

end % END

%% ------------------------------------------------------------------------
function [out, index] = mergeRegions(index, win)
% Merge overlapping [pos1,pos2] intervals per chromosome in table `index`
% Overlap rule: intervals merge if next.pos1 <= current running max(pos2)

index_cols = colnames(index);
chr_col = ismember(lower(index_cols), ["chr", "chrom"]);
chr_col = index_cols(chr_col);

if nargin > 1
    pos_col = find(ismember(lower(index_cols), ["bp", "pos", "genpos"]));
    index.pos1 = max(1, index.(pos_col) - win);
    index.pos2 = index.(pos_col) + win;
end

index.idx = (1:height(index))'; % to keep track of merged loci
T = index(:, [chr_col ,"pos1", "pos2", "idx"]);
T.(chr_col) = string(T.(chr_col));                  % normalize
T = sortrows(T, [chr_col ,"pos1", "pos2"]); % single sort (O(n log n))

[chrU, ~, gi] = unique(T.(chr_col), 'stable');

parts = cell(numel(chrU),1);
for g = 1:numel(chrU)
    m = gi == g;
    s = T.pos1(m);
    e = T.pos2(m);
    ids = T.idx(m); % provenance indices

    % linear sweep within this chromosome (O(k))
    runMax = cummax(e);
    % new block starts when the next start is AFTER previous max end
    sep = [true; s(2:end) > runMax(1:end-1)];
    blockId = cumsum(sep);

    mergedStart = s(sep);
    mergedEnd   = accumarray(blockId, e, [], @max);
    mergedIdx   = accumarray(blockId, ids, [], @(v){v(:).'});

    parts{g} = table(repmat(chrU(g), numel(mergedStart),1), ...
                     mergedStart, mergedEnd, mergedIdx, ...
                     'VariableNames', [chr_col ,"pos1", "pos2", "idx"]);
end

out = vertcat(parts{:});

end % END

%% ------------------------------------------------------------------------
function T = mergeFinemapTabs(df, tags)

assert(iscell(df) && ~isempty(df), 'df must be a non-empty cell array of tables');
tags = string(tags);
tags = tags.erase(textBoundary("start") + ".");
dot_idx = tags ~= "";
tags(dot_idx) = "." + tags(dot_idx);
assert(numel(tags)==numel(df), 'tags must match df in length');
if iscolumn(tags), tags = tags'; end

baseVars = ["Index","Locus","CHR","POS","SNP","variant","Consequence","nearestGene","A1","A2","P","key"];
measVars = ["PIP","Beta","CS","CS_size","SD","CS_LOG10BF","lead_r2"];

% --- normalize and ensure key exists in each table ---

% create a dict for union of columns and their data type
dt = dictionary(string.empty, string.empty);
for i = 1:numel(df)
    T = df{i};
    cols = colnames(T);  % your custom function
    for j = 1:numel(cols)
        col = cols(j);
        if ~isKey(dt, col)
            dt(col) = string(class(T.(col)));
        end
    end
end

all_cols = cellfun(@colnames, df, uni=false);
all_cols = unique(horzcat(all_cols{:}));
for i = 1:numel(df)
    Ti = df{i};
    v = colnames(Ti);
    % Build key if missing
    if ~any(v=="key")
        Ti.key = Ti.SNP + ":" + Ti.A1 + ":" + Ti.A2;
    else
        Ti.key = string(Ti.key);
    end

    % missing cols
    miss_cols = setdiff(all_cols, v);
    if ~isempty(miss_cols)
        for k = 1:numel(miss_cols)
            Ti.(miss_cols(k))(:) = "";
            Ti = convertvars(Ti, miss_cols(k), dt(miss_cols(k)));
        end
    end
    df{i} = Ti;
end

% --- base table: take first occurrence per key across all tables ---
baseAll = vertcat(df{:});
% keep only the baseVars that actually exist (defensive)
keepBase = baseVars(ismember(baseVars, colnames(baseAll)));
baseAll = baseAll(:, keepBase);
baseAll.key = string(baseAll.key);

% first non-missing record per key (stable)
[~, ia] = unique(baseAll.key, 'stable');
base = baseAll(ia, :);

% --- start result with base ---
T = base;

% --- add tagged measurement columns from each table via full outer join ---
for i = 1:numel(df)
    Ti = df{i};
    avail = measVars(ismember(measVars, colnames(Ti)));
    if isempty(avail), continue; end

    add = Ti(:, ["key", avail]);
    % rename e.g. PIP -> PIP.tag
    add = renamevars(add, avail, avail + tags(i));

    % full outer join on key; keep one merged key column
    T = outerjoin(T, add, Keys="key", MergeKeys=true, Type="full");
end

% --- final column order exactly as requested ---
want = [baseVars, ...
        compose("PIP%s", tags), ...
        compose("Beta%s", tags), ...
        compose("CS%s", tags), ...
        compose("CS_size%s", tags), ...
        compose("SD%s", tags),...
        compose("CS_LOG10BF%s", tags), ...
        compose("lead_r2%s", tags)];
% keep only those that exist (in case some inputs lacked a var)
want = want(ismember(want, colnames(T)));
T = movevars(T, setdiff(want, colnames(T),'stable')); % (no-op safeguard)
T = T(:, want);
T.key = [];

T = sortrows(T, "POS");
dt = underlyingType(T.CHR);
T.CHR = string(T.CHR);
[~, idx] = natsort(T.CHR);
T = T(idx, :);
T = convertvars(T, "CHR", dt);

end % END

%% ------------------------------------------------------------------------
function T = susieinf_cs_log10bf(T)
% SUSIEINF_CS_LOG10BF  Compute per-variant CS log10 Bayes factors (BF)
%   from a SuSiE-inf flat table, matching SuSiE-R’s internal logic.
%
%   This function reconstructs component-level log10 BFs for each
%   single-effect component (ℓ) and assigns that value to all variants
%   in the corresponding credible set (CS).
%
%   THEORY:
%     - Variant-level lnBF (SuSiE-R: res$lbf → fit$lbf_variable[l,],
%     SuSiE-inf: lbf_variable{ℓ}) gives evidence for SNP j being the causal
%     SNP *if effect ℓ exists*.
%
%     - Component-level lnBF (SuSiE-R: res$lbf_model, fit$lbf[ℓ])
%       aggregates across all SNPs in that effect using priors:
%
%           lbf_component(ℓ) = logsumexp_j( lbf_variable_{jℓ} + log π0_j )
%
%       This is the evidence for *including effect ℓ at all*,
%       i.e. model with CS ℓ vs. model without CS ℓ.
%
%     - Each CS is constructed from one effect’s PIPs (alpha{ℓ}),
%       so its log10BF is simply that effect’s component log10BF.
%
%   INPUT:
%     T : MATLAB table with columns
%         - alpha1..alphaL         : per-effect per-variant PIPs (α_{jℓ})
%         - lbf_variable1..lbf_variableL : per-variant lnBF (natural log)
%         - cs                     : credible set ID (-1 if variant not in any CS)
%
%   OUTPUT:
%     CS_LOG10BF : p×1 vector (p = #variants)
%       - For variants in a CS: log10 Bayes factor of that CS
%         (same for all variants in that CS).
%       - For variants with cs = -1: NaN
%
%   NOTES:
%     - Variant log10BF vs CS log10BF:
%
%       Variant log10BF (per-SNP):
%           log10BF_{jℓ} = lbf_variable_{jℓ} / ln(10)
%           → evidence SNP j is causal for effect ℓ, IF that effect exists.
%
%       Component/CS log10BF (per-effect):
%           log10BF_ℓ = (1/ln(10)) * logsumexp_j(lbf_variable_{jℓ} + log π0_j)
%           → evidence that effect ℓ (and its CS) should be included at all.
%
%       PIP_{jℓ} = softmax_j( lbf_variable_{jℓ} + log π0_j )
%           → posterior probability SNP j is the causal SNP for effect ℓ.
%
%     - FinnGen’s definition:
%           cs_log10bf = log10 BF of comparing
%                        model with this CS included
%                        vs model with this CS removed.
%       This is exactly the component log10BF assigned to that CS.


vn = colnames(T);
if any(vn == "CS_LOG10BF"), return; end % susieR

% Grab exactly alpha1..alphaL (exclude scalar 'alpha') and lbf_variable* in numeric order
alpha_cols = find(vn.startsWith("alpha") & vn ~= "alpha");
if isempty(alpha_cols), return; end
[~, ia] = sort(double(erase(vn(alpha_cols), "alpha")));
alpha_cols = alpha_cols(ia);

lbfv_cols = find(vn.startsWith("lbf_variable"));
[~, ib] = sort(double(erase(vn(lbfv_cols),"lbf_variable")));
lbfv_cols = lbfv_cols(ib);

Alpha  = T{:, alpha_cols};        % p x L  (per-effect PIPs)
LBFVAR = T{:, lbfv_cols};         % p x L  (per-variant ln-BF, natural log)
cs     = T.cs;                    % p x 1  (CS id; -1 if none)

[p,L] = size(Alpha);
assert(size(LBFVAR,2)==L, 'alpha and lbf_variable column counts disagree.');

% Component ln-BF via SuSiE aggregation:
%   lbf_component(ℓ) = logsumexp_j( LBFVAR(j,ℓ) + log π0_j )
% Uniform prior if not provided: log π0_j = -log(p)
X = LBFVAR - log(p);              % add log prior
m = max(X,[],1);                  % 1 x L (stability)
lbf_comp = m + log(sum(exp(X - m),1));   % 1 x L (natural log)
log10bf_comp = lbf_comp ./ log(10);      % convert to base-10

% Map each CS to its generating effect: the one with max α-mass inside the CS
T.CS_LOG10BF = nan(p,1);
ucs = unique(cs); ucs = ucs(ucs >= 0);
for k = ucs.'
    in_cs = (cs == k);
    [~, ell] = max(sum(Alpha(in_cs, :), 1));
    T.CS_LOG10BF(in_cs) = log10bf_comp(ell);
end

% fprintf("\t[INFO]: computed log10 BF for each single effect in finemap-inf.\n")

end % END

%% ------------------------------------------------------------------------
function writeFlamesCSfiles(res, opts)

res(res.CS < 0, :) = [];
assert(~any(ismissing(res.CS)), "Missing CS!")

% sanity check: CS size should be equal to number of variants picked (i.e.
% no filtering was done)
res.iid = res.Index + ":" + res.CS;
gs1 = groupsummary(res, "iid");
gs2 = res(:, ["iid", "CS_size"]);
gs2 = unique(gs2, "rows");
[~, idx] = ismember(gs1.iid, gs2.iid); gs2 = gs2(idx, :);
assert(all(gs1.GroupCount == gs2.CS_size), "Something went wrong (extra filtering happened?)")
clear gs1 gs2

% get CS:locus pairs
res.iid = res.Locus + ":" + res.CS;
loc = unique(res.iid);
t1 = tic;

msg = ifelse(opts.flames, "FLAMES", "CALENDRA");
fprintf("\t[INFO]: writing %d %s CS files.\n", numel(loc), msg)


% index table for FLAMES
if opts.flames
    index = table();
    [index.Filename, index.Annotfiles] = deal(strings(numel(loc), 1));
    opts.annot_path = fullfile(fileparts(opts.cs_path), "annot");
    if ~isfolder(opts.annot_path)
        mkdir(opts.annot_path);
    end
end

alpha_cols = colnames(res);
alpha_cols = alpha_cols(alpha_cols.contains("alphax_" + digitsPattern));

if opts.caldera
    out = cell(numel(loc), 1);
end

for k = 1:numel(loc)
    df = res(res.iid == loc(k), :);
    name1 = df.Locus(1); name2 = "_cs" + df.CS(1);
    name1 = unique(name1.split(","));
    name1 = name1.join("_");

    % use alphaX instead of PIP (see the notes on using prob per CS instead
    % of overall PIP)
    [~, idx] = max(sum(df{:, alpha_cols}, 1));
    alpha_col = alpha_cols(idx);
    
    if opts.flames
        index.Filename(k) = makeWSLpath(fullfile(opts.cs_path, name1 + name2 + ".txt"));
        index.Annotfiles(k) = makeWSLpath(fullfile(opts.annot_path, name1 + name2 + ".txt"));
    
        name = fullfile(opts.cs_path, name1 + name2 + ".txt");
        if isfile(name), continue; end

        df.index = df.CHR + ":" + df.POS + ":" + df.A1 + ":" + df.A2;
        df = df(:, ["index", alpha_col]);
        df = renamevars(df, alpha_col, "PIP");
        assert(isempty(duplicates(df.index)), "variants in CS must be unique!")
        df = renamevars(df, ["index", "PIP"], ["cred1", "prob1"]);
        sum_alpha = sum(df.prob1);

    else
        df = df(:, ["CHR", "POS", "SNP", alpha_col]);
        df.locus(:) = name1 + name2;
        df = renamevars(df, ["CHR", "POS", alpha_col], ...
            ["chr", "bp", "pip"]);
        sum_alpha = sum(df.pip);
    end

    
    if sum_alpha > 1
        fprintf("\t\t[WARNING]: cumulative prob (%.2f) > 1 for %s\n", sum_alpha, loc(k))
    else
        fprintf("\t\t[CUM PROB]: %s: %.2f\n", loc(k), sum_alpha)
    end
    
    if opts.flames
        writetable(df, name, Delimiter="\t", WriteVariableNames=true)
        dos2unix(name); clear name
    else
        out{k} = df;
    end

    clear df name1 name2
end

if opts.caldera
    out = vertcat(out{:});
    ifile = fullfile(opts.cs_path, "cs_file.txt");
    writetable(out, ifile, Delimiter="\t")
    dos2unix(ifile)
end

t1 = seconds(toc(t1));
fprintf('\t\tDone in %.2f min.\n', minutes(t1))

if opts.caldera, return; end

% write index file to upper dir
ifile = fullfile(fileparts(opts.cs_path), "indexfile.txt");
writetable(index, ifile, Delimiter="\t")
dos2unix(ifile)

end % END
