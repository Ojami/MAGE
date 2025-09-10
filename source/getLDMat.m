function index = getLDMat(index, opts)
% getLDMat constructs LD matrices (npz/gz) files for fine-mapping purposes.
% Oveis Jamialahmadi, University of Gothenburg, 17AUG2025

arguments
    index {mustBeA(index, "table")} % a table with required columns of chr, pos to define regions (e.g. from clumping or conditional analyses)
    opts.win (1,1) double = 1e6 % window around the lead variants in 'index' in bp (default is 1 mb).
    opts.path {mustBeTextScalar} = "LDcache" % output dir to write LD files
    opts.eid double {mustBeVector} % vector of eids for which LD matrices are to be generated (default: unrelated Brits from UKBB)
    opts.genodir {mustBeFolder} % path to bgen files 
    opts.parallel (1,1) logical = true
    opts.threads (1,1) double = 30 % only if 'parallel' is true
    opts.exclude_mhc (1,1) logical = true % remove MHC region on chr 6
    opts.method {mustBeMember(opts.method, ["susie-inf", "polyfun", "R"])} = "R" % target fine-mapping tool/env (R: susie R package)
    opts.sumstat {mustBeFile} % only if method is "R"/"susie-inf"
    opts.tril (1,1) logical = true % full:false, triangular:true. Only applicable if method is susie-inf
    opts.dtype {mustBeMember(opts.dtype, ["single", "double"])} = "double" % data type used for getInsampleLD
    opts.covars {mustBeA(opts.covars, "table")} = table() % optional covariate table with "eid" column to be regressed out from LD matrix (applicable only opts.metho == "R")
    opts.merge (1,1) logical = false % to merge overlapping regions (similar to FINNGEN pipeline)
end

% assert index has chr,pos cols
assert(all(ismember(["chr", "pos"], colnames(index))),...
    "chr and pos columns must be present in index table!")
index = index(:, ["chr", "pos"]);

% remove MHC region
if opts.exclude_mhc
    % PolyFun based (covers both GRCh37/38)
    idx = index.chr == "6" & index.pos >= 28000000 & index.pos <= 34000000;
    if any(idx)
        fprintf("[MHC region filter]: %d region(s) were excluded\n", nnz(idx))
    end
    index(idx, :) = [];

end

% construct regions
index.pos1 = max(1, index.pos - opts.win);
index.pos2 = index.pos + opts.win;
index.chr = regexprep(string(index.chr), "^chr", "");
if opts.merge
    index = mergeRegions(index);
end
index.region = "chr" + index.chr + "_" + index.pos1 + "_" + index.pos2;

% set up parallel pool
if opts.parallel && isempty(gcp("nocreate"))
    parpool("Processes", opts.threads);
end

% create output dir
if ~isfolder(opts.path), mkdir(opts.path); end

if ~isfield(opts, "eid")
    opts.eid = getQCEID(3, false);
end

% print info
disp("------------- getLDMat arguments -------------")
fprintf("\tWindow around each locus: +/- %.2f mb\n", opts.win/1e6)
fprintf("\tSample size for LD calculation: %d\n", numel(opts.eid))
fprintf("\tMethod: %s\n", opts.method)
disp("----------------------------------------------")

if opts.method ~= "polyfun"
    % REGENIE support only for now
    opts.sumstat = tabularTextDatastore(opts.sumstat, ...
        TextType="string", ...
        VariableNamingRule="preserve");

    if any(ismember(opts.sumstat.SelectedVariableNames, "LOG10P"))
        sumstat_mlog10 = true;
        p_col = "LOG10P";
    else
        sumstat_mlog10 = false;
        p_col = "P";
    end

    opts.sumstat.SelectedVariableNames = ["CHROM", "GENPOS", "ID", "ALLELE1", ...
        "ALLELE0", "BETA", "SE", p_col, "N"];
    opts.sumstat = tall(opts.sumstat); 
    if sumstat_mlog10
        opts.sumstat.P = 10.^-opts.sumstat.LOG10P;
        opts.sumstat.LOG10P = [];
    end

end

if opts.method ~= "R"
    opts.spar = py.importlib.import_module("scipy.sparse");
    opts.np = py.importlib.import_module("numpy");
end

for k = 1:height(index)
    
    t1 = tic;

    tag = index.region(k);
    fprintf('[%d of %d] calculating LD mat for region %s\n', k, height(index), tag)

    % check if region already exists
    if opts.method == "R"
        targetFile = fullfile(opts.path, tag + ".bin");
        is_sumstat_chunk = isfile(fullfile(opts.path, tag + ".txt"));
    else
        targetFile = fullfile(opts.path, tag + ".npz");
        is_sumstat_chunk = true;
    end

    if isfile(targetFile) && is_sumstat_chunk
        fprintf(" (already exists)")
    else
        snps = struct(snp=index.pos1(k) + "-" + index.pos2(k), chr=index.chr(k));
        bgi = bgireader(snps, ...
            bgenhome=opts.genodir, ...
            parallel=false, ...
            verbose=false);
        bgi = vertcat(bgi{:});

        if opts.method ~= "polyfun"

            exportBinFiles(bgi, opts.sumstat, index(k, :), ...
                fullfile(opts.path, tag), ...
                opts);
        else

            variants = unique([bgi.snp]');
            chrs = repmat(string(bgi(1).chr), numel(variants), 1);
            ld = getInsampleLD(variants, ...
                chrs, ...
                parallel=true, ...
                bgenhome=opts.genodir, ...
                ldsample=opts.eid, ...
                a2freq=true, ...
                maf=1e-10, ...
                dtype=opts.dtype, ...
                removenan=false); 
        
            rstab = table(ld.snp, ld.chr, ld.pos, ld.a2, ld.a1, ...
                VariableNames=["rsid", "chromosome", "position", "allele1", "allele2"]); % <-- A1 is effect allele, however in getInsampleLD A2 is effect allele
            writetable(rstab, fullfile(opts.path, tag + ".txt"), ...
                Delimiter="\t", FileType="text");
            movefile(fullfile(opts.path, tag + ".txt"), ...
                fullfile(opts.path, tag))
            dos2unix(fullfile(opts.path, tag))
            gzip(fullfile(opts.path, tag))
            delete(fullfile(opts.path, tag))
            
            R = pyrun("R = np.tril(ld_arr).astype(np.float64)", "R", np=opts.np, ld_arr=ld.ld);
            pyrun("np.fill_diagonal(R, np.diag(R)/2.0)", R=R, np=opts.np)
            pyrun("R = sparse.coo_matrix(R)", "R", R=R, sparse=opts.spar);
            pyrun("sparse.save_npz(npz_file, R, compressed=True)", ...
                sparse=opts.spar,npz_file=fullfile(opts.path, tag))

        end

    end
    
    t1 = seconds(toc(t1));
    fprintf('[done in %.2f min]\n\n', minutes(t1))
end

% for debugging
% R = pyrun("R = spar.load_npz(npz).toarray()", "R", spar = spar, npz = ldregion.files(i));

% % these files are  not available on borad institute servers anymore
% links = fullfile(links(1), fetchfiles);
% for i = 1:numel(fetchfiles)
%     websave(fullfile(opts.cachedir, fetchfiles(i)), links(i));
% end

end % END

%% subfunctions ===========================================================
function exportBinFiles(bgi, gss, index, targetFile, opts)

if opts.method == "R"
    is_ldfile = isfile(targetFile + ".bin");
else
    is_ldfile = isfile(targetFile + ".npz");
end


% subset sumstat for this region only
dtype = tall.getClass(gss.CHROM);
index = convertvars(index, "chr", str2func(dtype));

matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
chunk_idx = gss.CHROM == index.chr & gss.GENPOS >= index.pos1 & gss.GENPOS <= index.pos2;
df = gather(gss(chunk_idx, :));
matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);

% remove multiallelic variants (indels)
alMat = sort(df{:, ["ALLELE0", "ALLELE1"]}, 2);
df.tmp = df.CHROM + ":" + df.GENPOS + ":" + alMat(:, 1) + ":" + alMat(:, 2);
dup_vars = duplicates(df.tmp);
if ~isempty(dup_vars)
    dup_idx = ismember(df.tmp, dup_vars);
    fprintf("\t\t[multiallelic indels]: removed %d variants\n", nnz(dup_idx))
    df(dup_idx, :) = [];
end

df.tmp = [];

% remove variants not present in sumstat file
bgi = struct2table(bgi);
bgi(~ismember(bgi.snp, df.ID), :) = [];

% check if 'covariates' should be regressed out
if isempty(opts.covars)
    getCorr = true;
else
    getCorr = false;
end

if ~is_ldfile
    ld = getInsampleLD(bgi.snp, ...
        bgi.chr, ...
        parallel=true, ...
        bgenhome=opts.genodir, ...
        ldsample=opts.eid, ...
        a2freq=true, ...
        maf=1e-10, ...
        dtype=opts.dtype, ...
        removenan=false, ...
        matchid=df.CHROM + ":" + df.GENPOS + ":" + df.ALLELE0 + ":" + df.ALLELE1, ...
        getCorr=getCorr); % <-- A2 is treated as effect allele in getInsampleLD 
    
    % regress out covariates per SuSiE's recommendations (Also:
    % https://www.nature.com/articles/s41586-023-06592-6)
    if ~isempty(opts.covars)
        % match eid
        [idx1, idx2] = ismember(ld.eid, opts.covars.eid);
        
        if any(~idx1)
            ld.eid(~idx1) = [];
            ld.bed(~idx1, :) = [];
        end
    
        opts.covars = opts.covars(idx2(idx1), :);
        opts.covars.eid = [];
        ld = rmfield(ld, "eid");
    
        fprintf("\t\tRegressing out (%d,%d) covariate matrix\n", ...
            height(opts.covars), width(opts.covars))
    
        ld.bed = residualize_genotypes(ld.bed, opts.covars{:, :});
        ld.ld = corr(ld.bed);
        ld = rmfield(ld, "bed");
    end

    % final sanity check
    assert(all(df.ID == ld.snp & df.ALLELE1 == ld.a2 & df.ALLELE0 == ld.a1), ...
        "failed to match variants in sumstat to LD block!")

    % write chunk files to output
    ld = ld.ld;

end


% bin/shape files for susieR
if opts.method == "R"

    if ~is_ldfile
        m = size(ld,1); n = size(ld,2);
        
        fid = fopen(targetFile + ".bin",'w','ieee-le');     % little-endian doubles
        fwrite(fid, ld, 'double'); 
        fclose(fid);
        
        % shape metadata (space-delimited)
        writematrix([m n], targetFile + ".shape", ...
            Delimiter=" ", ...
            FileType="text",...
            QuoteStrings="none");
    end
    
    df = renamevars(df, ...
        ["CHROM", "GENPOS", "ID", "ALLELE1", "ALLELE0"],...
        ["CHR","POS", "SNP", "A1","A2"]);
    writetable(df, targetFile + ".txt", Delimiter="\t")

elseif opts.method == "susie-inf"
    
    if ~is_ldfile
        if opts.tril
            R = pyrun("R = np.tril(ld_arr).astype(np.float64)", "R", np=opts.np, ld_arr=ld);
        else
            R = pyrun("R = np.array(R, dtype=np.float64)", "R", np=opts.np, R=ld);
            pyrun("np.fill_diagonal(R, 1.0)", np=opts.np, R=R)
        end
        pyrun("np.fill_diagonal(R, 1.0)", np=opts.np, R=R);
        pyrun("np.savez_compressed(out, R)", np=opts.np, out=targetFile, R=R);
    end

    old_cols = ["CHROM", "GENPOS", "ID", "ALLELE1", "ALLELE0", "BETA", "SE", "P"];
    new_cols = ["chromosome", "position", "rsid", "allele1", "allele2", "beta", "se", "p"];
    n_eff = df.N(1);
    df.N = [];
    df = renamevars(df, old_cols, new_cols);
    df = df(:, new_cols);
    writetable(df, targetFile + ".txt", Delimiter="\t")
    save(targetFile + ".n_eff.mat", "n_eff")
    dos2unix(targetFile + ".txt");
end

fprintf("\t\tRegion has %d variants\n", height(df))

end % END

%% ------------------------------------------------------------------------
function X = residualize_genotypes(X, Z)
% X        : n x p genotype matrix (double; same samples as GWAS)
% Z    : n x q covariates used in GWAS (NO intercept column)


% Input checks
if ~isfloat(X)
    X = double(X);
end

if ~isfloat(Z)
    Z = double(Z);
end
[n, ~] = size(X);

if size(Z,1) ~= n, error('Zbase must have the same #rows as X.'); end

Z = [ones(n,1,'like',X), Z];

% Drop any NaN columns (just in case)
bad = any(~isfinite(Z),1);
if any(bad)
    Z(:,bad) = [];
end

% QR with column pivoting for numerical stability
[Q, Rz] = qr(Z, 0); 
r = rank(Rz);
Q = Q(:,1:r);  % orthonormal basis for col(Z)

% Residualize: (I - Q*Q') * X
X = X - Q * (Q' * X);

end % END

%% ------------------------------------------------------------------------
function out = mergeRegions(index)
% Merge overlapping [pos1,pos2] intervals per chromosome in table `index`
% Overlap rule: intervals merge if next.pos1 <= current running max(pos2)

T = index(:, {'chr','pos1','pos2'});
T.chr = string(T.chr);                  % normalize
T = sortrows(T, {'chr','pos1','pos2'}); % single sort (O(n log n))

[chrU, ~, gi] = unique(T.chr, 'stable');

parts = cell(numel(chrU),1);
for g = 1:numel(chrU)
    m = gi == g;
    s = T.pos1(m);
    e = T.pos2(m);

    % linear sweep within this chromosome (O(k))
    runMax = cummax(e);
    % new block starts when the next start is AFTER previous max end
    sep = [true; s(2:end) > runMax(1:end-1)];
    blockId = cumsum(sep);

    mergedStart = s(sep);
    mergedEnd   = accumarray(blockId, e, [], @max);

    parts{g} = table(repmat(chrU(g), numel(mergedStart),1), ...
                     mergedStart, mergedEnd, ...
                     'VariableNames', {'chr','pos1','pos2'});
end

out = vertcat(parts{:});

end % END
