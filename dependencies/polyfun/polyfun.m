function polyfun(file, opts)
% A wrapper for polyfun package (limited)
% implements the workflow represented here:
% https://github.com/omerwe/polyfun/wiki
% 
% Note: since we take care of duplicated indels in the input sumstat file,
% and the fact that we use --allow-swapped-indel-alleles in finemapper to
% swap indels (unique indels but flipped A1/A2, this is OK when ldscore and
% ld matrix files are in-sample), we manually set allow_duplicates to True
% in polyfun_utils.set_snpid_index() function to avoid duplicate errors in
% LDmatrix file (happenes with --allow-swapped-indel-alleles). We no longer
% care for duplicated variants in LDmatrix file, because our sumstat
% doesn't have such duplicated [indel] variants.
% 
% Note: as mentioned here: https://github.com/omerwe/polyfun/pull/149,
% rpy2 throws an error and can be fixed by changing this line in finemapper
% function:
% self.susie_dict = {key:np.array(susie_obj.rx2(key), dtype=object) for key in list(susie_obj.names)}
% to:
% self.susie_dict = {key:np.array(susie_obj.rx2(key), dtype=object) for key in (susie_obj.names).tolist()}

arguments
    file {mustBeFile}
    opts.conda {mustBeTextScalar} = "polyfun" % conda env name
    opts.prior (1,1) {mustBeMember(opts.prior, ["precomputed", "S-LDSC", "non-param"])} = "S-LDSC"
    opts.polyfunwd {mustBeFolder} = fullfile(fileparts(which("polyfun.m")), 'polyfun/')
    opts.annwd {mustBeFolder} = fullfile(fileparts(which("polyfun.m")), 'baselineLF2.2.UKB/')
    opts.method {mustBeMember(opts.method, ["susie", "finemap", "both"])} = "susie"
    opts.max_num_causal (1,1) double {mustBeGreaterThan(opts.max_num_causal, 0)} = 10
    opts.pvalue_cutoff (1,1) double {mustBeInRange(opts.pvalue_cutoff, 0, 1)} = 1e-7 % for fine-mapping
    opts.finemap_exe {mustBeFile} = fullfile(fileparts(which("polyfun.m")), 'finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64')
    opts.ldstore2 {mustBeFile} = fullfile(fileparts(which("polyfun.m")), 'ldstore_v2.0_x86_64/ldstore_v2.0_x86_64')
    opts.incl_samples {mustBeVector, mustBeNonmissing} % included sample files for LD matrix calculations (for ldstore2)
    opts.cachedir {mustBeTextScalar} % LD cache folder % = fullfile(fileparts(which("polyfun.m")), 'LDcache') % for LD matrices
    opts.genodir {mustBeFolder} % dir to bgen/bed files for calculating LD matrix
    opts.workers (1,1) double = 40
    opts.output {mustBeTextScalar}

    % non-native arguments
    opts.calcLDmat (1,1) logical = true % calculate LDmat in MATLAB and store them as npz files, 'geno' must be used.
    opts.index {mustBeA(opts.index, 'table')} % use index variants in this table for defining the regions. By default, it uses create_finemapper_jobs to define fine-mapping regions
    opts.win (1,1) double = 1.5e6 % window around the lead variants in 'index' in bp (default is 1.5 mb).
    opts.removeCachedir (1,1) logical = false % to remove the folder 'cachedir'
end

% apply changes to polyfun_utils and finemapper (see above)
pt = readlines(fullfile(opts.polyfunwd, "polyfun_utils.py"));
idx = pt.startsWith("def set_snpid_index(df, copy=False, allow_duplicates=False");
if any(idx)
    pt(idx) = "def set_snpid_index(df, copy=False, allow_duplicates=True, allow_swapped_indel_alleles=False):";
    writelines(pt, fullfile(opts.polyfunwd, "polyfun_utils.py"))
    dos2unix(fullfile(opts.polyfunwd, "polyfun_utils.py"))
end

pt = readlines(fullfile(opts.polyfunwd, "finemapper.py"));
idx = pt.startsWith(textBoundary("start") + whitespacePattern + "self.susie_dict = {key:np.array(susie_obj.rx2(key), dtype=object) for key in list(susie_obj.names)}");
if any(idx)
    pt(idx) = replace(pt(idx), "list(susie_obj.names)", "(susie_obj.names).tolist()");
    writelines(pt, fullfile(opts.polyfunwd, "finemapper.py"))
    dos2unix(fullfile(opts.polyfunwd, "finemapper.py"))
end

% get weights/ld/annotation files (SLDSC)
opts.refld = regexprep(opts.annwd, "(\\|/)$", ''); 
opts.refld = split(opts.refld, "/"|"\");
opts.refld(opts.refld == "") = [];
opts.refld = string(opts.refld(end)) + ".";
% weight files should be started with weights
opts.weight = string({dir(opts.annwd).name}.');
opts.weight = natsort(opts.weight(startsWith(opts.weight, "weight")));
opts.weight = strshare(opts.weight);

% make WSL path
opts.pwd = pwd;
opts.pwdwsl = makeWSLpath(opts.pwd) + "/";
opts.polyfunwd = makeWSLpath(opts.polyfunwd);
opts.annwd = makeWSLpath(opts.annwd);
if ~endsWith(opts.annwd, "/"), opts.annwd = opts.annwd + "/"; end
if ~endsWith(opts.polyfunwd, "/"), opts.polyfunwd = opts.polyfunwd + "/"; end
opts.finemap_exe = makeWSLpath(opts.finemap_exe);
opts.ldstore2  = makeWSLpath(opts.ldstore2 );

% check geno files
if isfield(opts, "genodir")
    opts.geno = getfilenames(opts.genodir, "bgen", 'fullpath', true).bgen;
    if isempty(opts.geno) % bed files are there?
        if ~isfield(opts, 'bed')
            opts.geno = getfilenames(opts.genodir, "bed", 'fullpath', true).bed;
            if isempty(opts.geno)
                error('no bgen/bed files were found in %s', opts.bgenhome)
            end
        end
    end

    % construct file string
    [~, name, ext] = fileparts(opts.geno);
    genofiles = name + ext;
    if any(contains(lower(genofiles), ["_cx_", "_cy_", "_cxy_"])) % sex chromosomes
        chrPos = strshare(genofiles, 'pattern', true, 'patternstr', "(?<=_c)(.*?)(?=_)");
    else
        chrPos = strshare(genofiles, 'pattern', true);
    end
    opts.genopattern = replace(chrPos, "%d", "%s"); 

end

if isempty(gcp('nocreate'))
    parpool("local", opts.workers);
end

% check sample size
n = bfilereader(file, "header", true, "summary", "only");
opts.n = double(n.N(1));

% gzip the file
[pth, name] = fileparts(file);
if pth == "", pth = pwd; end

if isfield(opts, 'output')
    [pth, name, ext] = fileparts(string(opts.output));

    if pth == "", pth = pwd; end

    if ~isfolder(pth)
        mkdir(pth)
    end
    output = fullfile(pth, name + ext);

else
    output = fullfile(pth, name + "_polyfun");
end

% samples to be used for LD matrix calculations (to be used with --geno and
% --ldstore2 arguments)
if isfield(opts, 'incl_samples')
    opts.incl_samples = string(opts.incl_samples);
    writematrix(opts.incl_samples, ...
        fullfile(pth, "polyfun_ldstore_samples.txt"), ...
        QuoteStrings="none", FileType="text");
    opts.incl_samples = fullfile(pth, "polyfun_ldstore_samples.txt");
    dos2unix(opts.incl_samples)
end

if ~isfile(output + ".parquet")
    dt = prunesumstat(file, output, opts);

    % ss2fuma(file, output=output)
    % dt = ".gz";
    
    % munging
    outputw = makeWSLpath(output);
    fprintf('step 0: munging...\n')
    mungcmd = "python3 " + opts.polyfunwd + "munge_polyfun_sumstats.py " + ...
        "--sumstats " + outputw + dt + " --out " + outputw + ...
        ".parquet --min-info 0.001 --min-maf 0.001";
    runbash(mungcmd, "polyfun_bash", "conda", opts.conda, ...
        "parallel", false, "verbose", true)
    
end

% if isfile(output + dt), delete(output + dt); end
% check indels file (those indels with both A1/A2 and A2/A1 present in
% sumstats file, only 1 with the lowest p is kept).
if isfile(output + ".indels.txt")
    opts.rem_indels = readtable(output + ".indels.txt", TextType="string");
end

inname = output + ".parquet";
if isfile(output + ".parquet.log"), delete(output + ".parquet.log"); end
% fprintf('--------------------------------------------------------------\n')

%% 1. Computing prior causal probabilities with PolyFun
opts.name = inname; opts.out = output; 
opts = compPrior(opts);

%% 2. Functionally informed fine mapping with finemapper
if opts.method == "both"
    % perform two fine-mapping, one with SuSiE and one with FINEMAP
    opts.method = "finemap";
    opts = finemapper(opts);
    opts.method = "susie";
    opts = finemapper(opts);
else
    opts = finemapper(opts);
end

%% 3. Estimating functional enrichment using S LDSC
% runSLDSC(opts);

% doesn't work with polyfun:
% %% 4. Estimate heritability (LDSC.files directory must be present)
% runLDSC(opts); 

if isfield(opts, 'incl_samples'), delete(opts.incl_samples); end

if opts.removeCachedir && isfield(opts, 'cachedir')
    rmdir(opts.cachedir, "s");
end

end % END

%% subfunctions ===========================================================
function opts = compPrior(opts)
opts.priordir = opts.out + "." + opts.prior;
if ~isfolder(opts.priordir)
    mkdir(opts.priordir)
end

[pth, name, ext] = fileparts(opts.out);
opts.out = name + ext;

% check if files are already there
gzfiles = getfilenames(opts.priordir, 'gz').gz;
if ~isempty(gzfiles)
    constidx = endsWith(gzfiles, '_constrained.gz');
    if any(constidx)
        opts.snpvar = makeWSLpath(opts.priordir) + "/" + natsort(gzfiles(~constidx));
        opts.snpvarcons = makeWSLpath(opts.priordir) + "/" + natsort(gzfiles(constidx));
        opts.out = fullfile(pth, opts.out);
        return
    end
end

if opts.prior == "precomputed"
    % 1.1.approach 1: Using precomputed prior causal probabilities based on a
    % meta-analysis of 15 UK Biobank traits
    cmd = "python3 " + opts.polyfunwd + "extract_snpvar.py " + ...
        "--sumstats " + makeWSLpath(opts.name) + ...
        " --out " + makeWSLpath(opts.priordir) + ...
        "/" + makeWSLpath(opts.out) + ".snpvar.gz";
    runbash(cmd, "polyfun_bash", "conda", opts.conda, ...
        "parallel", false, "verbose", true)
    
elseif opts.prior == "S-LDSC"
    % -------------------------------------------------------------------------
    % 1.2.PolyFun approach 2: Computing prior causal probabilities via an
    % L2-regularized extension of S-LDSC
    cmd = "python3 " + opts.polyfunwd + "polyfun.py " + ...
        "--sumstats " + makeWSLpath(opts.name) + " --allow-missing --output-prefix " + ...
        makeWSLpath(opts.priordir)  + "/" + makeWSLpath(opts.out) + ...
        " --compute-h2-L2 --no-partitions --ref-ld-chr " + ...
        opts.annwd + opts.refld + " --w-ld-chr " + ...
        opts.annwd + opts.weight;
    runbash(cmd, "polyfun_bash", "conda", opts.conda, ...
        "parallel", false, "verbose", true)

else
    error("this method has not yet been implemented")
    
end

% pick snpvar gz files
gzfiles = getfilenames(opts.priordir, 'gz').gz;
constidx = endsWith(gzfiles, '_constrained.gz');
opts.snpvar = makeWSLpath(opts.priordir) + "/" + natsort(gzfiles(~constidx));
opts.snpvarcons = makeWSLpath(opts.priordir) + "/" + natsort(gzfiles(constidx));

opts.out = fullfile(pth, opts.out);
end 

%% ------------------------------------------------------------------------
function opts = finemapper(opts)
% implements GWAS fine mapping
opts.finemapdir = opts.out + ".finemap." + opts.method;
if ~isfolder(opts.finemapdir)
    mkdir(opts.finemapdir)
end

if ~isfolder(opts.cachedir)
    mkdir(opts.cachedir)
end
opts.cachedirwsl = makeWSLpath(opts.cachedir);
opts.finemapdirwsl = makeWSLpath(opts.finemapdir);
if endsWith(opts.finemapdirwsl, "/")
    opts.finemapdirwsl = regexprep(opts.finemapdirwsl, "/$", "");
end

if isfield(opts, 'index')
    cols = lower(colnames(opts.index));
    opts.index.Properties.VariableNames = cols;
    opts.index.chr = string(opts.index.chr);
end

% step 1: Creating region-specific jobs (3Mb-long regions)
fprintf('finemapper: creating regions...\n')
opts.snpvarcons = natsort(opts.snpvarcons);
codes  = cell(numel(opts.snpvarcons), 1);
for i = 1:numel(opts.snpvarcons)
    filewd = fullfile(opts.finemapdir, "polyfun_all_jobs" + i + ".txt");

    if isfield(opts, 'index')
        % define regions around lead variant in the index table
        [~, thischr] = fileparts(opts.snpvarcons(i)); 
        thischr = erase(regexp(thischr, "[.](\d+)[.]", "match"), ".");
        regions = opts.index(opts.index.chr == thischr, :);
        if isempty(regions)
            continue
        end
        
        regions.start = max(0, regions.bp - opts.win);
        regions.stop = regions.bp + opts.win;
        regions.str = "chr" + regions.chr + "_" + regions.start + "_" + regions.stop;
        regions.out = regions.snp.replace(":","_") + "_" + regions.str;

        cmd = "python3 " + opts.polyfunwd + "/finemapper.py --chr " + ...
            regions.chr + " --start " + regions.start + " --end " + regions.stop + ... 
            " --out " + opts.finemapdirwsl + "/polyfun_all." + ...
            regions.out + ".gz --ld https://UKBB_LD/" + regions.str + ...
            " --method " + opts.method + " --sumstats " + ...
            opts.snpvarcons(i) + " --n " + opts.n + " --finemap-exe " + ...
            opts.finemap_exe + " --memory 1 --max-num-causal " + opts.max_num_causal;
        writelines(cmd, filewd)

    else
        cmd = "python3 " + cc + "create_finemapper_jobs.py " + ...
            " --sumstats " + opts.snpvarcons(i) + ...
            " --n " + opts.n + ...
            " --pvalue-cutoff " + opts.pvalue_cutoff + ...
            " --method " + opts.method + ...
            " --max-num-causal " + opts.max_num_causal + ...
            " --out-prefix " + opts.finemapdirwsl + "/polyfun_all" + ...
            " --finemap-exe " + opts.finemap_exe + ...
            " --jobs-file " + opts.finemapdirwsl + "/polyfun_all_jobs" + i + ".txt";
    
        if ~isfile(filewd)
            runbash(cmd, "polyfun_bash", "conda", opts.conda, ...
                "parallel", false, "verbose", true)
        end
    end

    % add cache option
    file = readlines(filewd); file(file == "") = [];
    if ~isempty(file)

        if isfield(opts, 'geno') && ~opts.calcLDmat % use geno file to compute LD matrices
            file = erase(file, "--ld " + wildcardPattern(1, inf) + whitespaceBoundary);
            file = file + " --cache-dir " +  opts.cachedirwsl;

            % find geno files
            chr = strtrim(extractBetween(file, "--chr", "--"));
            genofiles = compose(opts.genopattern, chr);
            genofiles = arrayfun(@makeWSLpath, fullfile(opts.genodir, genofiles));
            file = file + " --geno " + genofiles;
            
            % if bgen, --sample-file is needed
            bgenidx = endsWith(genofiles, ".bgen");
            if any(bgenidx)
                samplefile = regexprep(genofiles, ".bgen$", ".sample");
                file(bgenidx) = file(bgenidx) + " --sample-file " + ...
                    samplefile(bgenidx);
            end

            % compute LD matrix only for these samples
            if isfield(opts, 'incl_samples')
                file = file + " --incl-samples " + makeWSLpath(opts.incl_samples);
            end

            if isfield(opts, 'ldstore2')
                file = file + " --ldstore2 " + opts.ldstore2;
            end

        else
            % otherwise all LD matrices should be locally available
            % (@11MARCH2022: online LD npz files are unavailable, so
            % fetchLDmat cannot download any new files)
            opts.link = extractBetween(file, "--ld ", whitespacePattern);
            opts = fetchLDmat(opts);
            file = replace(file, opts.link, opts.cachedirwsl + "/" + opts.ldmat);
        end

        codes{i} = file;
        % writematrix(file, filewd, 'QuoteStrings', false, 'FileType', 'text')
    end
    % delete(filewd)
    
end
fprintf('--------------------------------------------------------------\n')

% step 2: Invoking the region-specific jobs
idx = cellfun(@isempty, codes);
codes(idx) = [];
codes = vertcat(codes{:});
codes = codes + " --allow-swapped-indel-alleles";
t1 = tic;
conda = opts.conda;
for k = 1:numel(codes)
    fprintf('running job %d of %d\n', k, numel(codes))
    gzfile = extractBetween(codes(k), "--out ", " ");
    if ~isfile(makeWSLpath(gzfile, true))
         runbash(codes(k), "polyfun_bash", "conda", conda, ...
            "parallel", false, "verbose", true)
    end
    fprintf('job %d is done.\n', k)
    disp("---------------------------------------------------------------")
end
% if mean(contains(codes, "--geno")) >= 0.3 % don't invokde parfor, since LDStore is intensive
%     for k = 1:numel(codes)
%         fprintf('running job %d of %d\n', k, numel(codes))
%          runbash(codes(k), "polyfun_bash", "conda", conda, ...
%             "parallel", false, "verbose", true)
%         fprintf('job %d is done.\n', k)
%     end
% else
%     parfor (k = 1:numel(codes), opts.workers)
%         fprintf('running job %d of %d\n', k, numel(codes))
%          runbash(codes(k), "polyfun_bash", "conda", conda, ...
%             "parallel", false, "verbose", false)
%         fprintf('job %d is done.\n', k)
%     end
% end

% step 3: Aggregating the results
idx = find(~idx); % only aggregate these jobs (the rest is empty)

fprintf('aggregating jobs...')
if ~isfield(opts, 'index')
    for i = 1:numel(idx)
        if isfile(fullfile(opts.finemapdir, "polyfun_agg." + idx(i) + ".txt.gz"))
            continue
        end
        cmd = "python3 " + opts.polyfunwd + "aggregate_finemapper_results.py " + ...
                " --sumstats " + opts.snpvarcons(idx(i)) + ...
                " --pvalue-cutoff " + opts.pvalue_cutoff + ...
                " --out-prefix " + opts.finemapdirwsl + "/polyfun_all" + ...
                " --out " + opts.finemapdirwsl + "/polyfun_agg." + idx(i) + ".txt.gz" + ...
                " --allow-missing-jobs";
        runbash(cmd, "polyfun_bash", "conda", conda, ...
                "parallel", false, "verbose", false)
    end
end

% aggregate into a single table
if ~isfile(fullfile(opts.finemapdir, "polyfun_agg.mat"))
    if isfield(opts, 'index')
        gz = extractBetween(codes, "--out ", " ");
        [~, region, ext] = fileparts(arrayfun(@(x)makeWSLpath(x, true), gz));
        gz = region + ext;
        region = regexprep(region, "^polyfun_all.", "");
        region = split(region, "_chr", 2);
        region = array2table(region, VariableNames=["index", "name"]);
        region.name = replace(region.name, "_", ":");

    else
        gz = getfilenames(opts.finemapdir, "gz").gz;
        gz = gz(startsWith(gz, 'polyfun_agg.'));
    end

    gz = fullfile(opts.finemapdir, gz);
    res = cell(numel(gz), 1);
    for k = 1:numel(gz)
        if ~isfile(gz(k)), continue; end
        res{k} = bfilereader(gz(k), 'header', true, 'verbose', 'off');
        if isfield(opts, 'index')
            res{k}.Region(:) = region.name(k);
            res{k}.Index(:) = region.index(k);
        end
    end
    res(cellfun(@isempty, res)) = [];
    res = vertcat(res{:});
    save(fullfile(opts.finemapdir, "polyfun_agg.mat"), 'res');
    % gz = cellstr(gz);
    % delete(gz{:})
end

t1 = seconds(toc(t1));
fprintf('\b\b Done in %.2f min\n', minutes(t1))

end

%% ------------------------------------------------------------------------
function opts = fetchLDmat(opts)

[~, opts.ldmat] = fileparts(opts.link);
fetchfiles = opts.ldmat + [".gz", ".npz"];
fetchfiles = fetchfiles(:);

doexist = isfile(fullfile(opts.cachedir, fetchfiles));
fetchfiles(doexist) = [];

if isempty(fetchfiles)
    % opts = removeIndelsFromLDmat(opts, fullfile(opts.cachedir, fetchfiles));
    return
end

fetchfiles = unique(regexprep(fetchfiles, [".gz$", ".npz$"], ""));

if isfield(opts, "incl_samples")
    eids = readmatrix(opts.incl_samples);
else
    eids = getQCEID(1, false);
end

region = fetchfiles.split("_", 2);
region(:, 1) = regexprep(region(:, 1), "^chr", "");

terminate(pyenv)
spar = py.importlib.import_module("scipy.sparse");
np = py.importlib.import_module("numpy");

for k = 1:numel(fetchfiles)

    t1 = tic;
    fprintf('(%d of %d) calculating LD mat: %s', k, numel(fetchfiles), fetchfiles(k))

    snps.snp = join(region(k, 2:end), "-");
    snps.chr = region(k, 1);
    bgi = bgireader(snps, bgenhome=opts.genodir, parallel=false, verbose=false);
    bgi = vertcat(bgi{:});
    
    variants = unique([bgi.snp]');
    chrs = repmat(string(bgi(1).chr), numel(variants), 1);
    ld = getInsampleLD(variants, chrs, parallel=true, ...
        bgenhome=opts.genodir, ldsample=eids, a2freq=true, maf=1e-10); 

    rstab = table(ld.snp, ld.chr, ld.pos, ld.a1, ld.a2, ...
        VariableNames=["rsid", "chromosome", "position", "allele1", "allele2"]);
    writetable(rstab, fullfile(opts.cachedir, fetchfiles(k) + ".txt"), ...
        Delimiter="\t", FileType="text");
    movefile(fullfile(opts.cachedir, fetchfiles(k) + ".txt"), ...
        fullfile(opts.cachedir, fetchfiles(k)))
    dos2unix(fullfile(opts.cachedir, fetchfiles(k)))
    gzip(fullfile(opts.cachedir, fetchfiles(k)))
    delete(fullfile(opts.cachedir, fetchfiles(k)))

    R = pyrun("R = np.tril(ld_arr).astype(np.float64)", "R", np = np, ld_arr=ld.ld);
    pyrun("np.fill_diagonal(R, np.diag(R)/2.0)", R=R, np=np)
    pyrun("R = sparse.coo_matrix(R)", "R", R=R, sparse=spar);
    pyrun("sparse.save_npz(npz_file, R, compressed=True)", ...
        sparse=spar,npz_file=fullfile(opts.cachedir, fetchfiles(k)))
    
    t1 = seconds(toc(t1));
    fprintf(' (done in %.2f min)\n', minutes(t1))
end

% for debugging
% R = pyrun("R = spar.load_npz(npz).toarray()", "R", spar = spar, npz = ldregion.files(i));

% % these files are  not available on borad institute servers anymore
% links = fullfile(links(1), fetchfiles);
% for i = 1:numel(fetchfiles)
%     websave(fullfile(opts.cachedir, fetchfiles(i)), links(i));
% end

end

%% ------------------------------------------------------------------------
% function opts = removeIndelsFromLDmat(opts, ldfiles)
% % check if rem_indels are present in LDmatrix files and remove them. This
% % does not affect the LDmatrix files themselves, because we copy them to
% % another backup folder and work with modified files. 
% if ~isfield(opts, 'rem_indels') 
%     return
% end
% 
% 
% end

%% ------------------------------------------------------------------------
% function opts = runSLDSC(opts)
% opts.enrichdir = fullfile(opts.pwd, opts.out + ".enrichment");
% if ~isfolder(opts.enrichdir)
%     mkdir(opts.enrichdir)
% end
% 
% opts.enrichdirwsl = makeWSLpath(opts.enrichdir);
% 
% cmd = "wsl python3 " + opts.polyfunwd + "ldsc.py" + ...
%     " --h2 " + opts.name + ...
%     " --ref-ld-chr " + opts.annwd + opts.refld + ...
%     " --w-ld-chr " + opts.annwd + opts.weight + ...
%     " --out " + opts.enrichdirwsl + "/enrichment" + ...
%     " --overlap-annot --not-M-5-50" + ...
%     " --chisq-max 9999999.0";
% 
% tic
% fprintf('performing functional enrichment with S-LDSC...')
% [~, ~] = system(cmd);
% res = readtable(fullfile(opts.enrichdir, 'enrichment.results'), ...
%     'TextType', 'string', 'FileType', 'text', ...
%     'VariableNamingRule', 'preserve');
% res = sortrows(res, 'Enrichment_p', 'ascend');
% save(fullfile(opts.enrichdir, 'enrichment.results.mat'), 'res')
% delete(fullfile(opts.enrichdir, 'enrichment.results'))
% fprintf('\b\b Done.\n')
% toc
% end

%% ------------------------------------------------------------------------
function dt = prunesumstat(file, output, opts)
% we read GWAS summary stat and match it to ld score files and flip the
% beta sign of those with A2/A1 matching. To read ld score files, we build
% index based on set_snpid_index module in parse.py function from polyfun
% package. 

% MATLAB writeparquet causes this error in pandas:
% UnicodeDecodeError: 'utf-8' codec can't decode byte 0xf6 in position 7: invalid start byte

if isfile(output + ".ss.parquet")
    dt = ".ss.parquet";
    fprintf('found summary stat file: %s\n', output + dt)
    return
elseif isfile(output + ".gz")
    dt = ".gz";
    fprintf('found summary stat file: %s\n', output + dt)
    return
end

t2 = tic;
fprintf('reading/matching GWAS summary stat file to ld score files...')

ss = readGWASfile(file, parallel=true, tall=false, full=true, ...
    light=false, n=true, workers=opts.workers, legacy=false);
ss.Properties.VariableNames = createGWASheader(colnames(ss));

% remove duplicated indel variants. This is because if we use
% allow_swapped_indel_alleles flag, set_snpid_index will throw an error
% that there are duplicated variants. One could use allow_duplicates flag
% as well, but it seems to be problematic when matching to ldmat file. I
% dedup these variants, and keep the one with lowest P-value.
idx_dup_indel = ss.A1 > ss.A2 & (strlength(ss.A1) > 1 | strlength(ss.A2) > 1);
index = ss.A1;
index(idx_dup_indel) = ss.CHR(idx_dup_indel) + "." + ss.BP(idx_dup_indel) + "." + ss.A1(idx_dup_indel) + "." + ss.A2(idx_dup_indel);
index(~idx_dup_indel) = ss.CHR(~idx_dup_indel) + "." + ss.BP(~idx_dup_indel) + "." + ss.A2(~idx_dup_indel) + "." + ss.A1(~idx_dup_indel);
index_dups = duplicates(index);
if ~isempty(index_dups)
    idx_dup = ismember(index, index_dups);

    % @27MARCH2023: for the time being, I decided to exclude all duplicated
    % indels since it's easier when matching to LDmatrix files. In future
    % this should be further implemeneted as there may be causal variants
    % among these indels.
    % idx_dup = find(ismember(index, index_dups));
    % idx_dup_logical = false(numel(index), 1);
    % index = index(idx_dup);
    % p_dup = string(ss.P(idx_dup));
    % p_dup_max = groupsummary(double(p_dup), index, @max);
    % p_dup_max_idx = ismember(p_dup, string(p_dup_max));
    % idx_dup_logical(idx_dup(p_dup_max_idx)) = true; % to be removed
    % 
    % % write removed indels to a file for LD matrix calculation. This helps
    % % to filter variants in LD matrix
    % rem_indels = ss(idx_dup_logical, :); % for LD matrix
    % [indels_pth, indels_file] = fileparts(output);
    % indels_file = fullfile(indels_pth, indels_file + ".indels.txt");
    % writetable(rem_indels, indels_file, Delimiter="\t")

    % ss(idx_dup_logical, :) = [];
    % clear idx_dup_logical idx_dup_indel
    
    minP = min(ss.P(idx_dup));
    ss(idx_dup, :) = [];
end

matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);

% read annotated variants
pq = getfilenames(makeWSLpath(opts.annwd, true), "parquet", fullpath=true).parquet; % only works with parquet files for now
apq = pq(pq.endsWith(".annot.parquet"));
ds = tall(parquetDatastore(apq, ...
    SelectedVariableNames=["CHR", "BP", "A1", "A2"], ...
    VariableNamingRule="preserve"));

ds.id = ds.CHR + "." + ds.BP + "." + ds.A1 + "." + ds.A2;
id_f = ss.CHR + "." + ss.BP + "." + ss.A1 + "." + ss.A2;
ss.id_b = ss.CHR + "." + ss.BP + "." + ss.A2 + "." + ss.A1;

% match overlapping variants
% id_f = gather(id_f);
idx_1 = ismember(id_f, ds.id);
if ~all(idx_1)  % match non-overlapping variants
    % ss = gather(ss);
    idx_2 = ismember(ss.id_b, ds.id);
    idx_2 = gather(idx_2);
    fprintf('\b\b Done (%.0f sec)\n', toc(t2))

    % if ~isempty(index_dups)
    %     fprintf('%d indels with duplicated ID (both A1/A2 and A2/A1 were present) were removed\n', numel(index_dups))
    % end
    if ~isempty(index_dups)
        fprintf('\t%d duplicated indels (A1/A2 and A2/A1 were present) were removed (min P = %.2g)\n', sum(idx_dup), minP)
    end

    % flip A1/A2 for these variants to match variants in ld2.score files
    fprintf('\t%d variants were flipped (A1/A2) to match ld2.score files\n', sum(idx_2))
    A1_flip = ss.A1(idx_2);
    ss.A1(idx_2) = ss.A2(idx_2);
    ss.A2(idx_2) = A1_flip;
    ss.A1FREQ(idx_2) = 1 - ss.A1FREQ(idx_2);
    ss.BETA(idx_2) = -ss.BETA(idx_2);
    ss.id_b = [];
    
    clear idx_1 idx_2 id_f A1_flip

    fprintf('writing to modified GWAS summary stats...')
    t1 = tic;
    % subset to each chr
    uchr = unique(ss.CHR);
    if isstring(uchr), uchr = natsort(uchr); end    
    outss = cell(numel(uchr), 1);
    for k = 1:numel(uchr)
        idx_2 = ss.CHR == uchr(k);
        outss{k} = ss(idx_2, :);
        ss(idx_2, :) = [];
    end

    % write to chunk files
    parfor k = 1:numel(outss)
        if k == 1
            wt = true;
        else
            wt = false;
        end
        writetable(outss{k}, output + "." + uchr(k) + ".txt", ...
            Delimiter="\t", WriteVariableNames=wt)
    end

    tmpsumstats = output + "." + (1:numel(outss)) + ".txt";
    mergeStatFiles(tmpsumstats, output)
    arrayfun(@delete, tmpsumstats)
    gzip(output)
    delete(output)

    % parquetwrite(output, ss)

    fprintf('\b\b Done (%.0f sec)\n', toc(t1))
    % dt = ".ss.parquet";
    dt = ".gz";

else
    % no need to flip
    ss2fuma(file, output=output)
    fprintf('\b\b Done (%.0f sec)\n', toc(t2))
    dt = ".gz";
end

matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);

end

% %% ------------------------------------------------------------------------
% function index = addSumstatIndex(ss)
% A1_first = ss.A1 < ss.A2 | strlength(ss.A2) > 1 | strlength(ss.A1) > 1;
% index = ss.A1;
% index(A1_first) = ss.CHR(A1_first) + "." + ss.BP(A1_first) + "." + ss.A1(A1_first) + "." + ss.A2(A1_first);
% index(~A1_first) = ss.CHR(~A1_first) + "." + ss.BP(~A1_first) + "." + ss.A2(~A1_first) + "." + ss.A1(~A1_first);
% end

%% ------------------------------------------------------------------------
% function opts = runLDSC(opts)
% opts.h2dir = fullfile(opts.pwd, opts.out + ".h2");
% if ~isfolder(opts.h2dir)
%     mkdir(opts.h2dir)
% end
% 
% opts.h2dirwsl = makeWSLpath(opts.h2dir);
% refwd = makeWSLpath(fullfile(pwd, 'LDSC.files', 'eur_ref_ld_chr/'));
% wldwd = makeWSLpath(fullfile(pwd, 'LDSC.files', 'eur_w_ld_chr/'));
% 
% cmd = "wsl python3 " + opts.polyfunwd + "ldsc.py" + ...
%     " --h2 " + opts.name + ...
%     " --ref-ld-chr " + refwd + ...
%     " --w-ld-chr " + wldwd + ...
%     " --out " + opts.h2dirwsl + "/h2" + ...
%     " --chisq-max 9999999.0";
% 
% tic
% fprintf('(step 1) estimating heritability with LDSC...')
% [~, ~] = system(cmd);
% res = readtable(fullfile(opts.enrichdir, 'enrichment.results'), ...
%     'TextType', 'string', 'FileType', 'text', ...
%     'VariableNamingRule', 'preserve');
% res = sortrows(res, 'Enrichment_p', 'ascend');
% save(fullfile(opts.enrichdir, 'enrichment.results.mat'), 'res')
% delete(fullfile(opts.enrichdir, 'enrichment.results'))
% fprintf('\b\b Done.\n')
% toc
% end