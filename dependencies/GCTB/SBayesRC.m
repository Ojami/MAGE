function SBayesRC(gctab, opts)
%@08AUG2024: a wrapper for GCTB SBayesRC:
% https://cnsgenomics.com/software/gctb/#SBayesRCTutorial

arguments
    gctab {mustBeFile}
    opts.ldm {mustBeFolder} = fullfile(fileparts(which("SBayesRC.m")), "ukbEUR_Imputed")% a folder comprising the eigen-decomposition data for each block.
    opts.annot {mustBeFile} = fullfile(fileparts(which("SBayesRC.m")), "annot_baseline2.2.txt") % annotation file, with columns being SNP ID, Intercept (a column of one), and annotation values. 
    opts.parallel (1,1) logical = true % for readGWASfile only
    opts.workers (1,1) double = 30 % for readGWASfile only
    opts.threads (1,1) double = 15 % for GCTA
    opts.gctb {mustBeFile} = fullfile(fileparts(which("SBayesRC.m")), "gctb") % gctb executable binary file
    opts.output {mustBeTextScalar} % output file name, if left empty, will be same as the input 'gctab'

end

% check parallel pool
if opts.parallel && isempty(gcp('nocreate'))
    if isnan(opts.workers)
        % number of physical cores
        maxNumCompThreads('automatic'); % reset (only if was changed in current session)
        opts.workers = maxNumCompThreads; 
    end
    parc = parcluster('Processes');
    parc.NumWorkers = opts.workers;
    parpool(parc); % should I delete in the end? 
end

% read summary stat file
if isfield(opts, 'output')
    output = opts.output;
else
    [hm, output] = fileparts(gctab);
    output = string(fullfile(hm, output));
end
opts.outdir = string(fileparts(output));
if opts.outdir == "", opts.outdir = pwd; end
if ~endsWith(opts.outdir, filesep), opts.outdir = opts.outdir + filesep; end

total_time = tic;

name = getRandomName("gctb"); % temp GCTB files

gctab = readGWASfile(gctab, 'parallel', opts.parallel, ...
    'full', true, 'light', false, 'n', true);

% Since p-value cannot go under realmin, modify 0 p-values based on
% the corresponding Z-scores 
f_zero = gctab.p <= realmin;
if any(f_zero)
    Z = abs(gctab.se(f_zero)./gctab.beta(f_zero)); % 1/Z
    gctab.p(f_zero) = realmin.*(Z./max(Z));
end

gctab.pos = []; % not needed 

chru = unique(gctab.chr); % must be numeric (i.e. X chr should be 23)
newcols = ["A1", "A2", "freq", "b", "SNP", "N"];
gctab = renamevars(gctab, ["allele1", "allele0", "afreq", "beta", ...
    "annotate", "n"], newcols);
gctab = gctab(:, ["chr", "SNP", "A1", "A2", "freq", "b", "se", "p", "N"]);

fprintf('splitting the GWAS data into chromosome tables...')
t1 = tic;
chrtab = cell(numel(chru), 1);
for k = 1:numel(chru)
    idx = gctab.chr == chru(k);
    chrtab{k} = gctab(idx, :);
    gctab(idx, :) = [];
    chrtab{k}.chr = [];
end
fprintf('\b\b done (%.0f sec)\n', toc(t1))

fprintf('writing to temp .ma files...')
t1 = tic;
outdir = opts.outdir;
if opts.parallel
    parfor k = 1:numel(chru)
        writetable(chrtab{k}, fullfile(outdir, name + "." + chru(k) + ".ma"), FileType='text', Delimiter=char(9));
        dos2unix(fullfile(outdir, name + "." + chru(k) + ".ma"))
    end
else
    for k = 1:numel(chru)
        writetable(chrtab{k}, fullfile(outdir, name + "." + chru(k) + ".ma"), FileType='text', Delimiter=char(9));
        dos2unix(fullfile(outdir, name + "." + chru(k) + ".ma"))
    end
end

fprintf('\b\b done (%.0f sec)\n', toc(t1))

% get ldm blocks
linfo = tabularTextDatastore(fullfile(opts.ldm, "ldm.info"), ...
    TextType="string", VariableNamingRule="preserve", ...
    SelectedVariableNames=["Block", "Chrom"], FileExtensions=".info");
linfo = readall(linfo);

opts.gctb = makeWSLpath(opts.gctb);
opts.annot = makeWSLpath(opts.annot);
opts.ldm = makeWSLpath(opts.ldm);

% Step 1: QC & imputation of summary statistics for SNPs in LD reference but not in GWAS data.
% gctb --ldm-eigen ldm --gwas-summary test.ma --impute-summary --out test --thread 4
fprintf('Running step 1: QC & imputation for %d blocks\n', height(linfo))
chru = unique(linfo.Chrom);
for k = 1:numel(chru)

    % make sure the summary stat for this chrom exists
    if ~isfile(fullfile(opts.outdir, name + "." + chru(k) + ".ma"))
        continue
    end

    sfile = makeWSLpath(opts.outdir) + name + "." + chru(k) + ".ma";

    % use --gwas-summary for each chrom and ldm blocks for that chrom
    % using linfo table
    blocksToUse = linfo.Block(linfo.Chrom == chru(k));
    cmd = strings(numel(blocksToUse), 1);

    t1 = tic;
    fprintf("\tchr %s (%d blocks)...", string(chru(k)), numel(blocksToUse))
    for j = 1:numel(blocksToUse)
        cmd(j) = opts.gctb + ...
            " --ldm-eigen " + opts.ldm +...
            " --gwas-summary " + sfile + ...
            " --impute-summary --out " + makeWSLpath(opts.outdir) + name +...
            " --block " + blocksToUse(j);
    end
    runbash(cmd, parallel=true, wait=true, verbose=false);
    delete(fullfile(opts.outdir, name + "." + chru(k) + ".ma"))
    fprintf('\b\b done (%.0f sec)\n', toc(t1))
   
end

t1 = tic;
imfiles = getfilenames(opts.outdir, "ma").ma;
imfiles(~imfiles.startsWith(name + ".block") | ~imfiles.endsWith(".imputed.ma")) = [];
imfiles = fullfile(opts.outdir, imfiles);
fprintf('Merging %d block files..', numel(imfiles))
mergeFiles(imfiles, output=fullfile(opts.outdir, name + ".imputed.ma"), delimiter=" ");
dos2unix(fullfile(opts.outdir, name + ".imputed.ma"));
arrayfun(@delete, imfiles)
fprintf('\b\b done (%.0f sec)\n', toc(t1))

t1 = tic;
fprintf('Step2: running sbayesRC..')
cmd = opts.gctb + " --ldm-eigen " + opts.ldm + " --gwas-summary " +...
    makeWSLpath(fullfile(opts.outdir, name + ".imputed.ma")) + ...
    " --sbayes RC --annot " + opts.annot + " --out " + ...
    makeWSLpath(fullfile(opts.outdir, name)) + " --thread " + opts.threads;
runbash(cmd, parallel=false, wait=true, verbose=true);

delete(fullfile(opts.outdir, name + ".imputed.ma"))
fprintf('\b\b done (%.0f sec)\n', toc(t1))

% rename outputs
ext = "." + ["badSNPlist", "parRes", "parSetRes", "snpRes"];
for k = 1:numel(ext)
    if isfile(fullfile(opts.outdir, name + ext(k)))
        movefile(fullfile(opts.outdir, name + ext(k)), opts.output + ext(k))
    end
end

fprintf("All jobs done in %s\n", string(seconds(toc(total_time)), "hh:mm:ss"))

end % END