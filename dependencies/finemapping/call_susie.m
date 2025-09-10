function call_susie(sumstats_file, ld_file_pref, opts)
% A wrapper for SuSiE fine-mapping.
% call_susie: Minimal FinnGen-style SuSiE-RSS runner (raw-binary LD by prefix)
% Reads GWAS summary stats (CSV/TSV), loads an in-sample LD correlation matrix
% from raw binary + shape files using a single prefix, and runs SuSiE-RSS.
% Extracts 95% credible sets, computes post-hoc purity (min/mean r^2), and
% writes ONE tab-delimited file = input sumstats with appended SuSiE columns.
%
% FinnGen reference:
% https://github.com/FINNGEN/finemapping-pipeline/blob/master/R/run_susieR.R
%
% REQUIRED sumstats columns: chr, pos, a1, a2, beta, se, n
% Oveis Jamialahmadi, University of Gothenburg, 28 Aug 2025.

arguments
    sumstats_file {mustBeFile} % path to summary statistics file for the fine-mapping region/locus.
    ld_file_pref {mustBeTextScalar} % path to LD prefix (e.g. /path/to/ld, where ld.bin and ld.shape are in /path/to/)
    opts.out_dir {mustBeTextScalar} = pwd() % output directory (will be created, default is pwd)
    opts.L (1,1) double = 10 % max single-effects
    opts.min_cs_corr (1,1) double {mustBeInRange(opts.min_cs_corr, 0, 1)} = 0.5   % |r| purity threshold at extraction (default 0.5)
    opts.good_cred_r2 (1,1) double {mustBeInRange(opts.good_cred_r2, 0, 1)} = 0.25 % post-hoc CS purity threshold on r^2 (default 0.25 ≈ |r|≥0.5) -> FINNGEN
    opts.estimate_residual_variance (1,1) logical = true % true is recommended for in-sample LD
    opts.estimate_prior_variance (1,1) logical = true % let SuSiE estimate prior variance
    opts.psd_check (1,1) logical = false % if TRUE, verify R is PSD after symmetrization and repair if not (nearPD first, else eigen-floor)
    opts.r_tol (1,1) double = 1e-8 % PSD eigenvalue tolerance for the check. Default 1e-8
    opts.rfun {mustBeFile} = fullfile(fileparts(which("call_susie.m")), "call_susie.R")
    opts.prior_weights {mustBeFolder} % path to polyfun estimated prior weights (*_constrained.gz files should be there) 
    
end

% output directory
opts.out_dir = string(opts.out_dir);
if ~isfolder(opts.out_dir), mkdir(opts.out_dir); end

% name of the region
t1 = tic;
[~, regname] = fileparts(sumstats_file);
fprintf("SuSiE fine-mapping for region: %s\n", regname)

% check prior_wieghts and append SNPVAR as additional column to summary
% stats. Note that only subset of variants in polyfun files will be kept
if isfield(opts, "prior_weights")
    sumstats_file = appendPriorWeights(sumstats_file, opts);
    opts = rmfield(opts, "prior_weights");
end

% make R-friendly paths
opts.out_dir = opts.out_dir.replace(filesep, "/");
sumstats_file = sumstats_file.replace(filesep, "/");
ld_file_pref = ld_file_pref.replace(filesep, "/");

% check if output file already exists
[~, sufile, ext] = fileparts(ld_file_pref);
sufile = fullfile(opts.out_dir, sufile + ext + ".susie.txt");
if isfile(sufile)
    fprintf("%s already exists! Skipped!\n", sufile)
    t1 = seconds(toc(t1));
    fprintf('\t[done in %.2f min]\n\n', minutes(t1))
    return
end

r = string;
r(1) = "setwd('" + opts.out_dir + "')";
r(2) = "source('" + opts.rfun.replace("\", "/") + "')";
ropts = rmfield(opts, "rfun");
ropts.ld_file_pref = ld_file_pref;
ropts.sumstats_file = sumstats_file;
tmp_rstr = struct2rfun(ropts,par=true, pretty=true);
tmp_rstr(1) = "call_susie" + tmp_rstr(1);
r = [r, tmp_rstr'];

code_name = getRandomName("susierss", 6);

MATLAB2Rconnector(fullfile(opts.out_dir, code_name + ".r"), code=r);

% make sure if output table is generated, otherwise inspect log file
if ~isfile(sufile)
    error("call_susie failed! inspect the log file!")
end

t1 = seconds(toc(t1));
fprintf('\t[done in %.2f min]\n\n', minutes(t1))

end % END

%% subfunctions ===========================================================
function file = appendPriorWeights(sfile, opts)
% uses prior weights from polyfun and append to sfile two additional
% columns: SNPVAR, idx; where latter is the index of variants present in
% polyfun files. This is to avoid any need to modify LD matrices and subset
% rows/cols when necessary, ergo new_R <- R[idx, idx]

pfile = getfilenames(opts.prior_weights, "gz").gz;
pfile(~pfile.endsWith("_constrained.gz")) = [];
assert(~isempty(pfile), "no constrained.gz files were found in prior_weights dir!")

% find chromosome
pchrom = pfile.extractBetween(".", ".snpvar_ridge_constrained.gz");
[pth, file, ext] = fileparts(sfile);
pchrom_idx = pchrom == file.extractBetween(textBoundary("start") + "chr", "_");
pfile(~pchrom_idx) = [];
assert(~isempty(pfile), "no constrained.gz files was found for this chromosome!")

file = fullfile(pth, file + ".snpvar" + ext); % new sfile

if isfile(file)
    fprintf("[SNPVAR]: found SNPVAR file for this sumstat (skipped)\n")
    return
end

% read polyfun file
gunzip(fullfile(opts.prior_weights, pfile), opts.prior_weights)
pfile = fullfile(opts.prior_weights, pfile.erase(".gz" + textBoundary("end")));

dopts = detectImportOptions(pfile, TextType="string", ...
    FileType="text", VariableNamingRule="preserve");
dopts.SelectedVariableNames = ["CHR", "SNP", "BP", "A1", "A2", "SNPVAR"];
pdf = readtable(pfile, dopts);
delete(pfile)

% read summary stat file and match to SNP/A1/A2 in polyfun file
sdf = readtable(sfile, TextType="string", VariableNamingRule="preserve");
sdf.idx = (1:height(sdf))';
sdf.idf = sdf.SNP + ":" + sdf.A1 + ":" + sdf.A2;
sdf.idb = sdf.SNP + ":" + sdf.A2 + ":" + sdf.A1;
pdf.idf = pdf.SNP + ":" + pdf.A1 + ":" + pdf.A2;

keep_idx = ismember(sdf.idb, pdf.idf); % flipped in pdf
assert(~any(ismember(sdf.idf(keep_idx), pdf.idf)), "failed to match variants!")
sdf.idf(keep_idx) = sdf.idb(keep_idx);
sdf.idb = [];
sdf.SNPVAR(:) = nan;
[idx1, idx2] = ismember(sdf.idf, pdf.idf); 
sdf.SNPVAR(idx1) = pdf.SNPVAR(idx2(idx1));


rm_idx = ismissing(sdf.SNPVAR);
if any(rm_idx)
    sdf(rm_idx, :) = [];
    fprintf("[INFO]: removed %d variants from sumstat without SNPVAR.\n", nnz(rm_idx))
end

sdf.idf = [];
writetable(sdf, file, Delimiter="\t");

end % END

