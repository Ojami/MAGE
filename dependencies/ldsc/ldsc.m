function ldsc(opts)
% A wrapper for ldsc python package: https://github.com/bulik/ldsc Note:
% follow the above link to install conda environment named "ldsc" In case
% of allele mapping between sumstats with 'rg', see:
% https://github.com/bulik/ldsc/issues/140#issuecomment-1031683241
% 
% Note: as explained here
% (https://github.com/bulik/ldsc/wiki/Partitioned-Heritability),
% --frqfile-chr is required when using --overlap-annot flag.
% 
% Note: for estimating the heritability and intercept, using partitioned
% LDSC is probably a better choice and better fit!
% https://www.nealelab.is/blog/2017/9/20/insights-from-estimates-of-snp-heritability-for-2000-traits-and-disorders-in-uk-biobank
% 
% Note: remove chi2 cap and at freq files when using partitioned LDSC:
% https://github.com/Nealelab/UKBB_ldsc_scripts/blob/master/ldsc_h2part_parallel_batch_v2_biomarkers.py
% 
% ./munge_sumstats.py \
% --sumstats /mnt/e/Incubator/polyfun/ALTn_INT_SNPxBMI.txt.gz \
% --out ALTn_INT_SNPxBMI
% --merge-alleles w_hm3.snplist
% 
% ./ldsc.py \
% --out ML_cluster \
% --h2 ML_cluster.sumstats.gz \
% --ref-ld-chr /mnt/e/Incubator/polyfun/LDSC.files/eur_ref_ld_chr/ \
% --w-ld-chr /mnt/e/Incubator/polyfun/LDSC.files/eur_w_ld_chr/ \
% --chisq-max 9999999.0 
% 
% ./munge_sumstats.py \
% --sumstats ${gwas} \
% --out PDFF
% --merge-alleles w_hm3.snplist
% 
% ldsc("h2", "ldsc_PDFF.munge.sumstats.gz", "out", "ldsc_PDFF", ...
% "ref_ld_chr", "F:\software\ldsc\eur_ref_ld_chr", "w_ld_chr", "F:\software\ldsc\eur_ref_ld_chr", "munge",false)
% 
% Oveis Jamialahmadi, GU, July 2022.

arguments
    opts.conda {mustBeTextScalar} = "ldsc" % name of the conda environment for ldsc
    opts.ldschome {mustBeFolder} % all codes should be present here. Also 'frqfile' and 'frqfile_chr' should be here
    opts.whome {mustBeFolder} % home directory for 'w-' arguments. If left empty, 'ldschome' is used
    opts.refhome {mustBeFolder} % home directory for 'ref-' arguments. If left empty, 'ldschome' is used
    opts.munge (1,1) logical = true % to munge sumstats
    opts.mungeonly (1,1) logical = false % only munge sumstats and ignore ldsc

    % munge_sumstats arguments
    opts.sumstats {mustBeFile} % can be > 1 file
    opts.mungemethod {mustBeTextScalar, mustBeMember(opts.mungemethod, ["ldsc", "matlab"])} = "ldsc" % method for munging of sumstats file 
    opts.N (1,1) double
    opts.N_cas (1,1) double
    opts.N_con (1,1) double
    opts.out {mustBeText}
    opts.info_min (1,1) double = 0.7
    opts.maf_min (1,1) double {mustBeInRange(opts.maf_min, 1e-12, 0.5)}
    opts.no_alleles (1,1) logical = false % Don't require alleles. Useful if only unsigned summary statistics are available and the goal is h2/partitioned h2 estimation rather than rg estimation.
    opts.merge_alleles (1,1) {mustBeFile} % e.g. w_hm3.snplist. Use this if INFO column is missing from 'sumstats'. see: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
    opts.chunksize (1,1) double = 500000
    opts.snp {mustBeTextScalar} % name of SNP column
    opts.N_col {mustBeTextScalar}
    opts.N_cas_col {mustBeTextScalar}
    opts.N_con_col {mustBeTextScalar}
    opts.a1 {mustBeTextScalar} % effect allele
    opts.a2 {mustBeTextScalar} % alternative allele
    opts.p {mustBeTextScalar}
    opts.signed_sumstats {mustBeTextScalar}
    opts.info {mustBeTextScalar}
    opts.a1_inc (1,1) logical = false % set to true only if beta/z/or is absent from sumstats and all betas are a1 oriented (i.e. either all are + or -). In this case Z sign is either + or -.
    opts.ignore {mustBeTextScalar}

    % ldsc (not the complete list of arguments)
    opts.annot {mustBeTextScalar} %  Filename prefix for annotation file for partitioned LD Score estimation. LDSC will automatically append .annot or .annot.gz to the filename prefix. See docs/file_formats_ld for a definition of the .annotformat.
    opts.thin_annot (1,1) logical = false % This flag says your annot files have only annotations, with no SNP, CM, CHR, BP columns.
    opts.maf (1,1) double {mustBeInRange(opts.maf, 0, 0.5)} %  Minor allele frequency lower bound. Default is MAF > 0.
    opts.h2 {mustBeFile} % Filename for a .sumstats[.gz] file for one-phenotype LD Score regression. --h2 requires at minimum also setting the --ref-ld and --w-ld flags.
    opts.h2_cts {mustBeFile} % Filename for a .sumstats[.gz] file for cell-type-specific analysis. --h2-cts requires the --ref-ld-chr, --w-ld, and --ref-ld-chr-cts flags.
    opts.rg {mustBeFile} % Comma-separated list of prefixes of .chisq files for genetic correlation estimation.
    opts.ref_ld {mustBeText} % Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD Score regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix.
    opts.ref_ld_chr {mustBeText} % Same as --ref-ld, but will automatically concatenate .l2.ldscore files split across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix. If the filename prefix contains the symbol @, LDSC will replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome numbers to the end of the filename prefix.Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gzExample 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz
    opts.w_ld {mustBeTextScalar} % Filename prefix for file with LD Scores with sum r^2 taken over SNPs included in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz.
    opts.w_ld_chr {mustBeTextScalar} %  Same as --w-ld, but will read files split into 22  chromosomes in the same manner as --ref-ld-chr.
    opts.overlap_annot (1,1) logical = false % This flag informs LDSC that the partitioned LD Scores were generates using an annot matrix with overlapping  categories (i.e., not all row sums equal 1), and prevents LDSC from displaying output that is meaningless with overlapping categories.
    opts.print_coefficients (1,1) logical = false % when categories are overlapping, print coefficients as well as heritabilities.
    opts.no_intercept (1,1) logical = false % If used with --h2, this constrains the LD Score regression intercept to equal 1. If used with --rg, this constrains the LD Score regression intercepts for the h2 estimates to be one and the intercept for the genetic covariance estimate to be zero.
    opts.intercept_h2 double {mustBeVector} %  Intercepts for constrained-intercept single-trait LD Score regression.
    opts.intercept_gencov double {mustBeVector} % Intercepts for constrained-intercept cross-trait LD  Score regression. Must have same length as --rg. The first entry is ignored.
    opts.two_step (1,1) double %  Test statistic bound for use with the two-step estimator. Not compatible with --no-intercept and --constrain-intercept. Note: The current package by default imposes a cutoff value of 30 on the test statistic in two-step estimation. We removed this limit by setting its default value to infinity. As a result, SNPs with high chi^2 values are no longer filtered out.
    opts.chisq_max (1,1) double % Max chi^2. Based on https://search.r-project.org/CRAN/refmans/bigsnpr/html/snp_ldsc.html, this is equivalent to 'two-step' argument.
    opts.ref_ld_chr_cts {mustBeText} % Name of a file that has a list of file name prefixes for cell-type-specific analysis.
    opts.print_cov (1,1) logical = false %  For use with --h2/--rg. This flag tells LDSC to print the covaraince matrix of the estimates.
    opts.chunk_size (1,1) double % Chunk size for LD Score calculation. Use the default.
    opts.not_M_5_50 (1,1) logical = false %  This flag tells LDSC to use the .l2.M file instead of the .l2.M_5_50 file.
    opts.no_check_alleles (1,1) logical = false % For rg estimation, skip checking whether the alleles match. This check is redundant for pairs of chisq files generated using munge_sumstats.py and the same argument to the --merge-alleles flag.
    opts.samp_prev (1,1) double % Sample prevalence of binary phenotype (for conversion to liability scale).
    opts.pop_prev (1,1) double % Population prevalence of binary phenotype (for conversion to liability scale).
    opts.frqfile {mustBeTextScalar} % For use with --overlap-annot. Provides allele frequencies to prune to common snps if --not-M-5-50 is not set.
    opts.frqfile_chr {mustBeTextScalar} % Prefix for --frqfile files split over chromosome.
end

% non-LDSC arguments
nativefis = ["munge", "mungemethod", "ldschome", "conda", "whome", "refhome", "mungeonly"];

% munge arguments
fd = ["N", "N_cas", "N_con", "out", "info_min", "maf_min", "no_alleles",...
    "merge_alleles", "chunksize", "snp", "N_col", "N_cas_col", "sumstats", ...
    "N_con_col", "a1", "a2", "p", "signed_sumstats", "info", "a1_inc", "ignore"];

% prepare inputs ----------------------------------------------------------
if opts.mungeonly, opts.munge = true; end % override 'munge' flag

% if opts.munge && ~isfield(opts, 'sumstats')
%     error("ldsc:'sumstats' cannot be left empty with 'munge' flag!")
% end

if ~isfield(opts, 'ldschome')
    opts.ldschome = pwd; 
end

if ~isfield(opts, 'whome'), opts.whome = opts.ldschome; end
if ~isfield(opts, 'refhome'), opts.refhome = opts.ldschome; end

% 'ref-' and 'w-' must be provided
if ~opts.mungeonly
    if all(~isfield(opts, ["w_ld", "w_ld_chr"]))
        error("ldsc: either of w-ld or w-ld-chr arguments should be used!")
    end
    if all(~isfield(opts, ["ref_ld", "ref_ld_chr"]))
        error("ldsc: either of ref-ld or ref-ld-chr arguments should be used!")
    end

    fi = ["ref_ld_chr_cts", "w_ld", "w_ld_chr", "ref_ld", "ref_ld_chr", "frqfile", "frqfile_chr"];
    for i = 1:numel(fi)
        if isfield(opts, fi(i))
            tmpfiles = opts.(fi(i));
            for k = 1:numel(tmpfiles)
                if fileparts(tmpfiles(k)) == ""
                    tmpfiles(k) = fullfile(opts.whome, tmpfiles(k));
                end
                
                % check LD file patterns
%                 if isfolder(tmpfiles(k)) % inputs can be in form of path/to/file/cell_type.1 pattern
%                     if startsWith(fi(i), "frqfile") 
%                         cpat = getfilenames(tmpfiles(k), "frq").frq;
%                     else
%                         cpat = getfilenames(tmpfiles(k), "gz").gz;
%                     end
%                 end
                
                if ~endsWith(tmpfiles(k), filesep) && isfolder(tmpfiles(k))
                    tmpfiles(k)  = tmpfiles(k) + filesep;
                end

                % convert relative paths to absolute paths in
                % ref_ld_chr_cts dictionary if necessary
                if fi(i) == "ref_ld_chr_cts"
                    cts = readmatrix(tmpfiles(k), Delimiter=["\t", ","], FileType="text", OutputType="string");
                    cts(:, 2:3) = arrayfun(@makeWSLpath, fullfile(opts.whome, cts(:, 2:3)));
                    [pth, name, ext] = fileparts(makeWSLpath(tmpfiles(k), true));
                    tmpfiles(k) = fullfile(pth, name + getRandomName("", 8) + ext);
                    cts(:, 2) = cts(:,2) + "," + cts(:,3);
                    cts(:, 3) = [];
                    writematrix(cts, tmpfiles(k), Delimiter="\t", QuoteStrings="none", FileType="text");
                    dos2unix(tmpfiles(k))
                end

%                 if isfolder(tmpfiles(k))
%                     if any(isnan(double(cpat.extractBefore(".")))) % chrom name
%                         tmpfiles(k) = fullfile(tmpfiles(k), strshare(cpat));
%                     end
%                 end
    
                tmpfiles(k) = makeWSLpath(tmpfiles(k));
            end
            opts.(fi(i)) = tmpfiles; % join(tmpfiles, ",");
        end
    end

%     fi = ["ref_ld", "ref_ld_chr"];
%     for i = 1:numel(fi)
%         if isfield(opts, fi(i))
%             if fileparts(opts.(fi(i))) == ""
%                 opts.(fi(i)) = fullfile(opts.refhome, opts.(fi(i)));
%             end
%             opts.(fi(i)) = makeWSLpath(opts.(fi(i)));
%             if ~endsWith(opts.(fi(i)), "/"), opts.(fi(i)) = opts.(fi(i)) + "/"; end
%         end
%     end
% 
%     fi = ["frqfile", "frqfile_chr"];
%     for i = 1:numel(fi)
%         if isfield(opts, fi(i))
%             opts.(fi(i)) = fullfile(opts.ldschome, opts.(fi(i)));
%             opts.(fi(i)) = makeWSLpath(opts.(fi(i)));
%             if ~endsWith(opts.(fi(i)), "/"), opts.(fi(i)) = opts.(fi(i)) + "/"; end
%         end
%     end
end

if ~opts.mungeonly % no 'sumstats' file
    if isfield(opts, 'h2'), opts.sumstats = opts.h2;
    elseif isfield(opts, 'h2_cts'), opts.sumstats = opts.h2_cts;
    elseif isfield(opts, 'rg'), opts.sumstats = opts.rg; 
    else, error("ldsc:any of h2, h2_cts or rg should be used!");
    end
end

if ~isfield(opts, 'out')
    [~, opts.out] = fileparts(opts.sumstats);
elseif opts.munge && (numel(opts.out) ~= numel(opts.sumstats))
    error("each element of out should correspond to each element in 'sumstats' with 'munge' option!")
elseif ~opts.munge && numel(opts.out) > 1
    opts.out = opts.out(1);
end

% start the analysis ------------------------------------------------------
if opts.munge
    for i = 1:numel(opts.sumstats)
        optsin = opts;
        optsin.sumstats = opts.sumstats(i);
        optsin.out = opts.out(i);
        opts.sumstats(i) = munge_sumstats(optsin, fd);
    end
end

if opts.mungeonly
    return
elseif ~isfield(opts, 'rg') && numel(opts.sumstats) > 1
    error("ldsc: only 'rg' option can be used with multiple sumstat files!")
end % ignore ldsc analysis

% run LDSC
if isfield(opts, 'h2')
    opts.out = opts.out + ".h2";
    if opts.munge, opts.h2 = opts.sumstats; end
elseif isfield(opts, 'h2_cts')
    opts.out = opts.out + ".h2_cts";
    if opts.munge, opts.h2_cts = opts.sumstats; end
elseif isfield(opts, 'rg')
    opts.out = opts.out + ".rg";
    if opts.munge, opts.rg = opts.sumstats; end
else
    error("ldsc:any of h2, h2_cts or rg should be used!")
end

tmpname = getRandomName("ldsc", 5);
fi = string(fieldnames(opts));
fi = union("out", setdiff(fi, union(fd, nativefis)));
cmd = makeCMD(opts, fi);
cmd = makeWSLpath(fullfile(opts.ldschome, "ldsc.py")) + " " + cmd;
runbash(cmd, tmpname, "conda", opts.conda, "parallel", false, "verbose", true)

if isfield(opts, 'ref_ld_chr_cts') % delete temp file
    delete(makeWSLpath(opts.ref_ld_chr_cts, true))
end

end % END

%% subfunctions ===========================================================
function sumstats = munge_sumstats(opts, fd)

opts.out = opts.out + ".munge";
tmpname = getRandomName("munge", 5);
hdr = bfilereader(opts.sumstats, "header", true, "summary", "firstline");
delsum = false;
if any(contains(hdr.lower, "log10p")) && ~isfield(opts, "p") % regenie
    fprintf("REGENIE file: converting log10p to P...\n")
    pcol = find(contains(hdr.lower, "log10p"));
    opts.p = hdr(pcol(1));
    opts.snp = "ID";
    opts.a1 = "ALLELE1";
    opts.a2 = "ALLELE0";
    sumstat_home = fileparts(opts.sumstats);
    if strcmp(sumstat_home, ""), sumstat_home = pwd; end
    newsumstats = makeWSLpath(fullfile(sumstat_home, opts.out + "_tmp"));

    % convert log10 P -> P
    runbash("awk 'NR!=1 {$"+ pcol(1) + "=10^-$" + pcol(1) + ...
        "}1' " + makeWSLpath(opts.sumstats) + " > " + newsumstats, ...
        tmpname, "wait", true);
    opts.sumstats = fullfile(sumstat_home, opts.out + "_tmp");
    delsum = true; % get rid of temprorary sumstat
end

cmd = makeCMD(opts, fd);

cmd = makeWSLpath(fullfile(opts.ldschome, "munge_sumstats.py")) + " " + cmd;
runbash(cmd, tmpname, "conda", opts.conda, "parallel", false, ...
    "verbose", true, "delete", true)

if delsum, delete(opts.sumstats); end % get rid of temprorary sumstat
sumstats = opts.out + ".sumstats.gz"; % munged
end

%% ------------------------------------------------------------------------
function cmd = makeCMD(opts, fd)
% fd: fields to keep

cmd = "";
for i = 1:numel(fd)
    if ~isfield(opts, fd(i)) || isempty(opts.(fd(i)))
        continue
    end

    if islogical(opts.(fd(i))) % FLAG
        if opts.(fd(i))
            cmd = join([cmd, "--" + replace(fd(i), '_', '-')], " ");
        end
    else
        if isstring(opts.(fd(i))) && all(isfile(opts.(fd(i)))) % files should be converted to linux path (if pc)
            tmpval = string(arrayfun(@makeWSLpath, opts.(fd(i))));
        else
            tmpval = string(opts.(fd(i)));
        end    
        
        for j = 1:numel(tmpval)
            if contains(tmpval(j), ["/", "\"]) % ? path to file
                tmpval(j) = '"' + tmpval(j) + '"';
            end
        end

        if numel(tmpval) > 1 % input files
            tmpval =  join(tmpval, ",");
        end

        cmd = join([cmd, "--" + replace(fd(i), '_', '-'), tmpval], ' ');
        
    end
end
end