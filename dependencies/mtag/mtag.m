function mtag(opts)
% a wrapper for mtag (Multi-Trait Analysis of GWAS) Python package.
% https://github.com/JonJala/mtag
% 
% Note: mtag requires a conda environment, default env is "ldsc" which uses
% Python 2.6 and also contains ldsc package. You need to install joblib
% (pip install joblib) within the 'ldsc' env.

arguments
    opts.sumstats {mustBeFile}
    opts.out {mustBeTextScalar} = "./mtag_results" %  Specify the directory and name prefix to output MTAG results. All mtag results will be prefixed with the corresponding tag. Default is ./mtag_results
    opts.make_full_path (1,1) logical = true %  option to make output path specified in -out if it  does not exist.
    opts.snp_name {mustBeTextScalar} %  Name of the single column that provides the unique identifier for SNPs in the GWAS summary statistics  across all GWAS results. Default is "snpid". This the index that will be used to merge the GWAS summary statistics. Any SNP lists passed to ---include or --exclude should also contain the same name.
    opts.z_name {mustBeTextScalar} % The common name of the column of Z scores across all input files. Default is the lowercase letter z.
    opts.beta_name {mustBeTextScalar} % overrides z_name
    opts.se_name {mustBeTextScalar}
    opts.n_name {mustBeTextScalar}
    opts.n_value double {mustBeVector} % Comma separated sample size values for each GWAS summary statistics files. This option is useful for GWAS input that does not include an N column, e.g. BOLT-LMM.
    opts.eaf_name {mustBeTextScalar} % effect allele freq
    opts.chr_name {mustBeTextScalar}
    opts.bpos_name {mustBeTextScalar}
    opts.a1_name {mustBeTextScalar}
    opts.a2_name {mustBeTextScalar}
    opts.p_name {mustBeTextScalar}

    % filter options
    opts.maf_min (1,1) double % set the threshold below SNPs with low minor allele frequencies will be dropped. Default is 0.01. Set to 0 to skip MAF filtering.
    opts.n_min (1,1) double = 0 %  set the minimum threshold for SNP sample size in input data. Default is 2/3*(90th percentile). Any SNP that does not pass this threshold for all of the GWAS input statistics will not be included in MTAG.
    opts.info_min (1,1) double %  Minimim info score for filtering SNPs for MTAG.
    opts.incld_ambig_snps (1,1) logical %  Include strand ambiguous SNPs when performing MTAG. by default, they are not used to estimate Omega or Sigma.

    % Special Cases:
    opts.use_beta_se (1,1) logical = false
    opts.no_overlap (1,1) logical % Imposes the assumption that there is no sample overlap between the input GWAS summary statistics. MTAG is performed with the off-diagonal terms on the residual covariance matrix set to 0.
    opts.perfect_gencov (1,1) logical % Imposes the assumption that all traits used are perfectly genetically correlated with each other. The off-diagonal terms of the genetic covariance matrix are set to the square root of the product of the heritabilities
    opts.equal_h2 (1,1) logical % Imposes the assumption that all traits passed to MTAG have equal heritability. The diagonal terms of the genetic covariance matrix are set equal to each other. Can only be used in conjunction with --perfect_gencov
   
    % misc
    opts.ld_ref_panel {mustBeFolder}
    opts.verbose (1,1) logical = true
    opts.stream_stdout (1,1) logical = true
    opts.force (1,1) logical = false % Force MTAG estimation even though the mean chi2 is small.

    % native options
    opts.conda {mustBeTextScalar} = "ldsc" % name of the conda environment for mtag
    opts.mtaghome {mustBeFolder} = fullfile(fileparts(which("mtag.m")), "mtag"); % all codes should be present here.
    opts.calcz (1,1) logical = true % at the moment, 'use_beta_se' has been deactivated, so Z should be calculated internally.
    opts.log10p (1,1) logical = true % effective only for REGENIE's sumstats
    opts.mungeonly (1,1) logical = false % only perform munging
    opts.use_beta_as_z (1,1) logical = false % use beta column as Z (if pre-processed files are present)
end

% non-mtag arguments
nativefis = ["conda", "mtaghome", "calcz", "log10p", "mungeonly", "use_beta_as_z"];

% prepare sumstats
for i = 1:numel(opts.sumstats)
    optsin = opts;
    optsin.sumstats = opts.sumstats(i);
    optsin = munge_sumstats(optsin);
    opts.sumstats(i) = optsin.sumstats;
    fi = ["snp", "z", "a1", "a2", "beta", "se", "p", "n", "eaf", "chr", "bpos"] + "_name";
    if i == 1 % all traits should have the same columns
        for j = 1:numel(fi)
            if isfield(optsin, fi(j)) && ~isfield(opts, fi(j))
                opts.(fi(j)) = optsin.(fi(j));
            end
        end
    end
end

% remove beta/se if Z is present
if isfield(opts, "z_name")
    opts = rmfield(opts, ["beta_name", "se_name"]);
end

if opts.mungeonly, return, end

tmpname = getRandomName("mtag", 5);
fi = string(fieldnames(opts));
fi = setdiff(fi, nativefis);
cmd = makeCMD(opts, fi);
cmd = makeWSLpath(fullfile(opts.mtaghome, "mtag.py")) + " " + cmd;
runbash(cmd, tmpname, "conda", opts.conda, "parallel", false, "verbose", true)


end % END

%% subfunctions ===========================================================
function opts = munge_sumstats(opts)

fid = fopen(opts.sumstats, 'r');
hdr = fgetl(fid); hdr = split(string(hdr));
fclose(fid);

if any(contains(hdr, "LOG10P")) || any(hdr == "ID")
    tool = "REGENIE";
    pcol = find(contains(hdr.lower, "log10p"));
    opts.p_name = hdr(pcol(1));
    opts.snp_name = "ID";
    opts.a1_name = "ALLELE1";
    opts.a2_name = "ALLELE0";
    opts.se_name = "SE";
    opts.chr_name = "CHROM";
    opts.bpos_name = "GENPOS";
    opts.n_name = "N";

elseif any(hdr == "A1FREQ")
    tool = "BOLT";
    opts.p = "P_BOLT_LMM";
    if any(hdr == "P_BOLT_LMM_INF")
        opts.p = "P_BOLT_LMM_INF";
    end
    opts.snp_name = "SNP";
    opts.chr_name = "CHR";
    opts.bpos_name = "BP";
    opts.a1_name = "ALLELE1";
    opts.a2_name = "ALLELE0";

elseif any(hdr == "AF_Allele2")
    tool = "SAIGE";
    opts.snp_name = hdr(contains(hdr.lowr, "rsid"));
    opts.chr_name = "CHR";
    opts.bpos_name = "POS";
    opts.p = "p_value";
    opts.a1_name = "Allele2";
    opts.a2_name = "Allele1"; 

else % unknown or already formatted
    return
end

opts.eaf_name = hdr(contains(hdr.lower, ["af_allele2", "a1freq", "a1_freq", "maf", "a2freq"]));
opts.beta_name = hdr(startsWith(hdr.lower, "beta"));

tmpname = getRandomName("munge", 5);
if tool == "REGENIE" && opts.log10p
    
    [~, name] = fileparts(opts.sumstats);
    newsumstats = makeWSLpath(fullfile(pwd, name + ".mtag.txt"));

    % find delimiter
    [~, sep] = bfilereader(opts.sumstats, summary="firstline");

    % convert log10 P -> P
    tt = tic;
    fprintf('converting log10p to p for REGENIE sumstats...')
    if ~isfile(fullfile(pwd, name + ".mtag.txt"))
        runbash("awk -F'" + sep.sep + "' -v OFS='" + sep.sep + ...
            "' 'NR!=1 {$"+ pcol(1) + "=10^-$" + pcol(1) + ...
            "}1' " + makeWSLpath(opts.sumstats) + " > " + newsumstats, ...
            tmpname, "wait", true);
    end
    fprintf('\b\bdone (%.2f sec)\n', toc(tt))
    opts.sumstats = fullfile(pwd, name + ".mtag.txt");

end

if opts.calcz % get Z-scores
    beta_col = find(hdr == opts.beta_name);
    se_col = find(hdr == opts.se_name);
    [~, name] = fileparts(opts.sumstats);
    newsumstats = makeWSLpath(fullfile(pwd, name + ".z.txt"));
    tt = tic;
    fprintf('calculating z-scores from beta/se...')
    if ~isfile(fullfile(pwd, name + ".z.txt"))
        runbash("awk -F'" + sep.sep + "' -v OFS='" + sep.sep + ...
            "' 'NR!=1 {$"+ beta_col + "=($" + beta_col + "/$" + se_col + ...
            ")}1' " + makeWSLpath(opts.sumstats) + " > " + newsumstats, ...
            tmpname, "wait", true);
    end
    fprintf('\b\bdone (%.2f sec)\n', toc(tt))    
    opts.sumstats = fullfile(pwd, name + ".z.txt");    
end

if opts.use_beta_as_z || opts.calcz
    opts.z_name = opts.beta_name;
    opts = rmfield(opts, ["beta_name", "se_name"]);
end

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
            cmd = join([cmd, "--" + fd(i)], " ");
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

        cmd = join([cmd, "--" + fd(i), tmpval], ' ');
        
    end
end
end