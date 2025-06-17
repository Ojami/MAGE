function gwasrunner(snps, opts)
% gwasrunner is the general function for genetic association analyses in UK
% Biobank (UKBB). It's mainly designed for few variants analysis in a
% phenome-wide manner; and to be efficient, all regression based analyses
% are performed on GPU. Currently it can perform, regression analysis for
% imputed UKBB variants; in case of rare-cases/or rare single variants,
% user can perform Firth's logistic regression. The function also can run
% whole-exome sequencing (WES) analysis for both single variant and
% rare-variants analyses (burden, SKAT and SKAT-O).
% INPUT:
%       -snps: list of input snps to check, in a struct with following
%        fields:
%               -imputed: mandatory fields of snp [list of snps] and chr
%                [chromosomes]; and optional minorAllele field (see
%                getbulkgeno function).
%                Example: snps.snp = "rs738409"; snps.chr = "22"
%               -WES: either gene [gene symbol, scalar] or variants [IDs
%                corresponding to UKBB WES BIM files] must be provided.
%                Example: snps.gene = "LPL";
%               -custom: a custom score (PRS, etc) can also be fed into the
%                function, with the same structure as output of getbulkgeno.
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, Jan 2021.
% NOTE: MATLAB R2019b or newer is needed!
% 
% [TODO]: add support for categorical covariates.
% 
% [EDITED]: rare variant tests were added (version 1.1), Feb 2021.
% [EDITED]: transformation functions for continuous traits were added, Apr 2021.
% [EDITED]: bug fixed in makedes subfunction with nan values in genotype vecotr, Apr 2021.
% [EDITED]: bug fixed with nan values in covariates, Apr 2021.
% [EDITED]: bugs fixed in callSKAT subfunction, Apr 2021
% [EDITED]: more than one gene can be queried for rare-variant analysis.
%           Example: snps.gene = ["LPL", "APOB"]; In this case, rare
%           variants on the queried genes are merged together.
% [EDITED]: collapse flag was implemented for rare variant analysis using
%           Firth logistic regression (version 1.1.1), Apr 2021
% [EDITED]: supports now parallele computing, see doc parfeval. Apr 2021
% [EDITED]: maf flag was added for MAF filtering. May 2021.
% [EDITED]: singlerare flag was added for association analysis of single 
%           rare variants from exome-sequencing data. Overrides 'collapse'
%           flag. Note that 'maf' flag will be used in this case instead of
%           'raremaf'. 14 May 2021.
% [EDITED]: 'error' flag was added for uninterrupted analysis if some
%           analyses failed. This currently works only for jobs sent to R 
%           i.e. SKAT and logist. If set to false, in case of error, empty
%           results will be printed to output. 18 May 2021.
% [EDITED]: 'maxTime' flag was added for R jobs. If R script runtime
%           exceeds this value, job will be interrupted incompletely, and
%           depending on 'error' flag, the function stops or ignores the
%           error. 18 May 2021.
% [EDITED]: now support parallel computations for single variant analysis.
%           also runModels function was modified for a better efficiency
%           and lower memory overhead. This can result in a more efficient
%           analysis using CPU-GPU hybrid. Note that number of threads 
%           should be selected carefully to avoid GPU memory issues. 21 May
%           2021 (version 1.2). 
% [EDITED]: "fitter" flag was added and now supports also "scoretest" which
%           is the implementation of REGENIE for score test (both linear
%           and logistic regressions). score test is much faster that "glm"
%           which uses MATLAB, however the default is still "glm".
%           Interaction terms and unadjusted p-values are not supported
%           with "scoretest" flag. 1 June 2021 (version 1.2.1)
% [EDITED]: for imputed variants, position range can be also included along
%           with variants. For instance input structure can be:
%                          query.snp = ["rs182549", "136608231-136616754"];
%                          query.chr = ["2", "2"];
%           position must correspond to genome version of bgi files, and
%           the format must be as pos1-pos2. 5 June 2021. 
% [EDITED]: online VEP annotation functionality was added for rare variants
%           filtering. It uses Ensembl REST API endpoints for VEP to
%           annotate and then filter variants based on different criteria.
%           Now only supports CADD-Phred, REVEL, SIFT, PolyPhen and
%           PrimateAI. 7 June 2021.
% [EDITED]: effPlotter was modified to add N (sample number per each
%           stratum) on Y-axis. 12 June 2021.
% [EDITED]: stratification analysis was modified to allow for any desired
%           trait (not limited to BMI anymore). 11 July 2021.
% [EDITED]: 'genesymbol' flag was added for gene/region-based analysis.
%           Default is true: gene names/symbols are provided to getWES
%           function; false: Ensembl (ENSG) identifiers in input query. 30
%           July 2021.
% [EDITED]: 'plot' flag now calls 'rarePlotter' for single rare variant
%           analysis ('singlerare' flag must be true) and 'effPlotter'
%           otherwise. 'raremaf' option can also be set, but is only
%           effective for 'rarePlotter' behvaior and not association
%           analysis (see notes above on 'singlerare' option). 3
%           August 2021.
% [EDITED]: 'vep' now also supports 'clinvar' flag for filtering based on
%           pathogenic variants in ClinVar database. 17 August 2021.
% [EDITED]: 'vep' filtering with 'dbnsfp' was added. 'dbnsfp' uses local
%           files parsed with dbNSFPparser function. Note that with
%           'dbnsfp', some filters may not be available,see callVEP
%           function for more details. 18 August 2021.
% [EDITED]: a bug fixed in call to getbulkgeno function. A bug fixed in
%           appending vep table to wesinfo table. 25 August 2021.
% [EDITED]: missingness filter from rare-variant analysis in
%           'prunewescalls' function was commented since now the collapsing
%           approaches both in SKAT and gwasrunner do not remove missing
%           samples. This filter origianally was used when setting the
%           missingness impute in SKAT package to false in order to keep
%           sample exclusion minimum. 27 August 2021.
% [EDITED]: a bug was fixed with opts.ld. 14 September 2021.
% [EDITED]: a bug was fixed with interaction term included in parallel
%           mode. 18 September 2021. 
% [EDITED]: 'plot' flag now calls interactionPlotter when 'interaction'
%           option is used. interactionPlotter function is still under
%           development. 29 September 2021.
% [EDITED]: 'barplot' flag was added for stratification analyses
%           (i.e. 'stratify' must be set). In addition, 'barfunc' can be
%           also set for each trait (can be vector with either median or
%           mean), which determines the function to be applied to the trait
%           in each stratum (e.g. for ALT one can use mean, and for
%           triglycerides median can be used). 30 September 2021.
% [EDITED]: 'defaultcovar' flag was added (default: true) to include
%           default covariates for UKBB analysis: first 10 genomic PC and
%           array batch. 15 October 2021.
% [EDITED]: 'getPheno' is now the default function for loading traits when
%           'trait' option is left empty. 15 October 2021.
% [EDITED]: a bug was fixed. 16 October 2021.
% [EDITED]: a bug fixed with unadjusted model fitting and missingness. 16
%           October 2021.
% [EDITED]: 'missingness' option was added for rare variant analysis. Any
%           rare variant with a missingness rate above this cutoff will be
%           excluded from the analyses. 'singlerare' flag overrides this
%           option. 25 October 2021.
% [EDITED]: 'callVEP' now also supports LOFTEE filtering. To diffrentiate
%           between 'loftee' flag from 'getWES', this option was renamed to
%           'veploftee'. 30 October 2021.
% [EDITED]: a bug was fixed in runModels. 17 December 2021.
% [EDITED]: 'ciround' option was added for rounding 95% CI. Similarly,
%           'betaround' was added for P, beta and SE. 3 January 2022.
% [EDITED]: 'lof' flag can be now combined with custom geno structures for
%           rare variant analysis. 3 January 2022.
% [EDITED]: some bugs were fixed. 10 January 2022.
% [EDITED]: a bug was fixed. 2 Februrary 2022.
% [EDITED]: getbulkgeno function can now accept 'eid' option, subsetting
%           directly bgen file and keeping only samples using 'qc', 'qcinc'
%           options. This imporves the efficiency of reading and memory
%           usage (less samples, in case of few qcinc ids). 17 Februrary
%           2022.
% [EDITED]: exclusion criteria can now be applied to covariates with
%           'exeid' field as well. 14 April 2022.
% [EDITED]: some bugs were fixed. 15 April 2022.
% [EDITED]: 'condition' option was included for conditioning on a list of
%           variants. These variants will be included in covariates and
%           removed from the main analyses. Now only works with 'regenie'
%           in gene-based tests. Example: if 'gene' field is
%           "gene1,gene2,gene3", and 'condition' is ["gene2", "gene3"],
%           then all variants on genes 2 and 3 will be included in
%           covariates. Alternatively, if 'condition' is "gene2", genes 1
%           and 3 will be merged (in practice, considered as a single
%           gene), and only variants on gene 2 will be condition on. 2 May
%           2022.
% [EDITED]: 'catCovar' option was added for categorical covariates. 22 May
%           2022.
% [EDITED]: 'snps' struct can have a 'tag' field to be used instead of
%           variant ids in the output. 7 June 2022.
% [EDITED]: now prints full summary stats when using 'interaction' option.
%           The function appends the rows similar to REGENIE's output and
%           adds a TEST column with ADD/REC/MOD, SNPxVAR and VAR rows,
%           where VAR is the 'interaction' argument. 8 August 2022. 
% [EDITED]: more 'regenie' otions were added to support gene-based tests
%           based on pre-prepared annotation/set/mask files (geneWrapper
%           output). 7 September 2022.
% [EDITED]: 'saveGeno' option was added to save pruned genotype calls for
%           gene-based analysis. Default is "no", options are "yes", "no",
%           and "only". "only" saves only genotype calls and skip the
%           rests. 'regenie' flag must be false. 17 October 2022.
% [EDITED]: 'qc' now accepts -1 for studies where their qced IDs have not
%           been defined or included in 'getQCEID' function. In such cases,
%           'qcinc' should be used and this field replaces 'qc', meaning
%           inclusion eids and qced eids should be the same. 09 February
%           2023.
% 13MAY2023: 'covarTab' option was implemeneted for additional covariates
%           to those already provided in 'covar'. 'covarTab' must be a
%           table with and 'eid' column. 
% 20MAY2023: negative values are no longer considered as missingness in
%           the 'runModels' subfunction.
% 
% 12MAR2024: 'force_robust' and 'rare_mac' options were added for
%            interaction analysis. These two options make it possible to
%            perform interaction using robust SE, HLM or mixture of both.
%            By default, 'force_robust' is false, meaning that for variants
%            > rare_mac, robust SE is used, while HLM is used for rare
%            variants (< rare_mac). If 'force_robust' is true, then HLM is
%            ignored (ergo 'rare_mac') and only robust SE (HC3 sandwich
%            estimator) is used. If one opts for only HLM, then
%            'force_robust' should be false (default) and 'rare_mac' should
%            be set to a fairly large value. Both work only when 'regenie'
%            is used.
% 
% 14MAR2024: 'debug' flag now plots diagnostic plots for GLM.
% 29APR2024: 'phenoFile' and 'covarFile' options added for REGENIE step 2
%            analysese. In this case, gwasrunner does not read
%            covariates/phenotypes anymore, and directly uses those files.
%            Note that at the moment, both files should be present, and
%            it's not possible to use either and left the other one empty.
% 09JULY2024: when 'phenoFile' and 'covarFile' for REGENIE are provided,
%             and 'qcinc' is non-empty, it subsets the individuals
%             automatically (this was a bug).
% 17NOV2024: if 'pred' is used without 'regenie' flag, the function
%            internally extracts step 1 (from REGENIE) genetic vector and 
%            adds it as an extra covariate (also for survival analysis).
% 26NOV2024: 'catCovarCell' option was added which accepts the same format
%            covarFile and phenoFile (i.e. cell of string vectors). This
%            is useful when categorical variables are different per each
%            pheno/covarfiles. This option is used only when
%            covarFile/phenoFile are used and in that case, overrides,
%            'catCovar' option.If there are multiple categorical
%             covariates for each pheno, they should be "," separated.
% 31MAY2025: 'scoretest' now supports 'firth' regression with exact
%             implementation of logistf based on profile penalized log
%             likelihood.

arguments
    snps {mustBeA(snps, 'struct')}
    opts.gpu (1,1) logical = true % use GPU for computations
    opts.bgenhome {mustBeFolder}
    opts.parallel (1,1) logical = false % use parallel computations for single variant analyses
    opts.workers (1,1) double = NaN % number of workers; default is number of physical cores available.
    opts.chunk (1, 1) double = 200 % max number of snps per calculations for "scoretest" fitter under GPU calculations. default for GPU computations is 200, and for CPU is 10-folde higher (i.e. 2000) because of memory limitations of GPU. In this case parallel and workers are not used, instead MATLAB built-in multithreading does the job.
    opts.infoscore (1,1) double {mustBeInRange(opts.infoscore, 1e-6, 1)} = 0.75 % INFO score cut-off: variants < INFO score are removed.
    opts.maf (1,1) double {mustBeInRange(opts.maf, 0, 0.5)} = 0 % MAF filter: removes MAF < this cut-off
    opts.dosage (1,1) logical = true % dosage or hard genotype calls
    opts.model {mustBeMember(opts.model, ["add", "dom", "rec"])} = "add"
    opts.ld (1,1) logical = false % computes r2
    opts.qc (1,1) double {mustBeMember(opts.qc, [-1, 0:6])} = 1 % see getQCEID function for more details. set -1 if custom qced IDs will be used in qcinc; in this case, 'qcinc' replaces 'qc' field, this is useful for non UKBB cohorts
    opts.qcinc (:, 1) {mustBeNonmissing} = [] % custom eids to be included; in this case intersect of qc and qcinc will be used for analyses.
    opts.exclude {mustBeFile} % exclusion file of variant IDs to remove from the analyses

    % options for stratification
    opts.stratify {mustBeTextScalar} = "" % if provided (default: don't do stratification), analysis is performed for strata of this phenotype
    opts.bins (1, :) double {mustBeGreaterThan(opts.bins, 0)} = [18.5 25 30 35]
    opts.stratifyby {mustBeMember(opts.stratifyby, ["mean", "median", "bins"])} = "bins" % stratification is performed based on either of mean or median. If left empty (if strata is queried), values in bin are used instead.
    opts.quantile (1, :) double = NaN % quantiles for stratification, e.g. for quartiles: [0.25, 0.5, 0.75]. If NaN (default), values in bins or 'strataby' are used. 
    opts.removeStrat (1, 1) logical = false % do not use the 'stratify' trait for analyses (i.e. remove it from covariates)
    
    opts.transform {mustBeMember(opts.transform, ["none", "log", "log1p", "log10", "irnt"])} = "irnt" % transformation function for continuous variables.
    opts.log (1,1) logical = false
    opts.phenopath {mustBeFolder} = fullfile(fileparts(which("gwasrunner.m")), "UKBB.Toolbox", "UKB_PHENO")
    opts.trait {mustBeText} = "" % if left empty, getPheno will fetch traits
    opts.plot (1,1) logical = false % plots the summary stats depending on the analysis scheme under test.
    opts.barplot (1,1) logical = false % plots grouped bars for stratification analysis 
    opts.barfunc {mustBeMember(opts.barfunc, ["mean", "median"]), mustBeVector} = "mean" % function to be applied to response in each stratum (on y-axis): can be a vector for each trait
    
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "jpg"; % format of output plots
    opts.resolution (1, 1) double = 400 % resolution of output plot
    opts.saveplot (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    opts.output {mustBeTextScalar} = "" % name of output file
    opts.interaction {mustBeTextScalar} = "" % adds SNP*trait interaction to the regression model. It's the name of MAT file and not the tag string.
    opts.covar {mustBeText, mustBeVector} = ["Age", "Sex", "BMI"]; % list of covariates, it automatically adds PC1-10 and array batch.
    opts.covarTab {mustBeA(opts.covarTab, "table")} % additional covariates to be added (a table with covariates + 'eid' column).
    opts.catCovar {mustBeText, mustBeVector} % list of categorical covariates
    opts.defaultcovar (1,1) logical = true % to automatically add PC1-10 and array batch to covariates.
    opts.wes (1,1) logical = false % for WES analysis. Array batch is excluded from covariates in this case.
    opts.sheet {mustBeTextScalar} = "Overall" % data-sheet name for results
    opts.firth (1,1) logical = false % performs firth regression via logistf R package.
    opts.flic (1,1) logical = false % see logistf package help for this, or https://github.com/georgheinze/flicflac
    opts.flac (1,1) logical = false % see logistf package help for this, or https://github.com/georgheinze/flicflac
    opts.des (1,1) logical = false % generates descriptive information, including overall counts/mean/median and stratified by genotype.
    opts.desonly (1,1) logical = false % doesn't perform association analysis, and only generates descriptive information, including overall counts/mean/median and stratified by genotype. Overrides opts.des flag.
    opts.cleanOutput (1,1) logical = true % keeps only one header and removes extra spaces between rows within same data sheet; useful when several traits for the same variant are being tested.
    opts.warning (1,1) logical = false % trun off warnings
    opts.error (1,1) logical = true; % if set to false, doesn't interrupt if fails running R functions, instead returns/saves NaN to the output. Useful when doesn't want unwanted errors interfere with the analyses.
    opts.maxTime (1,1) double = inf; % max time allowed to perform calculations in R (i.e. SKAT and logistf package functions)
    
    % getWES function specific inputs. For info see the doc within getWES
    % function
    opts.weshome {mustBeFolder} 
    opts.bed {mustBeTextScalar} = "" 
    opts.plinkdir {mustBeTextScalar} = "" % either PLINK 1.9 or PLINK 2 should be there 
    opts.annhome {mustBeFolder} % see getWES
    opts.annfiles {mustBeNonzeroLengthText} 
    opts.genesymbol (1, 1) logical = true; % true: gene symbol, false: Ensembl (ENSG) identifier.
    opts.genePos (:, 2) double = [NaN, NaN]; % gene boundary. Works only if 'variants' set to "ALL" (all variants on the gene)
    opts.lof (1,1) logical = false
    opts.loftee (1,1) logical = false % only subset those LoF with a true LOFTEE flag
    opts.verbose (1,1) logical = true % prints related info
    opts.debug (1,1) logical = false % detailed version of 'verbose' flag
    opts.ignore_anno (1,1) logical = false % ignores reading annotation file and only return 'variants'. Should be used only if 'variants' field is used with pre-annotated variants (e.g. with 'regenie' flag)
    
    % current approach: haploinsufficient model (presence/absence):
    % https://www.nature.com/articles/s41576-019-0177-4#Sec1 (also here:
    % https://www.medrxiv.org/content/10.1101/2021.07.12.21260400v1.full.pdf)
    opts.collapse (1,1) logical = false % doesn't perform rare-variant tests; instead, generates a collapsing score (first approach mentioned here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3277891/).
    opts.singlerare (1,1) logical = false % doesn't perform rare-variant tests; instead, generates summary statistics for each rare variant in the input query. Best if firth flag is set to true. Overrides collapse flag.
    opts.raremaf (1,1) double {mustBeInRange(opts.raremaf, 0, 0.5)} = 0.01 % MAF cut-off for rare variant analysis. Variants >= this MAF are excluded from analyses. If any value > 0.01 is inserted, rare+common variant analysis will be performed instead (see SKAT documentation).
    opts.missingness (1,1) double {mustBeInRange(opts.missingness, 0, 1)} = 0.15 % rare variants with missingness rate above this cut-off are excluded from analyses.
    opts.saveGeno {mustBeTextScalar, mustBeMember(opts.saveGeno, ["no", "yes", "only"])} = "no" % save pruned genotype calls as a mat file for the query gene. 'regenie' flag must be false. "only": skips the association tests and only saves the genotype calls. 

    % callVEP inputs for rare variant liter
    opts.vep (1,1) logical = false % uses VEP API (Ensembl REST) to further filter variants
    opts.cadd (1,1) double = 0; % filter variants with CADD-phred above this
    opts.revel (1,1) double = 0; % filter variants with REVEL above this
    opts.sift {mustBeMember(opts.sift, ["strict", "relaxed", "none"])} = "none";
    opts.polyphen {mustBeMember(opts.polyphen, ["strict", "relaxed", "none"])} = "none"; % probably_damaging is "D" and possibly_damaging "P" in dbNSFP 
    opts.primateai (1,1) logical = false; % if primateAI == "D"
    opts.spliceai (1,1) double = 0; % filter variants with spliceAI value above this. Hint: 0.2 (high recall) | 0.5 (recommended) | 0.8 (high precision)
    opts.clinvar (1,1) logical = false; % if true, matches ClinVar flags of "Likely_pathogenic" and "pathogenic"
    opts.veploftee (1,1) logical = false; % keeps HC variants
    opts.rule {mustBeMember(opts.rule, ["and", "or"])} = "or"; % how to merge filters
    opts.method {mustBeMember(opts.method, ["vep", "dbnsfp"])} = "dbnsfp"
    opts.dbnsfpHome {mustBeTextScalar}  % dbNSFP parsed annotations directory. Files should be available for each chromosome separately. only used when method is set to "dbnsfp"
    opts.dbnsfpFile {mustBeTextScalar} = "" % dbNSFP file pattern, this is helpful if there are non dbNSFP files in 'dbnsfpHome' dir. This can be part of the file name for all chromosomes, e.g. "dbNSFP4.2a_variant.chr" matches dbNSFP4.2a_variant.chr1.txt to dbNSFP4.2a_variant.chr22.txt. If left empty (default), only dbNSFP parsed files (txt) are present in the dir.
    
    % developer options
    opts.condition {mustBeText, mustBeVector} % list of variants to condition on (see doc for more details)
    opts.fitter {mustBeMember(opts.fitter, ["glm", "scoretest"])} = "glm" % "glm" uses fitglm function, while "scoretest" uses a faster approach of score test statistics for both logistic and linear regressions.
    opts.regenie (1,1) logical = false % use regenie software for regression analyses % UNDER DEVELOPMENT
    opts.ciround (1,1) {mustBeGreaterThan(opts.ciround, 1)} = 2 % round off for 95% CI
    opts.betaround (1,1) {mustBeGreaterThan(opts.betaround, 1)} = 3 % round off for p-value, beta and SE

    % regenie inputs
    opts.pred {mustBeText} % if left empty, typical linear/logisitc reg is performed
    opts.phenoFile {mustBeFile} % the same size as covarFile and pred
    opts.covarFile {mustBeFile} % the same size as phenoFile and pred
    opts.aaf_bins double {mustBeNonNan, mustBeNonnegative, mustBeVector} = [0.01, 0.005] % comma-separated list of AAF upper bounds to use when building masks [default is a single cutoff of 1%]
    opts.build_mask {mustBeTextScalar, mustBeMember(opts.build_mask, ["max", "sum", "comphet"])} = "max" % build masks using the maximum number of ALT alleles across sites ('max'; the default), or the sum of ALT alleles ('sum'), or thresholding the sum to 2 ('comphet')
    opts.skat_params (1,2) double {mustBeNonnegative, mustBeNonNan} = [1,25] % skat parameters
    opts.vc_maxAAF (1,1) double {mustBeInRange(opts.vc_maxAAF, 1e-12, 1)} = 0.01 % AAF upper bound to use for SKAT/ACAT-type tests [default is 100%]
    opts.anno_file {mustBeFile} % File with variant annotations for each set
    opts.set_list {mustBeFile} % File listing variant sets
    opts.mask_def {mustBeFile} % File with mask definitions using the annotations defined in --anno-file
    opts.rare_mac (1,1) double = 1000 % minor allele count (MAC) threshold below which to use HLM method for QTs [default is 1000]
    opts.force_robust (1,1) logical = false % use robust SE instead of HLM for rare variant GxE test with quantitative traits
    opts.write_mask (1,1) logical = false % write mask to PLINK bed format (does not work when building masks with 'sum')

    % survival analysis
    opts.checkDiagnostics (1,1) logical = false
    opts.logRank (1,1) logical = false
    opts.FGR {mustBeTextScalar, mustBeMember(opts.FGR, ["none", "crr", "fastcrr"])} = "none" % Fine-Gray for competing risks
    opts.prevalent (1,1) logical = false % exclude prevalent cases
    opts.timeOrigin {mustBeMember(opts.timeOrigin, ["TTE", "AOO"])} = "TTE"
    opts.LT (1,1) logical = false % no left-truncation
    opts.logRankFun {mustBeTextScalar, mustBeMember(opts.logRankFun, ["event", "cumhaz", "pct", "null"])} = "null" 
    opts.landmark (1,1) double = 0 % landmark analysis, individuals with a non-zero censoring status (including competing risks) with a tt <= landmark (yr) are excluded. 

    %@26NOV2024
    opts.catCovarCell {mustBeText, mustBeVector} % the same size as phenoFile and pred. overrides 'catCovar'. Only works with REGENIE for now.
end
% prepare inputs ----------------------------------------------------------
if ~opts.warning
    warning('off')
end

if opts.verbose; printLogo; end
opts.cleanOutputcheck = []; % for opts.cleanOutput
opts.sheet = string(opts.sheet);
if opts.desonly
    opts.des = true; % override opts.des flag
end

if opts.singlerare
    opts.collapse = false; % override opts.collapse flag
end

if opts.log 
    clc
    diary gwasrunner.log
end

if opts.gpu
    % check GPU memory
    g = gpuDevice(1);
    reset(g);
else
    opts.chunk = opts.chunk.*10; % only for "scoretest" fitter
end

% check parallel pool
if opts.parallel && isempty(gcp('nocreate'))
    if isnan(opts.workers)
        % number of physical cores usually gives the best efficiency for
        % regression analyses. Although this has not been thoroughly
        % tested.
        maxNumCompThreads('automatic'); % reset (only if was changed in current session)
        opts.workers = maxNumCompThreads; 
    end
    parc = parcluster("Processes");
    parc.NumWorkers = opts.workers;
    parpool(parc); % should I delete in the end? 
end

% getting samples ---------------------------------------------------------
if opts.qc >= 0 % otherwise qcinc should be used as custom qced ids
    qc_eid_out = getQCEID(opts.qc, opts.verbose, opts.phenopath); 
else
    qc_eid_out = opts.qcinc;
    if ~opts.regenie, opts.qcinc = []; end
end
if opts.verbose; fprintf('\nsample N = %d\n', numel(qc_eid_out)); end

if ~isempty(opts.qcinc) % custom samples (eid)
    qc_eid_out = intersect(opts.qcinc, qc_eid_out);
    if opts.verbose
        fprintf('sample N after applying custom samples = %d\n', numel(qc_eid_out))
    end
    if ~opts.regenie, opts.qcinc = []; end
end

% variant inclusion/exclusion criteria ------------------------------------
if isfield(opts, 'exclude')
    opts.exclude = readlines(opts.exclude);
    opts.exclude(opts.exclude == "") = [];
else
    opts.exclude = [];
end

% 'condition' option ------------------------------------------------------
if isfield(opts, 'condition')
    opts.condition = string(opts.condition);
    emptyidx = opts.condition == "" | ismissing(opts.condition);
    opts.condition(emptyidx) = [];
else
    opts.condition = [];
end

if opts.regenie
    opts.fitter = "regenie"; 
    checkPhenoCovarFile = isfield(opts, ["covarFile", "phenoFile"]);
    assert(all(checkPhenoCovarFile) | all(~checkPhenoCovarFile), "both phenoFile and covarFile must be present!")
end

% 'output' ----------------------------------------------------------------
if string(opts.output) == ""
    opts.output = "UKBBoutput." +...
        string(datetime('today', 'Format', 'dd.MM.yyyy'));
elseif ~isfolder(fileparts(opts.output))
        mkdir(fileparts(opts.output))
end

% get calls ---------------------------------------------------------------
if ~isfield(snps, 'minorAllele')
    snps.minorAllele = "";
end

if any(isfield(snps, {'gene', 'variants'})) % WES query ...................
    opts.wesin = true;
    if isfield(snps, 'gene') && isstruct(snps.gene)% check if raw WES struct is provided
        opts.customgeno = true;
    else
        opts.customgeno = false;
        if ~isfield(snps, 'gene')
            snps.gene = "";
            
            % apply variant inclusion/exclusion criteria
            if ~isempty(opts.exclude)
                idx = ismember(snps.variants, opts.exclude);
                snps.variants(idx) = [];
            else
                opts.exclude = "";
            end

        elseif ~isfield(snps, 'variants')
            snps.variants = [];
        end
        snps.gene = string(snps.gene);
    end
    opts.wes = true; 
    
    if opts.verbose
        disp("gene-/region-based analysis.")
        tic; 
    end
    
    if opts.customgeno
        geno = snps.gene.geno; snps.gene.geno = [];
        wesinfo = snps.gene.wesinfo; snps.gene.wesinfo = [];
        if ~iscell(geno); geno = {geno}; end
        if ~iscell(wesinfo); wesinfo = {wesinfo}; end

        % check if getWES loftee flag is requested and throw a warning
        if ~opts.vep && opts.loftee
            if opts.verbose
                fprintf("WARNING:'loftee' flag are ignored with custom geno structure!\n")
                fprintf("         use 'vep' option with veploftee flag\n")
            end
        end
    else
        
        if ~isempty(snps.variants) % reformat variants (to diffrentiate between variant sets)
            checkvars = snps.variants;
            if isstring(checkvars), checkvars = cellstr(checkvars); end
            if iscellstr(checkvars), checkvars = {checkvars}; end
            snps.variants = checkvars;
        else
            snps.variants = "";
        end

        if ~opts.regenie && (numel(snps.gene) > 1 || (iscellstr(snps.variants) && numel(snps.variants) > 1))
            error('gwasrunner:multipleGenes', 'multiple genes or set of variants are not supported without regenie flag!')
        end
        
        % if 'regenie' is true and anno/mask/set files are used, no need to
        % parse/fetch variants for the query gene. Only find the bed file
        % associated with the gene, and extract anno/set files based on the
        % query gene
        if ~opts.regenie && isfield(opts, 'set_list')
            if opts.verbose
                fprintf('WARNING: set_list will be overlooked when regenie is not used!\n')
            end
            opts.set_list = ""; opts.anno_file = ""; opts.mask_def = "";
        elseif ~isfield(opts, 'set_list')
            opts.set_list = ""; opts.anno_file = ""; opts.mask_def = "";
        end
        
        
        if isempty(opts.exclude), opts.exclude = ""; end
        outdir = fileparts(string(opts.output));
        if outdir == "", outdir = pwd; end
        [geno, wesinfo] = getWES('gene', snps.gene, 'variants', snps.variants, ...
            'home', opts.weshome, 'bed', opts.bed, 'plinkdir', opts.plinkdir,...
            'annfiles', opts.annfiles, ...
            'annhome', opts.annhome, 'lof', opts.lof, 'loftee', opts.loftee,...
            'verbose', opts.verbose, 'genesymbol', opts.genesymbol,...
            'genePos', opts.genePos, 'exclude', opts.exclude, ...
            'regenie', opts.regenie, 'outdir', outdir, ...
            'set_list', opts.set_list, 'anno_file', opts.anno_file, ...
            'mask_def', opts.mask_def, 'ignore_anno', opts.ignore_anno, ...
            'parallel', opts.parallel);
        if isempty(geno)
            if opts.verbose, fprintf("gene was not found!\n"); end
            return
        end

        if numel(geno) > 1 && opts.regenie
            % multiple genes without set/mask/anno files. Concatenate
            % different genes
            gtmp = cellfun(@(x)x.bim, geno);
            gfis = string(fieldnames(gtmp));
            if isrow(gtmp), gtmp = gtmp'; end
            geno = struct;
            for kk = 1:numel(gfis)
                geno.bim.(gfis(kk)) = {gtmp.(gfis(kk))};

                % keep idx of each gene (for call_regenie_in)
                if kk == 1
                    gidxtmp = cumsum(cellfun(@numel, geno.bim.(gfis(kk))));
                    geno.bim.idx = [[1, gidxtmp(1:end-1) + 1]', gidxtmp'];
                end
                geno.bim.(gfis(kk)) = vertcat(geno.bim.(gfis(kk)){:});        
            end
            wesinfo = vertcat(wesinfo{:});
        else
            geno = geno{1}; wesinfo = wesinfo{1};
        end
    end
    if opts.verbose; toc; end

    if opts.customgeno
        if numel(geno) > 1
            error("should be fixed!")
            % wesinfo = vertcat(wesinfo{:});
            % if ~istable(wesinfo), wesinfo = vertcat(wesinfo{:}); end
        else
            geno = geno{1}; wesinfo = wesinfo{1};
        end
    end
    
    % @03/01/2022:
    % check if 'lof' flag was requested with custom geno structure(s). In
    % this case only keep LoF variants using wesinfo table
    if opts.lof && opts.customgeno && ~opts.regenie
        keepIdx = contains(wesinfo.consequence, ["stop_gained"; "splice_donor_variant"; ...
            "frameshift_variant"; "splice_acceptor_variant"]);
        tmpLoFids = wesinfo.id(keepIdx);
        wesinfo(~keepIdx, :) = [];
        idx = find(~ismember(geno.bim.snp, tmpLoFids));
        flds = fieldnames(geno.bim);
        geno.bed(:, [2.*idx-1; 2.*idx]) = [];
        for j = 1:numel(flds)
            geno.bim.(flds{j})(idx) = [];
        end
    end
    
    % @08/03/2021: a aa.change was misspelled in getWES.
    wesinfo.Properties.VariableNames(ismember(wesinfo.Properties.VariableNames, 'aa.chage'))  = {'aa.change'};
    
    % VEP annotation 
    if opts.vep
        if opts.verbose; fprintf('\nVEP filtering\n'); end
        if opts.method == "dbnsfp" % don't filter LoF variants 
             lof_filters = ["stop_gained"; "splice_donor_variant"; ...
                    "frameshift_variant"; "splice_acceptor_variant"];
             lof_ids = wesinfo.id(contains(wesinfo.consequence, lof_filters));
             lof_idx = ismember(geno.bim.snp, lof_ids);
             bimfields = setdiff(fieldnames(geno.bim), "idx");
             bim = struct;
             if any(~lof_idx) % check if there are variants remained to be filtered
                 for k = 1:numel(bimfields)
                     bim.(bimfields{k}) = geno.bim.(bimfields{k})(~lof_idx);
                 end
             end
        else % 'vep'
            bim = geno.bim;
        end

        if ~isempty(setdiff(fieldnames(bim), "idx"))
        
            vep = callVEP(bim, 'filter', true, 'cadd', opts.cadd, ...
                'revel', opts.revel, 'sift', opts.sift, 'loftee', opts.veploftee,...   
                'polyphen', opts.polyphen, 'primateai', opts.primateai, ...
                'rule', opts.rule, 'spliceai', opts.spliceai,...
                'clinvar', opts.clinvar, 'method', opts.method,...
                'dbnsfpHome', opts.dbnsfpHome, ...
                'dbnsfpFile', opts.dbnsfpFile, 'verbose', opts.verbose);
            % keep only those remained after VEP filtering
            if opts.method == "dbnsfp"
                idx = find(~ismember(geno.bim.snp, [vep.id; lof_ids]));
            else
                idx = find(~ismember(geno.bim.snp, vep.id));
            end
            flds = setdiff(fieldnames(geno.bim), "idx");

            if ~opts.regenie, geno.bed(:, [2.*idx-1; 2.*idx]) = []; end
            for j = 1:numel(flds)
                geno.bim.(flds{j})(idx) = [];
            end
            
            % append VEP annotations to wesinfo table
            [idx1, idx2] = ismember(wesinfo.id, vep.id); idx2(idx2<1) = [];
            if opts.method == "dbnsfp" 
                vep.Properties.VariableNames{2} = 'rsid';
                appendvarnames = vep.Properties.VariableNames(2:end);
            else
                appendvarnames = vep.Properties.VariableNames(5:end-2); % don't keep overlapping names for 'vep'
            end
            
            if ~isempty(idx2)
                for ct = 1:numel(appendvarnames)
                    wesinfo.(appendvarnames{ct})(idx1) = vep.(appendvarnames{ct})(idx2);
                end
            end
            % those that are absent from wesinfo table
            vep(ismember(vep.id, wesinfo.id), :) = [];
            if ~isempty(vep)
                vep.LoF = []; vep.LoFMissense = [];
                vep.Properties.VariableNames{2} = 'rsid';
                if opts.method == "dbnsfp" 
                    [vep.pos, vep.gene, vep.ENSG, vep.consequence, vep.('aa.change')] = deal(strings(size(vep, 1), 1));
                    vep = movevars(vep, {'pos', 'gene', 'ENSG', 'consequence', 'aa.change'}, 'After', 'rsid');
                else
                    [vep.pos, vep.ENSG, vep.('aa.change')] = deal(strings(size(vep, 1), 1));
                    vep = movevars(vep, {'pos', 'gene', 'ENSG'}, 'After', 'rsid');
                end
                wesinfo = [wesinfo; vep];
            end
            if isnumeric(wesinfo.pos)
                wesinfo.pos(wesinfo.pos < 1) = nan;
            end
        end
            
        if opts.verbose; fprintf('\n'); end
    end
    
    % remove those ids absent from geno struct
    nonBimIdx = ~ismember(wesinfo.id, geno.bim.snp);
    wesinfo(nonBimIdx, :) = [];
    
else % imputed query.......................................................
    if opts.verbose, disp("single variant analysis."); end
    opts.wesin = false;
    if ~isfield(snps, 'snp') % variant calls input. It must conform to output of getbulkgeno
        opts.customgeno = true;
        geno = snps;
        geno = rmfield(geno, 'minorAllele');

        % check 'tag's
        fis = fieldnames(geno);
        for k = 1:numel(fis)
            if isfield(geno.(fis{k}), "tag") && all(geno.(fis{k}).tag ~= "")
                geno.(fis{k}).snp = geno.(fis{k}).tag;
            end
        end
    else

        % apply variant inclusion/exclusion criteria
        if ~isempty(opts.exclude)
            idx = ismember(snps.snp, opts.exclude);
            snps.snp(idx) = []; snps.chr(idx) = [];
            if isempty(snps.snp)
                fprintf('gwasrunner:no variant was remained after applying variant exclusion criteria!\n')
                return
            end
        end
        
        if ~isfield(snps, "tag"), snps.tag = ""; end
        if opts.regenie
            [~, geno] = getbulkgeno(strtrim(snps.snp), strtrim(snps.chr),...
                'dosage', opts.dosage, 'minorAllele', snps.minorAllele,...
                'model', opts.model, 'home', opts.bgenhome, ...
                'optionsOnly', true, 'infoscore', opts.infoscore, ...
                'verbose', opts.verbose, "tag", strtrim(snps.tag));

            if any(contains(snps.snp, digitsPattern+"-"+digitsPattern))
                % input is a region, fetch snps within the region
                bgid = bgireader(snps, bgenhome=opts.bgenhome, verbose=opts.verbose);
                if iscell(bgid), bgid = bgid{1}; end
                snps.chr = string({bgid.chr})';
                snps.snp = string({bgid.snp})';
            end
            opts.regenied.achr = snps.chr;
            opts.regenied.asnp = snps.snp; % no guarantee these snps are present in bgen files
            opts.regenied.atag = snps.tag;
            opts.regenied.abgenfile = geno.bgenfile;
            opts.regenied.abgenchr = geno.chrom;
        else
            
            geno = getbulkgeno(strtrim(snps.snp), strtrim(snps.chr),...
                'dosage', opts.dosage, 'minorAllele', snps.minorAllele,...
                'model', opts.model, 'home', opts.bgenhome,...
                'parallel', opts.parallel, 'infoscore', opts.infoscore, ...
                'eid', qc_eid_out, "tag", strtrim(snps.tag), 'verbose', opts.verbose);
        end
        opts.customgeno = false;
        try
            if all(snps.tag ~= "")
                fis = fieldnames(geno);
                for k = 1:numel(fis)
                    geno.(fis{k}).snp = geno.(fis{k}).tag;
                end
            end
        catch 
            % regenie 
        end
    end
end
% end of getting calls ----------------------------------------------------
if opts.regenie && opts.wesin
    opts.wesinfo = wesinfo; clear wesinfo
elseif opts.wesin% prune WES calls: remove missingness and keep only qced ids
    [geno, opts.wesinfo, qc_eid_out] = prunewescalls(geno, wesinfo, qc_eid_out, opts);
    
    % write genotype calls to a file
    if any(opts.saveGeno == ["yes", "only"])
        wesinfo = opts.wesinfo;
        save(opts.output + ".geno.mat", 'geno', 'wesinfo')
        if opts.saveGeno == "only", return; end
    end
    
    if isempty(geno.bed) % all missing!
        fprintf('ERROR:no variant calls remained after pruning!\n')
        return
    end
    clear wesinfo
    if opts.collapse || opts.singlerare % skip burden tests 
        % best guess impute: note that REGENIE constructs burden masks a bit
        % differently: https://github.com/rgcgithub/regenie/issues/259
        geno.bed = double(geno.bed);
        idx = isnan(geno.bed) | (geno.bed > 2) | (geno.bed < 0); % missingness
        geno.bed(idx) = nan;
        maf = mean(geno.bed, 1, 'omitnan')./2;
        p = size(geno.bed, 2);
        for k = 1:p, geno.bed(idx(:, k), k) = round(maf(k).*2); end

        if opts.collapse

            % use REGENIE build_mask argument
            if opts.build_mask == "max"
                geno.bed = max(geno.bed, [], 2);
            elseif opts.build_mask == "sum"
                geno.bed = sum(geno.bed, 2);
            elseif opts.build_mask == "comphet"
                 geno.bed = sum(geno.bed, 2);
                 geno.bed(geno.bed > 2) = 2;
            else
                error("gwasrunner:unknown build_mask %s", opts.build_mask)
            end
        end
        geno.eid = geno.fam;
        geno = rmfield(geno, 'fam');
        if opts.collapse
            geno.chr = string(geno.bim.chr(1));
            geno.snp = join(unique(opts.wesinfo.gene), ';');
            geno.a1 = ""; geno.a2 = "";
            geno.pos = min(geno.bim.pos) + "-" + max(geno.bim.pos);
        else
            geno.chr = geno.bim.chr;
            geno.snp = geno.bim.snp;
            geno.a1 = geno.bim.a1; geno.a2 = geno.bim.a2;
            geno.pos = geno.bim.pos;
        end
        geno = rmfield(geno, 'bim');
        geno.("chr"+geno.chr(1)) = geno;
        exfields = fieldnames(geno);
        exfields = setdiff(exfields, "chr"+geno.chr);
        geno = rmfield(geno, exfields);
        opts.wesin = false; opts.wes = true;
    end 
end
[skat, collapseres] = deal(cell(numel(opts.trait), 1)); 

if ~opts.wesin && ~opts.regenie % MAF filtering for single variant/collapsing analysis
    chrFields = fieldnames(geno);
    for i = 1:numel(chrFields)
        
        %@12AUG2024: skip is main predictor is categorical
        if iscategorical(geno.(chrFields{i}).bed), continue; end

        missedCols = any(geno.(chrFields{i}).bed < 0, 1);
        NmissedCols = find(~missedCols);
        missedCols = find(missedCols);
        maf = nan(size(geno.(chrFields{i}).bed, 2), 1);
        for j = 1:numel(missedCols)
            mafTmp =  geno.(chrFields{i}).bed(:, missedCols(j));
            mafTmp(mafTmp < 0) = [];
            maf(missedCols(j), 1) = mean(mafTmp)./2;
        end
        maf(NmissedCols) = mean(geno.(chrFields{i}).bed(:, NmissedCols))./2;
        
        mafidx = maf < opts.maf; % filter maf below this cut-off
        if any(mafidx)
            if opts.verbose
                fprintf('%s: %d variants with MAF < %.3g were removed.\n',...
                    chrFields{i}, sum(mafidx), opts.maf)
            end
            rmfields = {'I_A', 'snp', 'chr', 'pos', 'a1', 'a2'};
            geno.(chrFields{i}).bed(:, mafidx) = [];
            for j = 1:numel(rmfields)
                if isfield(geno.(chrFields{i}), rmfields{j}) % bedread output doesn't have I_A field
                    geno.(chrFields{i}).(rmfields{j})(mafidx) = [];
                end
            end
        end
        
        % check if any variant has remained after MAF or INFO score filtering
        if isempty(geno.(chrFields{i}).snp)
            if opts.verbose
                fprintf('WARNING:no variant has remained for %s.\n', chrFields{i})
            end
            geno = rmfield(geno, chrFields{i});
        end
    end
    
    if isempty(fieldnames(geno)) % no variant present!
        fprintf('ERROR: no variant was found in geno struct. this may be due to filtering criteria.\n')
        return
    end
end

% prepare cvoariate matrix ------------------------------------------------
if isfield(opts, "covarFile") && opts.regenie
    opts.covar = ""; 
    
    %@26NOV2024: 'catCovarCell' overrides 'catCovar' with REGENIE when
    %phenoFile and covarFile are used.
    if all(isfield(opts, ["catCovarCell", "catCovar"]))
        opts = rmfield(opts, "catCovar");
    end

end

% check stratification phenotype
if opts.stratify ~= ""
    opts.covar = union(opts.covar, opts.stratify);
    opts.strata = true;
else
    opts.strata = false;
end

% categorical covariates are unified with 'covar' list
if isfield(opts, 'catCovar')
    opts.catCovarOut = string(opts.catCovar);
    opts = rmfield(opts, "catCovar");
    opts.catCovarOut(ismissing(opts.catCovarOut) | opts.catCovarOut == "") = [];
    if all(isempty(opts.catCovarOut)), opts.catCovarOut = []; end
else
    opts.catCovarOut = [];
end

if ~isempty(opts.catCovarOut)
    opts.covar = union(opts.covar, opts.catCovarOut, "stable");
end

if opts.desonly
    if opts.strata % if strata flag is used, strata pheno is needed
        opts.covar = opts.stratify;
    else
        opts.covar = "";
    end
end

if opts.interaction ~= "" && ~opts.desonly% load interactor along with covariates
    if ~any(ismember(opts.covar , opts.interaction))
        opts.covar = [opts.covar, opts.interaction];
    end
    opts.interactionIdx = find(ismember(opts.covar, opts.interaction)) + 1;
else
    opts.interactionIdx = nan;
end
    
tfile = fullfile(opts.phenopath, opts.covar);
if ~isempty(opts.catCovarOut)
    opts.catCovarOut = ismember(opts.covar, opts.catCovarOut);
end
get_covars = struct;

covar_tags = strings(numel(tfile), 1);
if all(opts.covar ~= "")
    % validate input covariates
    mustBePheno(string(opts.covar), opts.phenopath) % validate input covariates
    for i = 1:numel(tfile)
        t = load(tfile{i});
        covar_tags(i, 1) = string(t.UKB_STRUCT_ALL.tag);
        get_covars.(matlab.lang.makeValidName(t.UKB_STRUCT_ALL.tag)) = t.UKB_STRUCT_ALL; 
    end

    if ~isnan(opts.interactionIdx)
        checkFields = fieldnames(get_covars);
        opts.interactiontag = string(checkFields(opts.interactionIdx-1));
    end

    if ~isempty(opts.catCovarOut)
        checkFields = fieldnames(get_covars);
        opts.catCovarTagOut = string(checkFields(opts.catCovarOut));
    else
        opts.catCovarTagOut = [];
    end
    
    strataIdx = find(ismember(opts.covar, opts.stratify)); % for stratification
    opts.stratify = covar_tags(strataIdx);
    [covar_matrix_out, opts] = getcovar(get_covars, qc_eid_out, opts);
    opts.covartag = string(fieldnames(get_covars));
    asz = numel(opts.covartag)+1:size(covar_matrix_out, 2);
    opts.covartag = [opts.covartag; "COV" + asz'];

    %@13MAY2023: additional covariates
    if isfield(opts, "covarTab")
        ctab = rmmissing(opts.covarTab); 
        if ~isempty(ctab)
            [f1, f2] = ismember(qc_eid_out, ctab.eid); f2(f2<1) = [];
            cols = setdiff(colnames(ctab), union(opts.covartag, "eid"), "stable");
            covar_matrix_out2 =nan(numel(qc_eid_out), numel(cols));
            covar_matrix_out2(f1, :) = ctab{f2, cols};
            covar_matrix_out = [covar_matrix_out, covar_matrix_out2];
            clear covar_matrix_out2 ctab
            opts.covartag = [opts.covartag; cols'];
        end
        % opts = rmfield(opts, "covarTab");
        opts.covarTab = true;
    end

    if opts.verbose
        fprintf('input covariates: %s\n', strjoin(opts.covartag, ','))
        if opts.defaultcovar
            fprintf('\tfirst 10 genomic PCs are included in covariate matrix.\n')
            if ~opts.wes
                fprintf('\tgenotyping array batch is included in covariate matrix.\n')
            end
        end
        
        %@15NOV2024
        if ~isempty(opts.catCovarTagOut)
            fprintf("\tcategorical covariate:%s\n", opts.catCovarTagOut)
        end

        if ~isnan(opts.interactionIdx)
            fprintf("\tinteraction term:%s\n", opts.interactiontag)
        end

    end
    
    f_rem = any(isnan(covar_matrix_out), 2);
    if any(f_rem) && opts.verbose
        fprintf('\t%d missing samples in covariate matrix\n', sum(f_rem))
    end
else
    covar_matrix_out = []; % no covariates are queried 
end

%@17NOV2024: supports 'pred' file without regenie: adding as an extra
%covariate
if isfield(opts, "pred") && ~opts.regenie
    fprintf("\tfound pred files, will add as extra covariates (similar to REGENIE step 2)\n")
    opts.pred_covars = readPredFiles(opts.pred, fieldnames(geno));
end

% loop over traits --------------------------------------------------------
offset = 1; offset_des = [1, 1];
if opts.verbose; tic; end

% check traits.
if opts.trait == ""
    opts.trait = getPheno(true, false, opts.phenopath);
else
    mustBePheno(string(opts.trait), opts.phenopath) % validate input traits
    opts.trait = string(opts.trait);
end

% % @27MAY2024: 'interactionIdx' was modified to avoid nonsensical results
% % when interaction terms is the same as outcome
% opts.interactionIdxOut = opts.interactionIdx;
% opts = rmfield(opts, "interactionIdx");

for n = 1:numel(opts.trait)

    covar_matrix = covar_matrix_out;
    qc_eid = qc_eid_out;
    
    if isstruct(opts.trait)
        D = opts.trait(n);
        opts.trait(n).eid = "";
        opts.trait(n).rawUKB = "";
    else
        dv = fullfile(opts.phenopath, opts.trait(n));
        if ~endsWith(dv, '.mat'); dv = dv + ".mat"; end
        D = load(dv, "UKB_STRUCT_ALL");
        try
            D = D.UKB_STRUCT_ALL;
        catch
            if opts.verbose
                fprintf('%s is not a valid phenotype struct!\n', opts.trait(n))
            end
            continue
        end
    end

    if opts.verbose
        fprintf('%d (of %d)-trait under test: %s\n',n, numel(opts.trait), string(D.tag))
    end

    % @17MAY2024: survival analysis
    if isfield(D, "tt")
        opts.surv = true;
    else
        opts.surv = false;
    end
    
    trait = string(D.tag);
    % if trait under test is also a covariate, remove it from covariates
    if all(opts.covar ~= "")
        if isfield(opts, "covarTab")
            CovarCols = find(contains(opts.covartag, matlab.lang.makeValidName(D.tag))); % @02APR2024: ismember was replace with contains to be more flexible
        else
            CovarCols = find(ismember(opts.covartag, matlab.lang.makeValidName(D.tag)));
        end
        
        if ~isempty(CovarCols) && opts.verbose
            fprintf('\tWARNING: %s is also a covariate, removed %s from covariates!\n', trait, join(opts.covartag(CovarCols), ","))
        end

        CovarCols = setdiff(1:size(covar_matrix, 2), CovarCols);
        if ~isempty(opts.catCovarTagOut)
            opts.catCovarTag = intersect(opts.catCovarTagOut, opts.covartag(CovarCols), "stable");
            opts.catCovar = ismember(opts.covartag(CovarCols), opts.catCovarTag);
        else
            opts.catCovarTag = []; opts.catCovar = [];
        end
    else
        CovarCols = [];
    end
    
    % @27MAY2024: check if interaction term is the same as outcome
    % if ~isnan(opts.interactionIdxOut) && ~ismember(opts.interactionIdxOut-1, CovarCols)
    if ~isnan(opts.interactionIdx) && ~ismember(opts.interactionIdx-1, CovarCols) && ~isempty(covar_matrix)
        error('Error: interaction term cannot be the same as outcome!')
    end
    
    if ~isfield(D, 'exeid')
        D.exeid = [];
        f_exeid = false(numel(qc_eid), 1);
    else
        f_exeid = ismember(qc_eid, D.exeid);
        if opts.verbose
            fprintf('\t%d samples were removed due to exclusion criteria.\n',...
                sum(f_exeid))
        end
    end

    if opts.surv % check for those that are missing from trait (dropped out in the latest basket?: getQCEID should be updated)
        f_notInTrait = ~ismember(qc_eid, D.eid) & ~f_exeid;
        if any(f_notInTrait) && opts.verbose
            fprintf('\t%d samples were removed because they are not in the time-to-event data.\n',...
            sum(f_notInTrait))
        end
        f_exeid = f_exeid | f_notInTrait;

        % covariates to be used
        opts.surv_covaridx = CovarCols;
    end
    qc_eid(f_exeid) = [];    
    if ~isempty(covar_matrix); covar_matrix(f_exeid, :) = []; end
   

    if D.numericFlag  || strcmpi(D.tag, 'sex')
        % for sex 0 is female and 1 is male. Note that get_covar also codes
        % Sex in the same way. So in descriptive stat, n(%) shows men %
        [f1, f2] = ismember(qc_eid, D.eid); f2(f2 < 1) = [];
        DV_des = double(D.rawUKB(f2));
        if strcmpi(D.tag, 'sex')
            opts.skiptransform = true;

            % for descriptive stats
            scats = unique(D.rawUKB);
            for k = 1:numel(scats)
                opts.("sex_" + (k-1)) = unique(D.termMeaning(D.rawUKB == scats(k)));
            end
        else
            opts.skiptransform = false;
            if opts.verbose
                fprintf('\ttransformation function: %s\n', opts.transform)
            end
        end
        DV = DV_des;
        if ~opts.des && ~opts.barplot; clear DV_des; end
    else
        opts.skiptransform = true;
        DV = zeros(numel(qc_eid), 1);

        if opts.surv
            eventIdx = D.censor == 1;
            DV(ismember(qc_eid, D.eid(eventIdx))) = 1; % 1:event, 0:censored
        else
            DV(ismember(qc_eid, D.eid)) = 1;
        end
        f1 = true(numel(qc_eid), 1);
        if opts.des || opts.barplot; DV_des = DV; end
    end
    if ~isempty(covar_matrix); covar_matrix(~f1, :) = []; end
    qc_eid(~f1) = [];
    if opts.gpu
        DV = gpuArray(DV);
        if ~isempty(covar_matrix); covar_matrix = gpuArray(covar_matrix); end
    end
    
    opts.binary = ~D.numericFlag;
    if opts.regenie
        opts.regenied.trait = trait;
        opts.qc_eid = qc_eid;
        if isfield(opts, 'pred') % prediction file
            opts.regenied.pred = opts.pred(n);
        end

        %@29APR2024: support for 'phenoFile' and 'covarFile' was added
        if isfield(opts, 'phenoFile')
            opts.regenied.phenoFile = opts.phenoFile(n);
            opts.regenied.covarFile = opts.covarFile(n);

            %@26NOV2024: 'catCovarCell': categorical covariates per each
            %pheno
            if isfield(opts, "catCovarCell")
                opts.regenied.catCovar = opts.catCovarCell(n);
            end
        end

    end
    
    % for WES redirect to rare variant analysis and skip the rest of loop
    if opts.wesin && ~opts.surv
        if opts.regenie
            opts.regenied.geno = geno;
            if numel(snps.gene) > 1 && ~geno.regenie_merger %@26JULY2024 to allow for region based and merging genes into a single region
                opts.regenied.geno.multi = true; % multiple genes
            else
                opts.regenied.geno.multi = false; % merger of multiple genes or single gene
            end
            [skat{n}, opts] = call_regenie_in(covar_matrix(:, CovarCols), DV, opts);
            opts.regenied.snplistAll{n, 1} = opts.regenied.snplist;
        else
            opts.westrait = trait;
            [skat{n}, collapseres{n, 1}] = callSKAT(geno, DV, covar_matrix(:, CovarCols), qc_eid, opts);
            
        end
        if opts.verbose; fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'); end
        continue
    end
    
    % loop over rows of res
    if ~opts.regenie
        chr = fieldnames(geno);
    else
        chr = cellstr(string(unique(opts.regenied.achr)));
    end
    [tabALL, tabDES] = deal(struct("Pheno", [], "SNP", [], "CHR", [], "POS", [], "A1", [], "A2", []));
    tabALL.ct = 1; tabALL.k = 1;
    opts.firthTime = 0; % time spent in logistf R package
            
%     if opts.regenie 
%         opts.regenied.trait = trait;
%         call_regenie_in(covar_matrix(:, CovarCols), DV, qc_eid, opts);
%         return
%     end

    for j = 1:numel(chr) % loop over chromosomes in raw geno ==============
        if opts.regenie
            covarIdx = (1:numel(qc_eid)).';
            idx = ismember(opts.regenied.achr, chr{j});
            opts.regenied.chr = opts.regenied.achr(idx);
            [chunk.snp, opts.regenied.snp] = deal(opts.regenied.asnp(idx));
            chunk.bed = false(numel(qc_eid), numel(chunk.snp)); % fake snp matrix
            
            idx = ismember(opts.regenied.abgenchr, chr{j});
            opts.regenied.bgenfile = opts.regenied.abgenfile(idx);

        else
            chunk = geno.(chr{j});
            [fkeep, covarIdx] = ismember(chunk.eid, qc_eid);
            chunk.bed(~fkeep, :) = [];
        end
        % chunk.bed = gpuArray(chunk.bed); % unnecessary memory overhead in parallel mode
        covarIdx(covarIdx < 1) = [];
        
        if opts.surv
            opts.fitter = "glm"; 
            [~, idx2] = ismember(qc_eid(covarIdx), D.eid); idx2(idx2 < 1) = [];
    
            % keep subset of trait df in opts
            opts.surv_data = table(D.eid(idx2), D.censor(idx2), D.tt(idx2), ...
                D.tt0(idx2), D.base(idx2), ...
                VariableNames=["eid", "censor", "tt", "tt0", "base"]);
            opts.qc_eid = qc_eid(covarIdx);
            if isfield(D, "ex")
                opts.surv_ex = D.ex;
            else
                opts.surv_ex = table;
            end
        else
            opts.qc_eid = qc_eid(covarIdx);
        end

        if ~opts.regenie
            tabALL.Pheno = [tabALL.Pheno; repmat(trait, numel(chunk.snp), 1)];
            tabALL.SNP = [tabALL.SNP; chunk.snp];
            tabALL.CHR = [tabALL.CHR; chunk.chr];
            tabALL.POS = [tabALL.POS; chunk.pos];
            tabALL.A1 = [tabALL.A1; chunk.a1];
            tabALL.A2 = [tabALL.A2; chunk.a2];
            
            if opts.des || opts.barplot
                tabDES.Pheno = [tabDES.Pheno; repmat(trait, numel(chunk.snp), 1)];
                tabDES.SNP = [tabDES.SNP; chunk.snp];
                tabDES.CHR = [tabDES.CHR; chunk.chr];
                tabDES.POS = [tabDES.POS; chunk.pos];
                tabDES.A1 = [tabDES.A1; chunk.a1];
                tabDES.A2 = [tabDES.A2; chunk.a2];
            end
        end
        
        if opts.verbose; fprintf('\t%d (of %d)-%s', j, numel(chr), chr{j}); end
        
        if opts.ld && ~opts.regenie % computes Pearson r2
            bed = gather(chunk.bed);
            bed(any(bed < 0, 2), :) = []; % remove missing calls
            opts.r2{j} = tril(corr(bed, 'type', 'Pearson').^2);
            opts.r2{j}(opts.r2{j} == 0) = nan;
            opts.r2{j} = round(opts.r2{j}, 4, 'significant');
            opts.r2{j} = array2table(tril(opts.r2{j}), 'VariableNames', chunk.snp);
            opts.r2{j}.Row = chunk.snp;
            opts.r2{j}.Properties.Description = string(chr{j}) + ".r2";

            clear bed
        end
        
        if opts.fitter == "scoretest"
            % only if there is no missingness in genotype/dosage data, otherwise it
            % should be calculated separately for each variant (currently not supported)
            % scoretestCond = ~opts.firth && ~any(chunk.bed < 0, 'all') && ~isempty(covar_matrix) && isnan(opts.interactionIdx);
            scoretestCond = ~any(chunk.bed < 0, 'all') && ~isempty(covar_matrix) && isnan(opts.interactionIdx);
            if ~scoretestCond % doesn't support firth/interaction/empty covariates/missingness in genotypes
                opts.fitter = "glm"; % reset model fitter 
            else
                opts.chunkIdx = unique([1:opts.chunk:numel(chunk.snp), numel(chunk.snp)+1]);
            end
        end
         
        if opts.strata
            %% stratification =============================================
            if isempty(covar_matrix)
                error('gwasrunner: stratification cannot be performed without covariates!')
            end
            covar_matrix_tmp = covar_matrix(covarIdx, :);
            CovarCols_tmp = CovarCols;
            
            % get strata values
            stratPheno = gather(covar_matrix_tmp(:, strataIdx));
            checkBinary = unique(stratPheno(~isnan(stratPheno)));
            if numel(checkBinary) == 2 && all(ismember(checkBinary, [0, 1])) % binary pheno
                strata = 1; % < 1 (i.e. 0) control, and >= 1 (i.e. 1) case
                CovarCols_tmp(CovarCols_tmp == strataIdx) = [];
                checkBinaryPheno = true;
                % fprintf('%s is a binary variable, removed from covariates to avoid rank deficiency.\n', opts.stratify)
            else
                if opts.removeStrat
                    CovarCols_tmp(CovarCols_tmp == strataIdx) = [];
                end
                checkBinaryPheno = false;
                if all(~isnan(opts.quantile)) % first check quantiles
                    strata = quantile(stratPheno(~isnan(stratPheno)), opts.quantile);
                elseif opts.stratifyby == "bins" % then check custom strata
                    strata = opts.bins;
                else % mean or median
                    strata = feval(opts.stratifyby, stratPheno(~isnan(stratPheno)));
                end
            end
            
            if opts.verbose
                if ~opts.parallel && (ismember(opts.fitter, ["glm", "regenie"]) || opts.desonly)% for parallel use different progress bars for each stratum
                    progressBar = progressGen(numel(chunk.snp)*(numel(strata) + 1));
                    opts.progressBar = progressBar;
                elseif opts.parallel && opts.desonly
                    progressBar = progressGen(numel(chunk.snp)*(numel(strata) + 1));
                elseif opts.fitter == "scoretest"
                    progressBar = progressGen((numel(opts.chunkIdx)-1)*(numel(strata) + 1));
                    opts.progressBar = progressBar;
                end
            end
            
            if opts.des || opts.barplot; DV_des_tmp = DV_des(covarIdx); end
            DV_tmp = DV(covarIdx);
            stratTag = string;
            for k = 1:numel(strata) + 1
                
                stableCt = tabALL.ct; % reset ct after each stratum run
                if k == 1
                    sidx = gather(stratPheno < strata(k));
                    if checkBinaryPheno
                        if strcmpi(opts.stratify, "sex")
                            stratTag(1, k) = opts.stratify + " " + opts.sex_0.lower;
                        else
                            stratTag(1, k) = opts.stratify + " (0)";
                        end
                    else
                        stratTag(1, k) = opts.stratify + "<" + round(strata(k), 2);
                    end
                elseif k > numel(strata)
                    if checkBinaryPheno
                        if strcmpi(opts.stratify, "sex")
                            stratTag(1, k) = opts.stratify + " " + opts.sex_1.lower;
                        else
                            stratTag(1, k) = opts.stratify + " (1)";
                        end
                    else
                        stratTag(1, k) = opts.stratify + "≥" + round(strata(k-1), 2);
                    end
                    sidx = gather(stratPheno >= strata(k-1));
                else
                    stratTag(1, k) = round(strata(k-1), 2) + "≤" + opts.stratify + "<" + round(strata(k), 2);
                    sidx = gather(stratPheno < strata(k) & stratPheno >= strata(k-1));
                end

                opts.stratTag = stratTag;
                
                if opts.verbose
                    if opts.parallel && opts.fitter == "glm" && ~opts.desonly % for parallel mode
                        fprintf('\n\tstratum: %d of %d ', k, numel(strata) + 1);
                        progressGen(numel(chunk.snp)); % only generate an empty progress bar
                    end
                end
                            
                tabALL.k = k;
                if ismember(opts.fitter, ["glm", "regenie"])
                    opts.strataCounter = size(chunk.bed, 2)*(k-1)+1:size(chunk.bed, 2)*k;
                elseif opts.fitter == "scoretest"
                    opts.strataCounter = (numel(opts.chunkIdx)-1)*(k-1)+1:(numel(opts.chunkIdx)-1)*k; % a different counter is needed in this case
                end

                if opts.regenie || opts.surv, opts.qc_eid = qc_eid(sidx); end
                    
                if ~opts.desonly
                    if ~strcmpi(D.tag, 'sex') && ~opts.binary && ~strcmp(opts.transform, 'none')
                        % transformation must be done per each stratum
                        DV_tmp_stratum = DV_tmp(sidx);
                        if opts.gpu; DV_tmp_stratum = gpuArray(DV_tmp_stratum); end

                        [tabALL, opts] = runModels(chunk.bed(sidx, :), ...
                            covar_matrix_tmp(sidx, CovarCols_tmp), DV_tmp_stratum,...
                            opts, tabALL);
                    else
                        
                        [tabALL, opts] = runModels(chunk.bed(sidx, :), ...
                            covar_matrix_tmp(sidx, CovarCols_tmp), DV_tmp(sidx),...
                            opts, tabALL);
                    end
                    
                end
                                 
                if opts.des || opts.barplot % descriptive statisticis 
                    len2 = size(tabDES.SNP, 1);
                    len1 = len2 + 1 - size(chunk.bed, 2);
                    len12 = len1:len2;
                    for k2 = 1:numel(len12)
                        tabDES = makedes(tabDES, chunk.bed(sidx, k2), DV_des_tmp(sidx), len12(k2), k, opts);
                        tabDES.Strata(len12(k2), k) = stratTag(1, k);
                        if opts.desonly && opts.verbose
                            progressGen(progressBar, opts.strataCounter(k2))
                        end
                    end
                end
                
                if k <= numel(strata) % don't reset counter for last stratum
                    tabALL.ct = stableCt;
                end
            end
            
            for k = stableCt:(tabALL.ct - 1) % add strata labels
                tabALL.Strata(k, :) = stratTag;
            end
            
            opts.strataCounter = 1;
            tabALL.k = 1; % reset (for future implementation of both stratification and overall analysis)
        
        
        else 
            %% overall ====================================================
            % initialize the progress bar
            if opts.verbose
                if ismember(opts.fitter, ["glm", "regenie"]) || opts.desonly
                    progressBar = progressGen(numel(chunk.snp));
                elseif opts.fitter == "scoretest"  % for "scoretest", progress bar length is different
                    progressBar = progressGen(numel(opts.chunkIdx)-1);
                end
                opts.progressBar = progressBar;
            end
            
            if ~opts.desonly
                if ~isempty(covar_matrix)
                    [tabALL, opts] = runModels(chunk.bed,...
                        covar_matrix(covarIdx, CovarCols),...
                        DV(covarIdx), opts, tabALL);
                else
                    [tabALL, opts] = runModels(chunk.bed, [],...
                        DV(covarIdx), opts, tabALL);
                end
            end
            
            if isfield(opts, 'strataCounter'); opts = rmfield(opts, 'strataCounter'); end
            
            if opts.des
                len2 = size(tabDES.SNP, 1);
                len1 = len2 + 1 - size(chunk.bed, 2);
                len12 = len1:len2;
                for k2 = 1:numel(len12)
                    tabDES = makedes(tabDES, chunk.bed(:, k2), DV_des(covarIdx), len12(k2), 1, opts);
                    if opts.desonly && opts.verbose
                        progressGen(progressBar, k2)
                    end
                end
            end    
        end
        
        if opts.verbose; fprintf('\n'); end
        clear chunk
        continue
    end
    
    if isfield(tabALL, 'ct')
        tabALL = rmfield(tabALL, {'ct', 'k'});
    end

    if opts.firth && opts.binary && ~opts.desonly && opts.verbose
        fprintf('\ttime spent with Firth logistic regression: %.2f s\n', opts.firthTime)
    end
    if opts.verbose; fprintf('\n'); end
    clear covar_matrix DV
    
    if opts.barplot % bar plot for stratified trait and genotypes
        opts.barfunc = string(opts.barfunc);
        if numel(opts.barfunc) == numel(opts.trait) % different barfunc for each trait
            cn = n;
        else
            cn = 1;
        end
        barPlotter(tabDES, 'barfunc', opts.barfunc(cn), ...
            'stratify', opts.stratify, 'output', opts.output, ...
            'saveplot', opts.saveplot, 'savefig', opts.savefig, ...
            'resolution', opts.resolution, 'format', opts.format)
    end
    
    if ~opts.surv %#SHOULD BE IMPLEMENTED
        if opts.plot && ~opts.desonly && isnan(opts.interactionIdx)
            if opts.singlerare % for rare variants use rarePlotter
                raretab = struct2table(tabALL);
                [idx1, idx2] = ismember(raretab.SNP, opts.wesinfo.id);
                raretab.consequence(idx1) = opts.wesinfo.consequence(idx2);
                rarePlotter(raretab, 'format', opts.format, ...
                    'resolution', opts.resolution,...
                    'output', opts.output + matlab.lang.makeValidName(raretab.Pheno(1)), ...
                    'save', opts.saveplot, 'savefig', opts.savefig, ...
                    'raremaf', opts.raremaf);
            else
                try
                    effPlotter(tabALL, opts)
                catch % empty summary stats
                    % should be taken care of 
                end
            end
        end
    end

     % write data to a file ===================================================
    if opts.strata && ~opts.desonly && ~opts.surv
        fnames = fieldnames(tabALL);
        for i = 1:size(tabALL.Strata, 2)
            tmp = struct;
            for j = 1:numel(fnames)
                if iscolumn(tabALL.(fnames{j}))
                    tmp.(fnames{j}) = tabALL.(fnames{j});
                else
                    tmp.(fnames{j}) = tabALL.(fnames{j})(:, i);
                end
            end
            tmp = toCleanTab(tmp, opts);
            tmp.Strata(2:end) = missing;
            if opts.sheet == "Overall"
                sheetName = tmp.Strata(1);
            else
                sheetName = opts.sheet + "." + tmp.Strata(1);
            end
            if opts.firth
                sheetName = sheetName + "(Firth)";
            end
            if n == 1
                writemode = 'overwritesheet';
            else
                writemode = 'inplace';
            end
            
            writetable(tmp, opts.output + ".xlsx", 'AutoFitWidth', true, ...
                    'WriteMode', writemode, 'WriteVariableNames', true,...
                    'Sheet', replace(trimSheetName(sheetName), ["?", "/", "\"], "_"), 'Range', "A"+offset, 'AutoFitWidth', true)
            opts.cleanOutputcheck = unique([opts.cleanOutputcheck; sheetName]); 
        end
        offset = offset + size(tmp, 1) + 2;
    elseif ~opts.desonly && ~opts.surv
        if n == 1
            writemode = 'overwritesheet';
        else
            writemode = 'inplace';
        end

        try
            tabALL = toCleanTab(tabALL, opts);
            tabALL = fillmissing(tabALL, "previous", "DataVariables", "Pheno"); %@27MAY2024
            writetable(tabALL, opts.output + ".xlsx", 'AutoFitWidth', true, ...
                        'WriteMode', writemode, 'WriteVariableNames', true,...
                        'Sheet', opts.sheet, 'Range', "A"+offset, 'AutoFitWidth', true)
            opts.cleanOutputcheck = unique([opts.cleanOutputcheck; opts.sheet]); 
            offset = offset + size(tabALL, 1) + 2;
        catch
            % with REGENIE can happen
            continue
        end
    end
    
    if opts.des % for descriptive statistics ------------------------------
        if opts.strata
            fnames = fieldnames(tabDES);
            for i = 1:size(tabDES.Strata, 2)
                tmp = struct; offset_des_tmp = offset_des(2);
                for j = 1:numel(fnames)
                    if iscolumn(tabDES.(fnames{j}))
                        tmp.(fnames{j}) = tabDES.(fnames{j});
                    elseif numel(size(tabDES.(fnames{j}))) > 2
                        tmp.(fnames{j}) = tabDES.(fnames{j})(:, i, :);
                    else
                        tmp.(fnames{j}) = tabDES.(fnames{j})(:, i);
                    end
                end
                [tmp, tmp2] = toCleanTabDES(tmp, opts);
                if height(tabDES) > 1; tmp.Pheno(2:end) = missing; end
                tmp.Strata(2:end) = missing;
                tmp = movevars(tmp, 'Strata', 'After', 'median(IQR)');
                sheetName = tmp.Strata(1) + ".desc";
                if n == 1
                    writemode = 'overwritesheet';
                else
                    writemode = 'inplace';
                end
                writetable(tmp, opts.output + ".xlsx", 'AutoFitWidth', true, ...
                        'WriteMode', writemode, 'WriteVariableNames', true,...
                        'Sheet', replace(trimSheetName(sheetName), ["?", "/", "\"], "_"), 'Range', "A"+offset_des(1), 'AutoFitWidth', true)
                opts.cleanOutputcheck = unique([opts.cleanOutputcheck;sheetName]);    
                sheetName = tmp.Strata(1) + ".geno.desc";    
                opts.cleanOutputcheck = unique([opts.cleanOutputcheck;sheetName]);
                for k = 1:numel(tmp2)
                    if n == 1 && k == 1
                        writemode = 'overwritesheet';
                    else
                        writemode = 'inplace';
                    end
                    writetable(tmp2{k}, opts.output + ".xlsx", 'AutoFitWidth', true, ...
                        'WriteMode', writemode, 'WriteVariableNames', true,...
                        'Sheet', replace(trimSheetName(sheetName), ["?", "/", "\"], "_"), 'Range', "A"+offset_des_tmp, 'AutoFitWidth', true)
                    offset_des_tmp = offset_des_tmp + 3;
                end
            end
            offset_des(2) = offset_des_tmp;
            offset_des(1) = offset_des(1) + size(tmp, 1) + 2;
        else % ------------------------------------------------------------
            % @ consider editing
            if opts.verbose
                fprintf('!!!! consider editing toCleanTabDES file:unify headers !!!!\n')
            end
            
            [tabDES, tabDES2] = toCleanTabDES(tabDES, opts);
            if height(tabDES) > 1; tabDES.Pheno(2:end) = missing; end
            if n == 1
                writemode = 'overwritesheet';
            else
                writemode = 'inplace';
            end
            writetable(tabDES, opts.output + ".xlsx", 'AutoFitWidth', true, ...
                    'WriteMode', writemode, 'WriteVariableNames', true,...
                    'Sheet', replace(trimSheetName(opts.sheet), ["?", "/", "\"], "_")+".desc", 'Range', "A"+offset_des(1), 'AutoFitWidth', true)
            opts.cleanOutputcheck = unique([opts.cleanOutputcheck; opts.sheet+".desc"]);
            offset_des(1) = offset_des(1) + size(tabDES, 1) + 2;
            for k = 1:numel(tabDES2)
                if n == 1 && k == 1
                    writemode = 'overwritesheet';
                else
                    writemode = 'inplace';
                end
                writetable(tabDES2{k}, opts.output + ".xlsx", 'AutoFitWidth', true, ...
                    'WriteMode', writemode, 'WriteVariableNames', true,...
                    'Sheet', replace(trimSheetName(opts.sheet), ["?", "/", "\"], "_")+".geno.desc", 'Range', "A"+offset_des(2), 'AutoFitWidth', true)
                offset_des(2) = offset_des(2) + 3;
            end
            opts.cleanOutputcheck = unique([opts.cleanOutputcheck; opts.sheet+".geno.desc"]);
        end
    end
        
    if opts.verbose; fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'); end
end % trait loop ----------------------------------------------------------

if opts.ld
    for i = 1:numel(opts.r2)
        opts.r2{i}.Properties.DimensionNames = {' ','Variables'};
        writetable(opts.r2{i}, opts.output + ".xlsx", 'WriteMode',...
            'overwritesheet', 'WriteRowNames', true, ...
            'Sheet', opts.r2{i}.Properties.Description)
    end
end

if opts.wesin
    writeWESresults(skat, collapseres, opts) % write SKAT results
elseif opts.cleanOutput % further prune tables already stored
    pruneMore(opts)
end

if (opts.collapse || opts.singlerare) && isfield(opts, 'wesinfo')  % write wesinfo
    writetable(opts.wesinfo, opts.output+".xlsx", 'AutoFitWidth', true, ...
    'Sheet', 'Info', 'WriteMode', 'overwritesheet')
end

% write survival results
if isfield(opts, "surv_results")
    writetable(opts.surv_results, opts.output + ".xlsx", ...
                    'WriteVariableNames', true,...
                    'Sheet', opts.sheet + "_surv", 'AutoFitWidth', true)
end

if opts.verbose; toc; end
if opts.verbose; fprintf('finished at %s\n\n', string(datetime('now'))); end
if opts.log 
    diary off
end

if ~opts.warning
    warning on
end
end % END


%% subfunctions =========================================================== 
function [tabALL, opts] = runModels(bed, covariates, DV, opts, tabALL)
% fits a glm to data and returns summary stat for snp in two models
% unadjusted and adjusted for covariates.

% create empty struct for results
tabEmpty = struct('A2FREQ', nan, 'p', nan, 'p_adj', nan, ...
    'BETA', nan, 'SE', nan, 'CI', "", 'N', nan, 'N_case', nan, ...
    'AF_case', nan, 'AF_control', nan);
if opts.firth && opts.binary
    tabEmpty.firthTime = 0;
end
if ~isnan(opts.interactionIdx)
    [tabEmpty.p_int, tabEmpty.p_int_adj] = deal(nan);
end
if opts.firth && opts.binary
    if opts.flic
        [tabEmpty.p_flic, tabEmpty.p_adj_flic, tabEmpty.BETA_flic] = deal(nan);
    end
    if opts.flac
        [tabEmpty.p_flac, tabEmpty.p_adj_flac, tabEmpty.BETA_flac] = deal(nan);
    end
end

% initialize progress bar for chunk/strata
tabC = ({}); firthTime = (0);
if ~isfield(opts, 'strataCounter') % for overall analysis
    if opts.fitter == "scoretest"
        opts.strataCounter = 1:(numel(opts.chunkIdx)-1);
    else
        opts.strataCounter = 1:size(bed, 2);
    end
% else
%     opts.strataCounter = opts.strataCounter(1):opts.strataCounter(2);
end

% use REGENIE for fitting 
if opts.regenie
    [tabC, opts] = call_regenie_in(covariates, DV, opts);

    if isempty(tabC)
        tabC = ({});
        tabC(1:size(bed, 2)) = {tabEmpty};
        for i = 1:size(bed, 2)
            tabALL = appendTab(tabALL, tabC{i}, tabALL.ct, tabALL.k, opts);
            tabALL.ct = tabALL.ct + 1;
        end
    else
        cols = string(tabC.Properties.VariableNames);
        ctrange = tabALL.ct:(tabALL.ct + height(tabC)) - 1;
        
        % remove empty fields (to avoid strings to be saved as NaNs)
        rmidx = structfun(@isempty, tabALL);
        tabFis = fieldnames(tabALL);
        tabALL = rmfield(tabALL, tabFis(rmidx));

        for i = 1:numel(cols)
            tabALL.(cols(i))(ctrange, tabALL.k) = tabC.(cols(i));
        end
        tabALL.ct = ctrange(end) + 1;
    end
    
    if opts.verbose
        for i = 1:size(bed, 2)
            progressGen(opts.progressBar, opts.strataCounter(i))
        end
    end

    return
end

if strcmp(opts.model, 'add')
    opts.denom = 2; % allele freq.
else
    opts.denom = 1; % geno freq. for rec/dom
end


interactStruct = struct;
if isfield(opts, "catCovar")
    opts.catCovar2 = find(opts.catCovar) + 1; % 1st predictor is genetic variant
    opts.catCovar1 = [];
    if ~isnan(opts.interactionIdx)
        covarlbls = join("x" + setdiff((2:size(covariates, 2)+1), opts.interactionIdx), '+');
        spec1 = "y ~ " + join(["x1", "x" + opts.interactionIdx], '*');
        spec2 = spec1 + "+" + covarlbls;
        interactStruct.interactTerm2 = join(["x1", "x" + opts.interactionIdx], ':');
        interactStruct.interactTerm1 = "x1:x2"; spec1 = "y ~ x1*x2";
        if any(ismember(opts.interactionIdx, opts.catCovar2))
            opts.catCovar1 = 2; % for spec1
            interactStruct.interactTerm1 = interactStruct.interactTerm1 + "_";
            interactStruct.interactTerm2 = interactStruct.interactTerm2 + "_";
        end
    else
        spec1 = "y~x1";
        spec2 = "y~" + join("x"+(1:size(covariates, 2)+1), '+');
    end
else
    opts.catCovar1 = []; opts.catCovar2 = [];
    spec1 = "y~x1";
    spec2 = "y~" + join("x"+(1:size(covariates, 2)+1), '+');
end

%@17NOV2024: check if 'pred' files have been provided
if isfield(opts, "pred") && ~opts.regenie

    % find pred matrix index
    tmp_pheno = matlab.lang.makeValidName(tabALL.Pheno(tabALL.ct));
    if opts.surv, tmp_pheno = tmp_pheno + "_TIME"; end
    [~, idx] = ismember(tmp_pheno, string(opts.pred_covars(:, 1)));
    idx(idx < 1) = [];
    if isempty(idx)
        fprintf("\twarning: cannot find predfile for %s\n", tabALL.Pheno(tabALL.ct))
    else
        
        pred_data = opts.pred_covars{idx, 2};
        pred_data = pred_data(:, ["eid", "chr" + tabALL.CHR(tabALL.ct)]);
        [f1, f2] = ismember(opts.qc_eid, pred_data.eid); f2(f2 < 1) = [];
        
        pred_covar = nan(numel(f1), 1);
        pred_covar(f1) = pred_data{f2, 2};
        clear pred_data
        if opts.gpu
            pred_covar = gpuArray(pred_covar);
        end
        
        % add pred vector to covariates
        if isempty(covariates)
            covariates = pred_covar;
        else
            covariates = [covariates, pred_covar];
        end
        
        % update spec2
        max_cov_num = max(double(spec2.erase(" ").split(["*", "+", ":", "~"] + "x")));
        spec2 = spec2 + "+x" + max_cov_num;

        % update covartag for survival analysis
        if opts.surv
            opts.covartag(opts.covartag == "COV_REGENIE_PRED") = [];
            n_index = numel(opts.covartag) + 1;
            opts.covartag = [opts.covartag; "COV_REGENIE_PRED"];
            if ~any(opts.surv_covaridx == n_index)
                opts.surv_covaridx = union(opts.surv_covaridx, n_index, "stable");
            end
        end
    end
end

if isempty(covariates)
    fnan = false(numel(DV), 1);
else
    fnan = any(isnan(covariates), 2); % nan values in covariates
end

% DVclean -> for adjusted model
% DV -> for unadjusted univariate model
DVclean = DV(~fnan);

if numel(DVclean) < 2 % numel DV < 2
    tabC(1:size(bed, 2)) = {tabEmpty};
    for i = 1:size(bed, 2)
        tabALL = appendTab(tabALL, tabC{i}, tabALL.ct, tabALL.k, opts);
        tabALL.ct = tabALL.ct + 1;
        if opts.verbose, progressGen(opts.progressBar, opts.strataCounter(i)); end
    end
    return % skip model fit
end

if ~strcmp(opts.transform, "none") && ~opts.skiptransform
    DVclean = feval(opts.transform, gather(DVclean));
    DV = feval(opts.transform, gather(DV));

    if opts.gpu
        DVclean = gpuArray(DVclean);
        DV = gpuArray(DV);
    end
end

% mean impute missing genotypes (similar to REGENIE -> github issue #70)
% @20MAY2023: negative values are no longer considered as missingness.
% bed(bed < 0) = nan; 

if ~iscategorical(bed)
    fmiss = isnan(bed);
    if any(fmiss, 'all')
        fprintf('missing genotypes will be mean-imputed.\n')
        cmiss = find(any(fmiss, 1));
        bedmean = mean(bed(:, cmiss), 1, 'omitnan');
        for i = 1:numel(cmiss)
            bed(fmiss(:, cmiss(i)), cmiss(i)) = bedmean(i);
        end
    end
end

% score test 
if opts.fitter == "scoretest"
    [tabALL, opts] = scoretest(DVclean, covariates, bed, opts, tabALL);
    return
end

% survival analysis
if opts.surv
    keepRange = tabALL.ct:numel(tabALL.Pheno);    
    opts.surv_preds = tabALL;
    rmfi = setdiff(fieldnames(opts.surv_preds), ["Pheno", "SNP", "CHR", "POS", "A1", "A2"]);
    opts.surv_preds = rmfield(opts.surv_preds, rmfi);
    fis = string(fieldnames(opts.surv_preds));
    for k = 1:numel(fis)
        opts.surv_preds.(fis(k)) = opts.surv_preds.(fis(k))(keepRange);
    end
    survFitResults = fitSurvCox(bed, DV, covariates, opts, plotFlag=opts.plot, ...
        checkDiagnostics=opts.checkDiagnostics, logRank=opts.logRank, ...
        FGR=opts.FGR, prevalent=opts.prevalent, timeOrigin=opts.timeOrigin, ...
        LT=opts.LT, logRankFun=opts.logRankFun, landmark=opts.landmark);
    if isfield(opts, "surv_results")
        opts.surv_results = [opts.surv_results; survFitResults];
    else
        opts.surv_results = survFitResults;
    end

    tabC(1:size(bed, 2)) = {tabEmpty};
    for i = 1:size(bed, 2)
        tabALL = appendTab(tabALL, tabC{i}, tabALL.ct, tabALL.k, opts);
        tabALL.ct = tabALL.ct + 1;
        if opts.verbose, progressGen(opts.progressBar, opts.strataCounter(i)); end
    end

    if opts.verbose
        for i = 1:size(bed, 2)
            progressGen(opts.progressBar, opts.strataCounter(i))
        end
    end
    return
                
end

if opts.parallel % parfor seems the best option
    parInfo.tabEmpty = tabEmpty;
    parInfo.denom = opts.denom;
    parInfo.spec1 = spec1; parInfo.spec2 = spec2;
    parInfo.interactStruct = interactStruct;
    parInfo.fnan = fnan;

    [tabC, opts] = runModelsPar(bed, covariates, DV, DVclean, opts, parInfo);

else % single-threaded ----------------------------------------------------
    for i = 1:size(bed, 2) % loop over snps
        snp = bed(:, i);

        if ~any(snp) || numel(DVclean) < 2 % either empty variant or numel < 2
            tabC{i, 1} = tabEmpty;
            continue
        end

        if opts.gpu; snp = gpuArray(snp); end

        if opts.binary
            if opts.firth % perform Firth logistic regression (to deal with imbalanced/separation issues)
                Fopts = opts;
                Fopts.snpIdx = 1; %Fopts.n = size(covariates, 2)-1;
                [mdl1, mdl2, firthTime(i)] = logistf(snp, covariates, DV, Fopts, fnan); % R logistf removes nan values when necessary
                if ~isnan(opts.interactionIdx)
                    interactStruct = struct;
                    interactStruct.interactTerm1 = mdl1.interactTerm;
                    interactStruct.interactTerm2 = mdl2.interactTerm;
                end

            else % binomial GLM or Cox regression

                if ~isnan(opts.interactionIdx)
                    mdl1 = fitglm([snp, covariates(:, opts.interactionIdx-1)], DV, spec1, 'Distribution', 'binomial', 'Link', 'logit', 'CategoricalVars', opts.catCovar1);
                else
                    mdl1 = fitglm(snp, DV, spec1, 'Distribution', 'binomial', 'Link', 'logit', 'CategoricalVars', opts.catCovar1);
                end

                if isempty(covariates)
                    mdl2 = mdl1;
                else
                    mdl2 = fitglm([snp(~fnan), covariates(~fnan, :)], DVclean, spec2, 'Distribution', 'binomial', 'Link', 'logit', 'CategoricalVars', opts.catCovar2); % , 'LikelihoodPenalty', "jeffreys-prior");
                end
                
            end
        else % continuous traits --------------------------------------------------
            if ~isnan(opts.interactionIdx)
                mdl1 = fitlm([snp, covariates(:, opts.interactionIdx-1)], DV, spec1, 'CategoricalVars', opts.catCovar1);
            else
                mdl1 = fitlm(snp, DV, spec1, 'CategoricalVars', opts.catCovar1);
            end

            if isempty(covariates)
                mdl2 = mdl1;
            else
                mdl2 = fitlm([snp(~fnan), covariates(~fnan, :)], DVclean, spec2, 'CategoricalVars', opts.catCovar2);
            end
        end

        if opts.debug
            format("compact")
            disp(mdl1)
            disp(mdl2)
            dfig = figure;
            dtile = tiledlayout(dfig, "flow", "TileSpacing", "compact");
            dax = nexttile(dtile);
            plotResiduals(mdl2, 'probability', "Parent", dax)
            dax = nexttile(dtile);
            plotResiduals(mdl2, 'histogram', "Parent", dax)
            dax = nexttile(dtile);
            plotResiduals(mdl2, 'fitted', "Parent", dax)
            dtile.Title.String = tabALL.Pheno;
            if isfield(opts, "stratTag")
                dtile.Title.String = dtile.Title.String + ":" + opts.stratTag(end);
            end
            try disp(vif(mdl2)); catch; end
            format("default")
        end

        tab = struct; % faster than table
        tab.A2FREQ = gather(mean(snp)./opts.denom); % freq calculated here is independent of missing values in covariates

        snp(fnan) = []; % effective sample in adjusted model
        tab.N = numel(snp); % N corresponds to effectize sample size used for adjusted model
        if opts.binary
            tab.N_case = gather(sum(DVclean)); % effective case N
            tab.AF_case = gather(mean(snp(DVclean == 1))./opts.denom);
            tab.AF_control = gather(mean(snp(DVclean == 0))./opts.denom);
        end

        if ~isnan(opts.interactionIdx)
            tabC{i, 1} = extractModelInfo(compact(mdl1), compact(mdl2), opts, tab, interactStruct);
        else
            try
                tabC{i, 1} = extractModelInfo(compact(mdl1), compact(mdl2), opts, tab);
            catch % for Firth's model
                tabC{i, 1} = extractModelInfo(mdl1, mdl2, opts, tab);
            end
        end
        if opts.verbose, progressGen(opts.progressBar, opts.strataCounter(i)); end % update progress bar 
        
        % draw interaction plots (doesn't support Firth regression now)
        if opts.plot && ~isnan(opts.interactionIdx) && ~opts.firth
            incnt = tabALL.ct+i-1;
            labels = [tabALL.SNP(incnt), opts.interactiontag];

            output = matlab.lang.makeValidName(labels);
            output = join(output, "X") + "_" + matlab.lang.makeValidName(tabALL.Pheno(incnt));
            output = fullfile(fileparts(opts.output), output);
            if isfield(opts, "stratTag")
                stratTag = replace(opts.stratTag, ...
                    [">", "≥", "<", "≤", "=="], ...
                    [".gt.", ".ge.", ".lt.", ".le.", ".eq."]);
                output = output + "_stratum." + stratTag(tabALL.k);
            end
            
            % get interaction p-value
            cfs = gather(mdl2.Coefficients);
            cfs_p = cfs{"x1:x"+opts.interactionIdx, "pValue"};
            cfs_p = compose("P_{int} = %.2g", string(cfs_p));

            genotag = [tabALL.A1(incnt)+tabALL.A1(incnt), ...
                tabALL.A1(incnt)+tabALL.A2(incnt),...
                tabALL.A2(incnt)+tabALL.A2(incnt)];
            interactionPlotter(gather(mdl2),...
                ["x1", "x"+opts.interactionIdx],...
                'labels', labels,...
                'trait', "Adjusted " + tabALL.Pheno(incnt), ...
                'ci', true, 'save', opts.saveplot, ...
                'res', opts.resolution, 'frmt', opts.format,...
                'savefig', opts.savefig, ...
                'xlab1', genotag, "output", output, ...
                "title", cfs_p);
        end
    end % variant loop
    
    if opts.firth && opts.binary
        opts.firthTime = sum(firthTime);
    end
end

for i = 1:size(bed, 2)
    tabALL = appendTab(tabALL, tabC{i}, tabALL.ct, tabALL.k, opts);
    tabALL.ct = tabALL.ct + 1;
end

end

%% ------------------------------------------------------------------------
function [tabALL, opts] = scoretest(DV, covariates, bed, opts, tabALL)

if isempty(covariates)
    fnan = false(numel(DV), 1);
else
    fnan = any(isnan(covariates), 2); % nan values in covariates
end
bed(fnan, :) = [];
covariates(fnan, :) = [];
% DVclean == DV

if ~isempty(opts.catCovar)
    nc = find(opts.catCovar);
    catcovar = covariates(:, nc);
    covariates(:, nc) = [];
    dcvar = cell(numel(nc), 1);
    for i = 1:numel(nc)
        dcvar{i} = dummyvar(categorical(catcovar(:, i)));
        dcvar{i}(:, 1) = [];
    end
    covariates = [covariates, horzcat(dcvar{:})];
end

% https://stackoverflow.com/questions/20562177/get-hat-matrix-from-qr-decomposition-for-weighted-least-square-regression
% https://stats.stackexchange.com/questions/139969/speeding-up-hat-matrices-like-xxx-1x-projection-matrices-and-other-as
% Regenie paper
Xcovar = [ones(size(covariates, 1), 1), covariates];
[N, p] = size(Xcovar);
S = struct;

if opts.firth && opts.binary
    S.batch = @(G) firthBatch(Xcovar,DV,G);   % returns [betaVec,seVec]
    S.CI    = @(a,b,se,dfe)[ b - norminv(1-a/2).*se , b + norminv(1-a/2).*se];

    % profile-PL p-value (same as logistf$prob)
    S.P  = @(LR) chi2pval(LR,1);

elseif opts.binary
    [B, ~, fit] = glmfit(Xcovar(:,2:end), DV, "binomial", "constant","on");
    
    mu = glmval(B, Xcovar(:,2:end), "logit", "constant","on");   % N×1
    w  = fit.wts;                                         % N×1
    tau = sqrt(w);                   % √W   (no zero because mu∈(0,1))

    Xw        = tau .* Xcovar;       % (√W)X
    [Qw,~]    = qr(Xw,0);            % N×p,   QwᵀQw = I
    
    Q1 = Qw ./ tau;                  % (√W)⁻¹ Qw     (N×p)
    Q2 = Qw .* tau;                  % (√W)    Qw    (N×p)
    
    % projection matrix P = X(XᵀWX)⁻¹XᵀW  realised as Q1*Q2ᵀ
    Gtilde = @(G) G - Q1*(Q2' * G);  % orthogonalised SNP(s), N×M

    ytilde = DV - mu;                           % score residuals, N×1

    S.U  = @(G) (Gtilde(G).' * ytilde).';         % 1×M
    S.V  = @(G) sum( (Gtilde(G).^2) .* w , 1 );   % 1×M   (always ≥0)
    
    S.beta = @(G)  S.U(G) ./ S.V(G);                              % Wald β̂  (1×M)
    S.se   = @(G)  sqrt(1 ./ S.V(G));                             % Wald SE  (1×M)
    S.P    = @(beta,se) chi2pval( (beta./se).^2 , 1 );          % p-value
    S.CI   = @(alpha, beta,se,dfe)[beta - norminv(1-alpha/2).*se ,  ...
                              beta + norminv(1-alpha/2).*se];

% ============================ LINEAR CASE =========================
else
    [Q,~]  = qr(Xcovar,0);
    ytilde = DV - Q*(Q'*DV);
    rss    = sum(ytilde.^2);
    sigma2 = rss / (N - p);           % residual variance

    S.Gtilde = @(G) G - Q*(Q' * G);                           % N×M
    Gnorm2   = @(G) sum( S.Gtilde(G).^2 , 1 );                % 1×M

    S.U     = @(G)           (S.Gtilde(G).' * ytilde).';      % 1×M
    S.V     = @(G) sigma2 .* Gnorm2(G);                       % 1×M
    S.beta  = @(G) S.U(G) ./ Gnorm2(G);                       % Wald β̂
    S.se    = @(G) sqrt( sigma2 ./ Gnorm2(G) );
    S.P  = @(beta, se) chi2pval( (beta./se).^2 , 1 );
    S.CI = @(alpha, beta, se, dfe)[beta - se.*tinv(1-alpha/2,dfe), beta + se.*tinv(1-alpha/2,dfe)];
end

S.dfe   = N - p - 2;      % +intercept +covars +SNP
clear stat

if ~opts.firth, clear Xcovar; end

% estimate parameters on each chunk. It seems MATLAB multithreading
% is faster than distributed calculations using parfor. I also
% avoid to implement this in a parfor loop because of memory
% overhead issues.
chunkIdx = opts.chunkIdx;

tab = struct('A2FREQ', [], 'p', [], 'p_adj', [], 'BETA', [],...
    'SE', [], 'CI', [], 'N', [], 'N_case', [], 'AF_case', [], ...
    'AF_control', []); % to be merged with tabALL columnwise or rowwise (for stratification)
for j = 1:numel(chunkIdx)-1 % loop over chunks of input bed for this chromosome
    idxRange = chunkIdx(j):chunkIdx(j+1)-1;
    bed_chunk = bed(:, idxRange);
    if opts.gpu, bed_chunk = gpuArray(bed_chunk); end

    if opts.firth && opts.binary
        [betaBlock, seBlock, LRblock] = S.batch(bed_chunk);
        tab.BETA = [tab.BETA ; betaBlock.' ];
        tab.SE   = [tab.SE   ; seBlock.'  ];
        tab.CI   = [tab.CI   ; S.CI(0.05,betaBlock.',seBlock.',S.dfe)];
        tab.p_adj= [tab.p_adj; S.P(LRblock.').'];
    else
        tab.BETA = [tab.BETA; S.beta(bed_chunk).'];
        tab.SE = [tab.SE; S.se(bed_chunk).'];
        tab.CI = [tab.CI; S.CI(0.05, tab.BETA(idxRange), tab.SE(idxRange), S.dfe)];
        tab.p_adj = [tab.p_adj; S.P(tab.BETA(idxRange), tab.SE(idxRange))];
    end



    tab.A2FREQ = [tab.A2FREQ; gather(mean(bed(:, idxRange))./opts.denom).']; % freq calculated here is independent of missing values in covariates
    if opts.binary
        tab.N_case = [tab.N_case; repmat(gather(sum(DV)), numel(tab.A2FREQ), 1)]; % effective case N
        tab.AF_case = [tab.AF_case; gather(mean(bed(DV == 1, idxRange), 1)./opts.denom).'];
        tab.AF_control = [tab.AF_control; gather(mean(bed(DV == 0, idxRange), 1)./opts.denom).'];
    end 

    if opts.verbose, progressGen(opts.progressBar, opts.strataCounter(j)); end % update progress bar 
end % end of loop over chunk of bed file

tab.N = repmat(size(bed, 1), size(bed, 2), 1); % N corresponds to effectize sample size used for adjusted model
tab.p = nan(numel(tab.p_adj), 1); % no need to estimate unadjusted p-value (use "glm" if desired)
if ~opts.binary
    [tab.N_case, tab.AF_case, tab.AF_control] = deal(nan(numel(tab.N), 1));
else
    tab.CI = exp(tab.CI); tab.BETA = exp(tab.BETA);
end

tab.CI = gather(tab.CI);
if opts.binary
    roundType = 'decimals';
else
    roundType = 'significant';
end
tab.CI = "[" + round(tab.CI(:, 1), opts.ciround, roundType) + ", " ...
    + round(tab.CI(:, 2), opts.ciround, roundType) + "]";

% merge tab structure with tabALL columnwise and rowwise (if
% stratification is requested)
fnames = fieldnames(tab);

for k = 1:size(bed, 2)
    for j = 1:numel(fnames)
        if ~isfield(tabALL, fnames{j}) % check if tabALL has these fields, if no, add 
            if strcmp(fnames{j}, 'CI')
                tabALL.(fnames{j}) = "";
            else
                tabALL.(fnames{j}) = [];
            end
        end
        tabALL.(fnames{j})(tabALL.ct, tabALL.k) = gather(tab.(fnames{j})(k));
    end
    tabALL.ct = tabALL.ct + 1;
end

end % END

%% ------------------------------------------------------------------------
function [beta, logLstar, ok] = firthLogit(X, y, wPrior, maxIter, tol, betaInit)
% Bias-reduced (Firth / Jeffreys) logistic regression.
%
% Required:
%   X         N×p  design (intercept included if wanted)
%   y         N×1  binary response
%
% Optional:
%   wPrior    N×1  prior weights              (default ones)
%   maxIter        max Newton iterations      (default 50)
%   tol            convergence tolerance      (default 1e-8)
%   betaInit  p×1  starting coefficients      (default zeros)
%
% Returns:
%   beta       p×1  coefficient vector
%   logLstar        penalised log-lik  ℓ + ½log|I|
%   ok              logical convergence flag

if nargin < 3 || isempty(wPrior),  wPrior  = ones(size(y)); end
if nargin < 4 || isempty(maxIter), maxIter = 50;            end
if nargin < 5 || isempty(tol),     tol      = 1e-8;         end
if nargin < 6 || isempty(betaInit)
    beta = zeros(size(X,2),1);
else
    beta = betaInit;
end

ok   = false;
epsW = 1e-12;                     % weight floor

for k = 1:maxIter
    eta = X*beta;
    mu  = 1 ./ (1 + exp(-eta));
    w   = max(wPrior .* mu .* (1-mu), epsW);

    % Weighted QR  (√W X = Q R)
    Xw     = bsxfun(@times, sqrt(w), X);
    [Q,R]  = qr(Xw,0);
    h      = sum(Q.^2, 2);                        % diag(H)

    % Firth-adjusted score
    Ustar  = X' * ( y - mu + 0.5*(1-2*mu).*h );

    % Newton step  (R'R) step = U*
    step   = R \ ( R' \ Ustar );
    beta   = beta + step;

    if norm(step,inf) < tol
        ok = true;
        break;
    end
end

logL     = sum( wPrior .* ( y.*log(mu) + (1-y).*log(1-mu) ) );
logDetI  = 2*sum( log( abs(diag(R)) ) );          % |I| = |R|²
logLstar = logL + 0.5*logDetI;

end % END

function [betaVec, seVec, LRvec] = firthBatch(X0, y, Gblock, wPrior)
% Firth profile-PL tests for a block of SNPs.
%   betaVec : 1×M  log-odds ratios
%   seVec   : 1×M  Wald SE (same shortcut logistf prints)
%   LRvec   : 1×M  penalised LR statistics (profile, 1 df)

if nargin<4,  wPrior = ones(size(y)); end

% ---------- null (covariates only) -----------------------------
[beta0, ~, ok0] = firthLogit(X0, y, wPrior);
if ~ok0,  error('Null Firth fit failed'); end

eta0 = X0*beta0;                         % fixed for every SNP
mu0  = 1./(1+exp(-eta0));
w0   = max(wPrior .* mu0 .* (1-mu0), 1e-12);

% log-lik part already in L0star; we still need |I_c| per SNP
logL0 = sum(wPrior .* (y.*log(mu0)+(1-y).*log(1-mu0)));

M        = size(Gblock,2);
betaVec  = NaN(1,M);
seVec    = NaN(1,M);
LRvec    = NaN(1,M);

for m = 1:M
    Gm = Gblock(:,m);
    if all(Gm==Gm(1)),  continue; end     % monomorphic

    % -------- full, unconstrained fit --------------------------
    X1         = [X0 , Gm];
    betaInit   = [beta0 ; 0];
    [b1, L1, ok1] = firthLogit(X1, y, wPrior, 40, 1e-8, betaInit);
    if ~ok1,   continue; end

    betaVec(m) = b1(end);

    % -------- Wald SE (same diag shortcut logistf uses) -------
    eta = X1*b1;  mu = 1./(1+exp(-eta));  w = max(wPrior.*mu.*(1-mu),1e-12);
    Xw = bsxfun(@times, sqrt(w), X1);
    [~,R] = qr(Xw,0);
    seVec(m) = 1/abs(R(end,end));

    % -------- penalised LR: full vs constrained (β_G = 0) ----
    % compute |I_c| with genotype column included but β_G = 0
    Xw_c   = bsxfun(@times, sqrt(w0), X1);  % note: w0 from null
    Rc     = triu(qr(Xw_c,0));
    logDetIc = 2*sum(log(abs(diag(Rc))));
    Lcstar   = logL0 + 0.5*logDetIc;

    LRvec(m) = 2*(L1 - Lcstar);            % 1 df
end

end % END

%% ------------------------------------------------------------------------
function [tabC, opts] = runModelsPar(bed, covariates, DV, DVclean, opts, parInfo)
% parallel version of runModels. 

% prepare progress text bar
if opts.verbose
    D = parallel.pool.DataQueue;
    parCounter = 1;
    progressBar = floor(linspace(1, size(bed, 2)+1, 11));
    afterEach(D, @parProgressBar);
    % afterEach(D, @(~)parProgressBar(F)); % with factor
end

interactStructOUT = parInfo.interactStruct;
tabEmpty = parInfo.tabEmpty;
denom = parInfo.denom;
spec1 = parInfo.spec1; spec2 = parInfo.spec2;
fnan = parInfo.fnan;

if opts.fitter == "scoretest"
     % scoretestFnc = parInfo.scoretestFnc; % DEPRECATED
     clear parInfo
end

tabC = ({}); firthTime = (0);

parfor i = 1:size(bed, 2) % loop over snps
    covariatesX = covariates;
    DVX = DV;
    DVXclean = DVclean;
    optsX = opts;

    interactStruct = interactStructOUT;
    snp = bed(:, i);

    if ~any(snp) || numel(DVXclean) < 2 % either empty variant or numel < 2
        tabC{i} = tabEmpty;
        continue
    end

    if optsX.gpu; snp = gpuArray(snp); end

    if optsX.binary
        if optsX.firth % perform Firth logistic regression (to deal with imbalanced/separation issues)
            Fopts = optsX;
            Fopts.snpIdx = 1; %Fopts.n = size(covariates, 2)-1;
            [mdl1, mdl2, firthTime(i)] = logistf(snp, covariatesX, DVX, Fopts, fnan); % R logistf removes nan values when necessary
            if ~isnan(optsX.interactionIdx)
                interactStruct.interactTerm1 = mdl1.interactTerm;
                interactStruct.interactTerm2 = mdl2.interactTerm;
            end

        else % binomial GLM
            if ~isnan(optsX.interactionIdx)
                mdl1 = compact(fitglm([snp, covariatesX(:, optsX.interactionIdx-1)], DVX, spec1, 'Distribution', 'binomial', 'Link', 'logit', 'CategoricalVars', optsX.catCovar1));
            else
                mdl1 = compact(fitglm(snp, DVX, spec1, 'Distribution', 'binomial', 'Link', 'logit', 'CategoricalVars', optsX.catCovar1));
            end
            mdl2 = compact(fitglm([snp(~fnan), covariatesX(~fnan, :)], DVXclean, spec2, 'Distribution', 'binomial', 'Link', 'logit', 'CategoricalVars', optsX.catCovar2));
        end
    else % continuous traits --------------------------------------------------
        if ~isnan(optsX.interactionIdx)
            mdl1 = compact(fitlm([snp, covariatesX(:, optsX.interactionIdx-1)], DVX, spec1, 'CategoricalVars', optsX.catCovar1));
        else
            mdl1 = compact(fitlm(snp, DVX, spec1, 'CategoricalVars', optsX.catCovar1));
        end
        mdl2 = compact(fitlm([snp(~fnan), covariatesX(~fnan, :)], DVXclean, spec2, 'CategoricalVars', optsX.catCovar2));

    end

    tab = struct; % faster than table
    tab.A2FREQ = gather(mean(snp)./denom); % freq calculated here is independent of missing values in covariates

    snp(fnan) = []; % effective sample in adjusted model
    tab.N = numel(snp); % N corresponds to effectize sample size used for adjusted model
    if optsX.binary
        tab.N_case = gather(sum(DVXclean)); % effective case N
        tab.AF_case = gather(mean(snp(DVXclean == 1))./denom);
        tab.AF_control = gather(mean(snp(DVXclean == 0))./denom);
    end

    if ~isnan(optsX.interactionIdx)
        tabC{i} = extractModelInfo(mdl1, mdl2, optsX, tab, interactStruct);
    else
        tabC{i} = extractModelInfo(mdl1, mdl2, optsX, tab);
    end
    if optsX.verbose; send(D, i); end % update progress bar 
end 

if opts.firth && opts.binary
    opts.firthTime = sum(firthTime);
end

function parProgressBar(~)
    progressTxt = [repmat('=', 1, sum(parCounter >= progressBar)),'>',...
         repmat(' ', 1, 10-sum(parCounter >= progressBar))];
     fprintf(repmat('\b', 1, 12))
     fprintf('%s]', progressTxt)
     parCounter = parCounter + 1;
end

end

%% ------------------------------------------------------------------------
function [out, opts] = call_regenie_in(cov, response, opts)

trait = opts.regenied.trait;
opts.regenied.trait = matlab.lang.makeValidName(opts.regenied.trait);
filenames.patt = getRandomName(opts.regenied.trait, 15);

if fileparts(opts.output) ~= ""
    if ~isfolder(fileparts(opts.output))
        mkdir(fileparts(opts.output))
    end
    absPath = string({dir(fileparts(opts.output)).folder});
    filenames.patt = fullfile(absPath(1), filenames.patt);
end

if isfield(opts.regenied, "phenoFile") % use already existing phenoFile
    filenames.pheno = opts.regenied.phenoFile;
    % new_covar_file = filenames.patt + ".covarFile.txt";
    % copyfile(opts.regenied.covarFile, new_covar_file)
    filenames.covar = opts.regenied.covarFile;
    % cmd = "awk 'BEGIN {FS=OFS=""\t""} {for (i=1; i<=NF; i++) {if (i>2) {$i=""C""(i-2)}}} 1' data.txt > temp && mv temp data.txt";

    % subset if 'qcinc' is provided'
    if isfield(opts, "qcinc") && ~isempty(opts.qcinc)
        response = readtable(filenames.pheno);
        response(~ismember(response.FID, opts.qcinc), :) = [];
        cov = readtable(filenames.covar);
        cov(~ismember(cov.FID, opts.qcinc), :) = [];
        filenames.pheno = filenames.patt + ".phenoFile.txt";
        writetable(response, filenames.pheno, 'Delimiter', '\t')
        filenames.covar = filenames.patt + ".covarFile.txt";
        writetable(cov, filenames.covar, 'Delimiter', '\t');
        
        dfe = height(cov) - 2 - width(cov);
        clear response cov
    else
        dfe = bfilereader(filenames.pheno, "summary", "linecount", "verbose", "off");
        dfe = dfe.linecount - 3;
        [~, cdfe] = system("wsl head -n 1 " + makeWSLpath(filenames.covar));
        cdfe = string(split(cdfe)); cdfe(cdfe == "") =[];
        dfe = dfe - numel(cdfe) - 2; % IID and FID
    end
    covFlag = true;

else
    eid = opts.qc_eid;
    response = gather(response);
    cov = gather(cov);
    fnan = any(isnan(cov), 2);
    response(fnan, :) = [];
    dfe = size(response, 1) - 2; % degree of freedom for 95% CI 
    if ~isempty(cov)
        cov(fnan, :) = [];
        dfe = dfe - size(cov, 2);
    end
    eid(fnan) = [];
end

% check binary input (firth also?) and calculate N
opts.firth_se = false;
opts.approx = false; %@07NOV2024
opts.spa = true;
if opts.binary
    opts.bt = true;
    if opts.firth
        opts.firth_se = true;
        opts.approx = true;
        opts.spa = false; % either SPA or Firth
    end
else
    opts.bt = false;
    opts.spa = false;
end

if ~strcmp(opts.transform, "none") && ~opts.skiptransform && ~isfield(filenames, "pheno")
    response = feval(opts.transform, gather(response));
end

% write response to pheno file


if ~isfield(filenames, "pheno")
    response = array2table(response, 'VariableNames', opts.regenied.trait);
    [response.FID, response.IID] = deal(eid);
    clear eid
    response = movevars(response, {'FID', 'IID'}, 'Before', opts.regenied.trait);
    
    filenames.pheno = filenames.patt + ".phenoFile.txt";
    writetable(response, filenames.pheno, 'Delimiter', '\t')
    
    response(:, 3:end) = []; % keep it for covariates
    
    if ~isempty(cov)
        cov = array2table(cov, 'VariableNames', "C" + (1:size(cov, 2)));
        response = [response, cov];
        filenames.covar = filenames.patt + ".covarFile.txt";
        writetable(response, filenames.covar, 'Delimiter', '\t');
        clear response cov
        covFlag = true;
    else
        covFlag = false;
    end
end

% check 'condition' option
condids = "";
if ~isempty(opts.condition)
    if opts.wesin % gene-based tests
        idx = ismember(opts.wesinfo.gene, opts.condition);
        if ~any(idx) % gene id?
            idx = ismember(opts.wesinfo.ENSG, opts.condition);
        end

        if any(idx)
            condids = opts.wesinfo.id(idx);
        else
            condids = opts.condition;
        end

    else % SVA
        condids = opts.condition;
    end
end
if isempty(condids), condids = ""; end

% write required files for either association or gene-based tests
if opts.wesin 
    if isfield(opts.regenied.geno, "anno_file") % use existing (geneWrapper) anno/mask/set files
        filenames.mask = opts.regenied.geno.mask_def;
        filenames.ann = opts.regenied.geno.anno_file;
        filenames.set = opts.regenied.geno.set_list;
    else
        % create mask/set/annotation files        
        wesinfotmp = opts.wesinfo;
        [~, fidx] = ismember(opts.regenied.geno.bim.snp, wesinfotmp.id);
        fidx(fidx < 1) = []; wesinfotmp = wesinfotmp(fidx, :);
        
        if ~opts.regenied.geno.multi
            gidx = [1, numel(opts.regenied.geno.bim.snp)];
        else
            gidx = opts.regenied.geno.bim.idx;
        end
        
        for j = 1:size(gidx, 1) % loop over genes
            gene = wesinfotmp.gene(gidx(j, 1));
            if ismember(gene, ["", "-"]), gene = "NA"; end
            if ~ismember(wesinfotmp.ENSG(gidx(j,1)), ["", "-"]), gene = gene + "_" + wesinfotmp.ENSG(gidx(j,1)); end
            varids = opts.regenied.geno.bim.snp(gidx(j, 1):gidx(j, 2));
            varids = setdiff(varids, condids); % same as opts.wesinfo.id
            annfile = compose("%s\t%s\t%s", varids, ...
                repmat(gene, numel(varids), 1), repmat("lof", numel(varids), 1));
            setfile = compose("%s\t%s\t%d\t%s", gene, ...
                string(opts.regenied.geno.bim.chr(gidx(j, 1))), 1, join(varids, ","));

            filenames.mask(j,1) = filenames.patt + j + ".mask.txt";
            writematrix(compose("%s\t%s", "m1", "lof"), filenames.mask(j), ...
                'QuoteStrings', false)
        
            filenames.ann(j,1) = filenames.patt + j + ".ann.txt";
            filenames.set(j,1) = filenames.patt + j + ".set.txt";
            writematrix(annfile, filenames.ann(j), 'QuoteStrings', false)
            writematrix(setfile, filenames.set(j), 'QuoteStrings', false)
        end
    end
else
    filenames.extract = filenames.patt + ".extract.txt";
    if isrow(opts.regenied.snp), opts.regenied.snp = opts.regenied.snp'; end
    writematrix(opts.regenied.snp, filenames.extract, 'QuoteStrings', false)
end

if ~isempty(opts.condition) && all(condids ~= "")
    filenames.condition = filenames.patt + ".condition.txt";
    writematrix(condids, filenames.condition, 'QuoteStrings', false)
    ropts.condition_list = filenames.condition;
end

% if opts.wesin % gene-based tests ------------------------------------------

ropts.step = 2;
if isfield(opts, "pred") && ~ismissing(opts.regenied.pred) && opts.regenied.pred ~= ""
    ropts.pred = opts.regenied.pred;

    %@07NOV2024: find null firth file 
    pred_path = fileparts(ropts.pred);
    if pred_path == "", pred_path = pwd; end
    firth_null_file = getfilenames(pred_path, "list", fullpath=true).list;
    idx_firth_null = firth_null_file.endsWith("_firth.list");
    if any(idx_firth_null) && opts.firth
        ropts.use_null_firth = firth_null_file(idx_firth_null);
    end

end
ropts.phenoFile = filenames.pheno;
ropts.minINFO = opts.infoscore;
ropts.bgenhome = opts.bgenhome;
ropts.bt = opts.bt;
ropts.firth = opts.firth;
ropts.spa = opts.spa;
ropts.firth_se = opts.firth_se;
ropts.approx = opts.approx;
ropts.debug = true;
ropts.phenopath = opts.phenopath;
ropts.out = filenames.patt;
ropts.verbose = false;
ropts.rare_mac = opts.rare_mac; % @12MARCH2023: for interaction (< this, HLM is used, above this, robust SE)
ropts.force_robust = opts.force_robust; % @12MARCH2023: for interaction (if ture, only robust SE is used)
ropts.write_mask = opts.write_mask;
if ~isnan(opts.workers), ropts.threads = opts.workers; end

if opts.interaction ~= ""
    if isfield(opts, "covarFile") %@29APR2024: use interaction pheno name directly
        [opts.interactiontag, ropts.interaction] = deal(opts.interaction);
    else
        ropts.interaction = "C" + (opts.interactionIdx - 1);
    end
    ropts.spa = false; % REGENIE --interaction conflicts with SPA or Firth
    ropts.firth = false;
end

%@16OCT2024: support for categorical covariates
if isfield(opts.regenied, "catCovar") % @26NOV2024: categorical covariates directly in the covar file
    ropts.catCovarList = opts.regenied.catCovar;
elseif isfield(opts, "catCovar") && any(opts.catCovar)
    ropts.catCovarList = join("C" + find(opts.catCovar), ",");
end

if opts.wesin
    ropts.anno_file = filenames.ann; ropts.set_list = filenames.set;
    ropts.mask_def = filenames.mask; ropts.aaf_bins = opts.aaf_bins;
    
    ropts.build_mask = opts.build_mask;
    ropts.skat_params = opts.skat_params;
    if isfield(opts, "vc_maxAAF")
        ropts.vc_maxAAF = opts.raremaf;
        ropts.vc_maxAAF = opts.vc_maxAAF;
    end
    ropts.weshome = opts.weshome;
else
    % check if sex chromosomes are queried. In this case "pgen" files
    % should be used: https://github.com/rgcgithub/regenie/issues/195
    if endsWith(opts.regenied.bgenfile, ".bed")
        ropts.bed = opts.regenied.bgenfile;
    elseif any(contains(lower(opts.regenied.bgenfile), reshape(["_c", "_chr"]' + ["x", "y", "xy"] + "_", 1, [])))
        ropts.pgen = regexprep(opts.regenied.bgenfile, ".bgen$", ".pgen");
    else
        ropts.bgen = opts.regenied.bgenfile; 
        ropts.sample = regexprep(opts.regenied.bgenfile, ".bgen$", ".sample");
    end
    
    ropts.extract = filenames.extract;
end
    
if opts.wesin && isfield(opts.regenied.geno, 'bfile') % use merged file for multi-gene mode
    ropts.bed = opts.regenied.geno.bfile;
end

if covFlag
    ropts.covarFile = filenames.covar;
end
ropts = namedargs2cell(ropts);
call_regenie(ropts{:});

% files to be removed in the end
if any(fileparts(filenames.patt) == "")
    rmfiles = getfilenames(pwd, "*", fullpath=false).x_;
else
    rmfiles = getfilenames(fileparts(filenames.patt), "*", fullpath=true).x_;
end
rmfiles(~rmfiles.startsWith(filenames.patt)) = [];

% if write_mask flag is on there are bed files for each pheno
if any(rmfiles.endsWith(".bed"))
    to_bed_path = fileparts2(filenames.patt, 2);
    these_plink_files = rmfiles(rmfiles.endsWith([".bed", ".bim", ".fam"]));
    [~, ~, bext] = fileparts(these_plink_files);
    to_bed_path = fullfile(to_bed_path, opts.regenied.trait + bext);
    arrayfun(@(x,y)movefile(x, y), these_plink_files, to_bed_path)
end

logfile = rmfiles(endsWith(rmfiles, ".log"));
logfile = readlines(logfile);

if opts.wesin
    % read snplists
    opts.regenied.snplist = readtable(filenames.patt + ".snplist.txt", 'Delimiter', "\t", ...
        'ReadVariableNames', false, 'TextType', 'string');
    if width(opts.regenied.snplist) < 2 % no variants remained for REGENIE
        out = table;
        if ~opts.debug % @30MARCH2025:otherwise keep all intermediate files
            delete(rmfiles{:})
        end
        return
    end

    try % REGENIE v3.5 onwards
        opts.regenied.snplist.Properties.VariableNames = {'mask', 'CHR', 'POS', 'variants'};
    catch
        opts.regenied.snplist.Properties.VariableNames = {'mask', 'variants'};
    end
    opts.regenied.snplist.Pheno = repmat(trait, height(opts.regenied.snplist), 1);
    opts.regenied.snplist.Pheno(2:end) = "";
    opts.regenied.snplist = movevars(opts.regenied.snplist, 'Pheno', 'Before', 1);
end

% read and modify summary stat ----------------------------------------
if opts.wesin
    out = load(filenames.patt + ".sumstat.step2.mat").res;
    out = vertcat(out{:});
    out(:, ["GENPOS", "CHISQ", "LOG10P"]) = [];
    out.ALLELE1 = out.ALLELE1 + "." + out.TEST; out.TEST = [];
    out.ALLELE0 = [];
    if opts.regenied.geno.multi
        out.ID(:) = join(unique(opts.wesinfo.gene), ",");
    end
else
    out = readtable(filenames.patt + ".sumstat.step2.txt", ...
        'VariableNamingRule', 'preserve', 'ReadVariableNames', true, ...
        'TextType', 'string');
    if isempty(out)
        if ~opts.debug % @30MARCH2025:otherwise keep all intermediate files
            delete(rmfiles{:})
        end
        return
    end
end

out.Pheno = repmat(trait, height(out), 1);
% if opts.wesin, out.Pheno(2:end) = ""; end
out = movevars(out, "Pheno", 'Before', 1);

% make table 'out' consistent over binary/continuous traits
bcol = find(strcmpi(out.Properties.VariableNames, 'beta'));
out.Properties.VariableNames(bcol) = "β|OR(SPA)";
if opts.firth, out.Properties.VariableNames(bcol) = "β|OR(Firth)"; end

% add 95% CI
if ~opts.wesin
    out.CI = [out.(bcol) - out.SE.*tinv(1-0.05/2, dfe),...
        out.(bcol) + out.SE.*tinv(1-0.05/2, dfe)];
end

if opts.binary % beta -> OR
    out.(bcol) = exp(out.(bcol));
    if ~opts.wesin, out.CI = exp(out.CI); end
end

if ~opts.wesin
    if opts.binary
        roundType = 'decimals';
    else
        roundType = 'significant';
    end
    out.CI = "[" + round(out.CI(:, 1), opts.ciround, roundType) + ", " ...
        + round(out.CI(:, 2), opts.ciround, roundType) + "]";
end
bcol = startsWith(out.Properties.VariableNames, 'β|OR');
bcol = string(out.Properties.VariableNames{bcol});
        
if ~any(strcmp(out.Properties.VariableNames, 'A1FREQ_CASES'))
    [out.A1FREQ_CASES, out.A1FREQ_CONTROLS] = deal(nan(height(out), 1));
end

% rename columns
if opts.wesin
    out = renamevars(out, ...
        ["A1FREQ_CASES", "A1FREQ_CONTROLS", "ALLELE1", "CHROM"], ...
        ["A1F.case", "A1F.control", "A1", "CHR"]);
else
    rmcols = intersect(["TEST", "INFO", "CHISQ", "EXTRA"], out.Properties.VariableNames);
    if any(isnan(opts.interactionIdx))
        out(:, rmcols) = [];
    else % interaction
        out(:, setdiff(rmcols, "TEST")) = [];
        out.TEST = regexprep(out.TEST, '(?<=.*-INT_SNPx)(.*)', opts.interactiontag);
    end
    pcol = startsWith(out.Properties.VariableNames, "LOG10P");
    pcol = string(out.Properties.VariableNames(pcol));
    out.(pcol) = 10.^-out.(pcol);

    out = renamevars(out, ...
        ["GENPOS", "ID", "A1FREQ_CASES", "A1FREQ_CONTROLS", "ALLELE0", "ALLELE1", "A1FREQ", "CHROM", pcol, "CI", bcol], ...
        ["POS", "SNP", "AF_case", "AF_control", "A1", "A2", "A2FREQ", "CHR", "p_adj", "CI", "BETA"]);
    out.p = nan(height(out), 1);
end

% add N.case/N.control columns
n.case = nan; n.control = nan; n.total = unique(out.N(~isnan(out.N)));
if opts.binary
    cc = regexp(logfile, "(\d+) cases and (\d+) controls", 'tokens');
    cc(cellfun(@isempty, cc)) = [];
    if ~isempty(cc)
        n.case = double(cc{1}{1}(1));
        n.control = double(cc{1}{1}(2));
    end
end

if opts.wesin
    out.("N.case") = n.case.*ones(height(out), 1);
    out = convertvars(out, ["A1F.case", "A1F.control"], "double");
    out.("MAC.case") = round(n.case*out.("A1F.case")*2).*ones(height(out), 1);
    out.("MAC.control") = round(n.control*out.("A1F.control")*2).*ones(height(out), 1);
else
    out.N_case = n.case.*ones(height(out), 1);
    out = convertvars(out, ["AF_case", "AF_control"], "double");
end

if opts.wesin
    out = out(:, ["Pheno", "ID", "A1", bcol, "SE", "P", "A1FREQ", ...
        "A1F.case", "A1F.control", "N", "N.case", "MAC.case", ...
        "MAC.control", "CHR", "EXTRA"]);
elseif any(isnan(opts.interactionIdx))
    out = out(:, ["Pheno", "SNP", "CHR", "POS", "A1", "A2", "A2FREQ", ...
        "p", "p_adj", "BETA", "SE", "CI", "N", "N_case", ...
        "AF_case", "AF_control"]);
else
    out = out(:, ["Pheno", "SNP", "CHR", "POS", "A1", "A2", "A2FREQ", ...
        "p", "p_adj", "BETA", "SE", "CI", "N", "N_case", ...
        "AF_case", "AF_control", "TEST"]);
end

% check SNP tags
if isfield(opts.regenied, "atag") && all(opts.regenied.atag ~= "")
    [f1, f2] = ismember(out.SNP, opts.regenied.asnp);
    out.SNP(f1) = opts.regenied.atag(f2(f1));
end
% ---------------------------------------------------------------------

if opts.verbose && opts.wesin
    disp(strjoin(logfile, '\n'))
end

if ~opts.debug % @30MARCH2025:otherwise keep all intermediate files
    delete(rmfiles{:})
end

end

%% ------------------------------------------------------------------------
function tab = extractModelInfo(mdl1, mdl2, opts, tabS, interactStruct)
% a helper subfunction for runModels function, extracting required
% information for report from GLMs. Serves both parallel and
% single-threaded functions.

coefCI = mdl2.coefCI;
coefCI = gather(coefCI(2, :));
beta = gather(mdl2.Coefficients.Estimate(2));

if opts.binary
    coefCI = exp(coefCI);
    beta = exp(beta);
end

tab = struct; % faster than table
tab.A2FREQ = tabS.A2FREQ;% freq calculated here is independent of missing values in covariates

if ~isnan(opts.interactionIdx)
    idx = startsWith(mdl1.CoefficientNames, interactStruct.interactTerm1);
    tab.p_int = gather(mdl1.Coefficients.pValue(idx));
    idx = startsWith(mdl2.CoefficientNames, interactStruct.interactTerm2);
    tab.p_int_adj = gather(mdl2.Coefficients.pValue(idx));
    tab.BETA_int = gather(mdl2.Coefficients.Estimate(idx));
    tab.SE_int = gather(mdl2.Coefficients.SE(idx));

    % add interaction summary stats for the second variable SNPxVAR
    idx = startsWith(mdl1.CoefficientNames, extractAfter(interactStruct.interactTerm1, ":"));
    tab.p_var = gather(mdl1.Coefficients.pValue(idx));
    idx = startsWith(mdl2.CoefficientNames, extractAfter(interactStruct.interactTerm2, ":"));
    tab.p_var_adj = gather(mdl2.Coefficients.pValue(idx));
    tab.BETA_var = gather(mdl2.Coefficients.Estimate(idx));
    tab.SE_var = gather(mdl2.Coefficients.SE(idx));

    if opts.binary
        tab.BETA_int = exp(tab.BETA_int);
        tab.BETA_var = exp(tab.BETA_var);
    end
end

tab.p = gather(mdl1.Coefficients.pValue(2));
tab.p_adj = gather(mdl2.Coefficients.pValue(2));
tab.BETA = beta;
if opts.firth && opts.binary
    if opts.flic
        tab.p_flic = mdl1.flic.Coefficients.pValue(2);
        tab.p_adj_flic = mdl2.flic.Coefficients.pValue(2);
        tab.BETA_flic = exp(mdl2.flic.Coefficients.Estimate(2));
    end
    if opts.flac
        tab.p_flac = mdl1.flac.Coefficients.pValue(2);
        tab.p_adj_flac = mdl2.flac.Coefficients.pValue(2);
        tab.BETA_flac = exp(mdl2.flac.Coefficients.Estimate(2));
    end
end
tab.SE = gather(mdl2.Coefficients.SE(2));
if opts.binary
    roundType = 'decimals';
else
    roundType = 'significant';
end
tab.CI = "[" + round(coefCI(1), opts.ciround, roundType) + ", " ...
    + round(coefCI(2), opts.ciround, roundType) + "]";

tab.N = tabS.N; % N corresponds to effectize sample size used for adjusted model
if opts.binary
    tab.N_case = tabS.N_case; % effective case N
    tab.AF_case = tabS.AF_case;
    tab.AF_control = tabS.AF_control;
end
end

%% ------------------------------------------------------------------------
function [mdl1, mdl2, ti] = logistf(snp, covariates, DV, opts, fnan)
% runs Firth's Logistic Regression with R package 'logistf'
% It also returns Firth's Logistic Regression With Intercept Correction
% (FLIC) and FLAC - Firth's Logistic Regression With Added Covariate

% NOTE: When rare variants are tested for association with a highly
% unbalanced traits (i.e. a trait that has low sample prevalence),
% quasi-complete separation can occur in logistic regression where none of
% the cases contain the minor allele which leads to very unstable estimates
% of the effect sizes and standard errors. One solution to address this is
% with Firth correction, where a penalty based on the observed Fisher
% information matrix is added in the log-likelihood. Taken from:
% https://www.biorxiv.org/content/10.1101/2020.06.19.162354v2.full

% TO DO: data pruning should be performed within R, also regressions for
% chunks of bed files must be performed within R since now analyses are
% single variant oriented, which is suboptimal for large number of
% variants.

snp = gather(snp); covariates = gather(covariates); DV = gather(DV);

[skipModel1, skipModel2] = deal(false);
if any(fnan)
    if numel(unique(snp)) < 2 % variance(snp) == 0
        disp([newline, 'var(snp) is 0, skipping Firth''s logistic regression.'])
        [skipModel1, skipModel2] = deal(true);
    else
        skipModel1 = false;
        if numel(unique(snp(~fnan))) < 2 
            disp([newline, 'var(snp) is 0, skipping adjusted Firth''s logistic regression.'])
            skipModel2 = true;
        else
            skipModel2 = false;
        end
    end
end

ncovar = size(covariates, 2); nsnp = size(snp, 2);
if ~isnan(opts.interactionIdx) % check interaction term
    spec1 = "y ~"  + join(["x"+opts.snpIdx, "c"+(opts.interactionIdx-1)], '*');
    covarlbls = "c"+(1:ncovar);
    covarlbls(covarlbls  == "c"+(opts.interactionIdx-1)) = [];
    covarlbls = join(covarlbls, '+');
    spec2 = spec1 + "+" + covarlbls;
    interactTerm = join(["x"+opts.snpIdx, "c"+(opts.interactionIdx-1)], ':');
else
    spec1 = "y ~"  + "x"+opts.snpIdx;
    spec2 = "y ~ " + "x"+opts.snpIdx + "+" + join("c"+(1:ncovar), '+');
end

if skipModel1 && skipModel2 % variance of snp is zero, no need to perform regression analysis.
    mdl1 = toMatlabGlm("logistf.matlab.RETURNME", false, spec1);
    mdl2 = toMatlabGlm("logistf.matlab.RETURNME", false, spec2);
    if opts.flic
        mdl1.flic = toMatlabGlm("logistf.matlab.RETURNME", false, spec1);
        mdl2.flic = toMatlabGlm("logistf.matlab.RETURNME", false, spec2);
    end
    if opts.flac
        mdl1.flac = toMatlabGlm("logistf.matlab.RETURNME", false, spec1);
        mdl2.flac = toMatlabGlm("logistf.matlab.RETURNME", false, spec2);
    end
    ti = 0;
    return
end

% generate random string. This is useful for parfeval call (asynchronously)
% or parallel (parfor) calling of gwasrunner 
charSet = ['a':'z',upper('a':'z'),'0':'9'];
randStr = string(charSet(randi(numel(charSet),1, randi([9 15], 1))));
maxloop = 1;
while exist("snp.logistf."+randStr+".mat", 'file')
    randStr = string(charSet(randi(numel(charSet),1, randi([9 15], 1))));
    maxloop = maxloop + 1;
    if maxloop > 5
        error('cannot generate unique random string for logistf package :(')
    end
end
save("snp.logistf."+randStr+".mat", 'snp')
save("covariates.logistf."+randStr+".mat", 'covariates')
save("DV.logistf."+randStr+".mat", 'DV')
clear snp covariates DV

currDir = strrep(pwd, '\', '/');
Rcode = string;
Rcode(1, 1) = "setwd('" + currDir + "')";
Rcode(2, 1) = "library(logistf)";
Rcode(3, 1) = "data <- data.frame(R.matlab::readMat('snp.logistf."+randStr+".mat')$snp,"+ ...
    " R.matlab::readMat('covariates.logistf."+randStr+".mat')$covariates,"+ ...
    " R.matlab::readMat('DV.logistf."+randStr+".mat')$DV)";
Rcode(4, 1) = "colnames(data)<-c(paste0('x', 1:"+nsnp+...
    "),paste0('c', 1:"+ncovar+"), 'y')";
% Rcode(3, 1) = "data = read.table('" + data + "',sep = '\t', header = T)";
Rcode(5, 1) = "lf <- logistf(formula = " + spec1 + ", data = data)";
Rcode(6, 1) = "r = data.frame(beta = lf$coefficients, se = sqrt(diag(vcov(lf))), p = lf$prob, ci_l = lf$ci.lower, ci_u = lf$ci.upper)";
Rcode(7, 1) = "write.table(r, 'logistf.matlab1"+randStr+".txt', sep = '\t', quote = F, row.names = T)";

if opts.flic
    Rcode(8, 1) = "lf1 <- flic(lf)";
    Rcode(9, 1) = "r = data.frame(beta = lf1$coefficients, se = sqrt(diag(lf1$var)), p = lf1$prob, ci_l = lf1$ci.lower, ci_u = lf1$ci.upper)";
    Rcode(10, 1) = "write.table(r, 'logistf.flic.matlab1"+randStr+".txt', sep = '\t', quote = F, row.names = T)";
end

if opts.flac
    Rcode(11, 1) = "lf2 <- flac(lf)";
    Rcode(12, 1) = "r = data.frame(beta = lf2$coefficients, se = sqrt(diag(lf2$var)), p = lf2$prob, ci_l = lf2$ci.lower, ci_u = lf2$ci.upper)";
    Rcode(13, 1) = "write.table(r, 'logistf.flac.matlab1"+randStr+".txt', sep = '\t', quote = F, row.names = T)";
end

if ~skipModel2
    Rcode(14, 1) = "lf <- logistf(formula = " + spec2 + ", data = data)";
    Rcode(15, 1) = "r = data.frame(beta = lf$coefficients, se = sqrt(diag(vcov(lf))), p = lf$prob, ci_l = lf$ci.lower, ci_u = lf$ci.upper)";
    Rcode(16, 1) = "write.table(r, 'logistf.matlab2"+randStr+".txt', sep = '\t', quote = F, row.names = T)";

    if opts.flic
        Rcode(17, 1) = "lf1 <- flic(lf)";
        Rcode(18, 1) = "r = data.frame(beta = lf1$coefficients, se = sqrt(diag(lf1$var)), p = lf1$prob, ci_l = lf1$ci.lower, ci_u = lf1$ci.upper)";
        Rcode(19, 1) = "write.table(r, 'logistf.flic.matlab2"+randStr+".txt', sep = '\t', quote = F, row.names = T)";
    end

    if opts.flac
        Rcode(20, 1) = "lf2 <- flac(lf)";
        Rcode(21, 1) = "r = data.frame(beta = lf2$coefficients, se = sqrt(diag(lf2$var)), p = lf2$prob, ci_l = lf2$ci.lower, ci_u = lf2$ci.upper)";
        Rcode(22, 1) = "write.table(r, 'logistf.flac.matlab2"+randStr+".txt', sep = '\t', quote = F, row.names = T)";
    end
end

writematrix(Rcode, "gwasrunner_firth_"+randStr+".r", 'QuoteStrings', false, 'FileType','text')
% asynchronous: https://www.mathworks.com/matlabcentral/answers/339952-force-matlab-to-continue-when-stuck-in-busy-mode#answer_267135
timeM2R = MATLAB2Rconnector("gwasrunner_firth_"+randStr+".r", 'delr', true, 'maxTime', opts.maxTime);
if opts.verbose; tic; end
delete("snp.logistf."+randStr+".mat")
delete("covariates.logistf."+randStr+".mat")
delete("DV.logistf."+randStr+".mat")

% depends on error flag, can generate an empty struct
if skipModel2; opts.error = false; end
checkLofistfFiles = ["logistf.matlab1"+randStr+".txt", "logistf.matlab2"+randStr+".txt"];
if ~all(arrayfun(@(x)exist(x, 'file'), checkLofistfFiles))
    if opts.error
        error('gwasrunner:logistf', 'Firth''s regression failed!')
    end
end
mdl1 = toMatlabGlm("logistf.matlab1"+randStr+".txt", opts.error, spec1);
mdl2 = toMatlabGlm("logistf.matlab2"+randStr+".txt", opts.error, spec2);
if opts.flic
    mdl1.flic = toMatlabGlm("logistf.flic.matlab1"+randStr+".txt", opts.error, spec1);
    mdl2.flic = toMatlabGlm("logistf.flic.matlab2"+randStr+".txt", opts.error, spec2);
end
if opts.flac
    mdl1.flac = toMatlabGlm("logistf.flac.matlab1"+randStr+".txt", opts.error, spec1);
    mdl2.flac = toMatlabGlm("logistf.flac.matlab2"+randStr+".txt", opts.error, spec2);
end

if ~isnan(opts.interactionIdx)
    mdl1.interactTerm = interactTerm;
    mdl2.interactTerm = interactTerm;
end

if opts.verbose
    ti = toc + timeM2R;
else
    ti = 0;
end
end 

%% ------------------------------------------------------------------------
function mdl = toMatlabGlm(rfile, errorFlag, spec)
vars = {'Estimate', 'SE', 'pValue', 'coefCI_l', 'coefCI_u'};
if ~errorFlag && ~exist(rfile, 'file')
    spec = split(spec, '~');
    spec = cellstr(split(strtrim(spec(2)), '+'));
    idxInteract = contains(spec, '*'); 
    specInteract = spec(idxInteract);
    specInteract = split(specInteract, '*');
    specInteract = [specInteract; join(specInteract, ':')];
    spec(idxInteract) = [];
    spec = vertcat('(Intercept)', specInteract, spec);
    mdl.CoefficientNames = spec;
    mdl.coefCI = nan(numel(spec), 2);
    [mdl.Coefficients.SE, mdl.Coefficients.Estimate, ...
        mdl.Coefficients.pValue, mdl.Coefficients.coefCI_l, ...
        mdl.Coefficients.coefCI_u] = deal(nan(numel(spec), 1));
    return
end
mdls = readtable(rfile, 'ReadRowNames', true, 'ReadVariableNames', false);
mdls.Properties.VariableNames = vars;
CoefficientNames = mdls.Row;
mdls = table2struct(mdls, 'ToScalar', true);
mdl.Coefficients = mdls;
mdl.coefCI = [mdl.Coefficients.coefCI_l, mdl.Coefficients.coefCI_u];
mdl.CoefficientNames = CoefficientNames;
delete(rfile)
end

%% ------------------------------------------------------------------------
function tab = toCleanTab(tab, opts)
tab = struct2table(tab);
dblvars = ["A2FREQ", "p", "p_adj", "BETA", "SE"];
if ~isnan(opts.interactionIdx) && ~opts.regenie
    dblvars = [dblvars, "p_int", "p_int_adj", "p_var", "p_var_adj", ...
        reshape(["BETA_", "SE_"] + ["int", "var"]', 1, [])];
end
if opts.binary
    dblvars = [dblvars, "AF_case", "AF_control"];
end
if opts.firth && opts.binary
    if opts.flic
        dblvars = [dblvars, "p_flic", "p_adj_flic", "BETA_flic"];
    end
    if opts.flac
        dblvars = [dblvars, "p_flac", "p_adj_flac", "BETA_flac"];
    end
end
% round to 3 significant digits
tab_double = varfun(@(x)round(x, opts.betaround, 'significant'), tab,...
    'InputVariables', dblvars);
for i = 1:numel(dblvars)
    tab.(dblvars{i}) = tab_double.(i);
end
clear tabALL_double

tab = renamevars(tab, ["N_case", "AF_case", "AF_control"], ...
    ["N.case", "AF.case", "AF.control"]);

if opts.binary
    tab = renamevars(tab, "BETA", "OR");
    if opts.firth
        tab = renamevars(tab, "OR", "OR.Firth");
        if opts.flic
            tab = renamevars(tab, "BETA_flic", "OR.FLIC");
        end
        if opts.flac
            tab = renamevars(tab, "BETA_flac", "OR.FLAC");
        end
    end
else
    tab = renamevars(tab, "BETA", "β");
end

if opts.firth && opts.binary
    tab = renamevars(tab, ["p", "p_adj"], ["P.Firth", "P-adj.Firth"]);
    oldVars = tab.Properties.VariableNames;
    if opts.flic
        tab.Properties.VariableNames(ismember(oldVars, 'p_flic')) = {'P(FLIC)'};
        tab.Properties.VariableNames(ismember(oldVars, 'p_adj_flic')) = {'P-adj(FLIC)'};
    end
    if opts.flac
        tab.Properties.VariableNames(ismember(oldVars, 'p_flac')) = {'P(FLAC)'};
        tab.Properties.VariableNames(ismember(oldVars, 'p_adj_flac')) = {'P-adj(FLAC)'};
    end
else
    tab = renamevars(tab, ["p", "p_adj"], ["P", "P-adj"]);
end

if ~isnan(opts.interactionIdx) && ~opts.regenie
    n = height(tab);
    tab = stack(tab, ["P", "p_int", "p_var"], "NewDataVariableName", "P", "IndexVariableName", "TEST");
    tab.TEST = repmat([upper(opts.model); "SNPx" + opts.interactiontag; opts.interactiontag], n, 1);
    idx = 1:3:height(tab);
    tab.("P-adj") = reshape([tab.("P-adj")(idx), tab.p_int_adj(idx), tab.p_var_adj(idx)]', height(tab), 1);

    if opts.binary, beta = "OR"; else, beta = "β"; end
    tab.(beta) = reshape([tab.(beta)(idx), tab.BETA_int(idx), tab.BETA_var(idx)]', height(tab), 1);
    tab.SE = reshape([tab.SE(idx), tab.SE_int(idx), tab.SE_var(idx)]', height(tab), 1);
    tab.CI(setdiff(1:height(tab), idx)) = missing;
    tab(:, ["p_int_adj", "p_var_adj", "BETA_int", "BETA_var", "SE_int", "SE_var"]) = [];
    tab = movevars(tab, ["TEST", "P"], "Before", "P-adj");
end

tab = renamevars(tab, "CI", "95% CI");
if ~any(isnan(double(tab.CHR)))
    tab.CHR = double(tab.CHR);
end
tab.Pheno(2:end) = missing;
end

%% ------------------------------------------------------------------------
function effPlotter(tab, opts)
% a modified version of Print2ExcellPlotter function 
% close all gcf
tab = toCleanTab(tab, opts);

if string(opts.output) == ""
    opts.output = "UKBBoutput." + matlab.lang.makeValidName(tab.Pheno(1)) +...
            "." + string(datetime('today', 'Format', 'dd.MM.yyyy'));
else
    if ~isfield(opts, 'sheet'); opts.sheet = ""; end
    opts.output = opts.output + "." + opts.sheet + "." + matlab.lang.makeValidName(tab.Pheno(1));
    opts.output = regexprep(opts.output, '\.{2,}', '.');
end

for i = 1:numel(tab.SNP)
    fig = figure('visible','on');
    if opts.binary
        if opts.firth
            eff = tab.('OR(Firth')(i, :).';
        else
            eff = tab.OR(i, :).';
        end
        lineBegin = 1;
    else
        eff = tab.("β")(i, :).';
        lineBegin = 0;
    end
    ci = tab.("95% CI")(i, :).';
    missedData = ismissing(ci) | isnan(eff);
    ci = cellfun(@(x)sscanf(x, '[%f, %f]').', ci, 'uni', false);
    ci = vertcat(ci{:});
    eff(missedData) = [];
    
    % upper CI must be > eff, the otherwise may happen because of roundoff
    % when effects are close to zero (OR ~ 1). In this case upper CI is
    % estimate again assuming enough sample size (BETA +/- 1.96SE)
    f_roundoff = ci(:, 2) <= eff | ci(:, 1) >= eff | any(abs(ci - 1) < 1e-6, 2) | any(abs(ci - 0) < 1e-6, 2);
    if any(f_roundoff)
        se = tab.SE(i, :).';
        se(missedData) = [];
        if opts.binary
            ci(f_roundoff, :) = [exp(log(eff(f_roundoff)) - se(f_roundoff).*1.96), ...
                exp(log(eff(f_roundoff)) + se(f_roundoff).*1.96)];
        else
            ci(f_roundoff, :) = [eff(f_roundoff) - se(f_roundoff).*1.96, ...
                eff(f_roundoff) + se(f_roundoff).*1.96];
        end
    end
    
    if opts.strata
        Qlabels = tab.Strata(i, :).';
        N = "N = " + string(tab.N(i, :).');
        Qlabels(missedData) = [];
        N(missedData) = [];
        Ylabels = Qlabels;
        % Ylabels = join(pad([Qlabels, N], 'right'), '\newline'); % doesn't allign well
%         % latex format (there is almost no control on font name in MATLAB,
%         % so note used now).
%         Ylabels = compose('\\begin{tabular}{r} %s \\\\ %s \\end{tabular}', Qlabels(:), N(:));
%         Ylabels = replace(Ylabels, '<', '$<$');
%         Ylabels = replace(Ylabels, '>', '$>$');
%         Ylabels = replace(Ylabels, '≤', '$\leq$');
%         Ylabels = replace(Ylabels, '≥', '$\geq$');
        
        faceColor = 'r';
        markerSize = 8;
        vertH = 0.05;
    else
        Ylabels = "Overall";
        faceColor = 'b';
        markerSize = 10;
        vertH = 0.025;
    end
    Qtitle = "";
    YTEXT = tab.Pheno(1);

    % f = figure('visible','on');
    h = gca;
    h.Box = 'on';
    if opts.strata
        h.YLim = [0.5 numel(eff)+0.5];
    else
        h.YLim = [.75 1.25];
    end
    minLim = min(min(ci(:)), lineBegin);
    maxLim = max(max(ci(:)), lineBegin);
    h.XLim = [minLim maxLim]; 
    h.XLim(1) = h.XLim(1) - 0.12*abs(maxLim - minLim);
    h.XLim(2) = h.XLim(2) + 0.12*abs(maxLim - minLim);
    h.YTick = 1:numel(eff);
    h.YTickLabel = Ylabels;
    h.FontWeight = 'bold';
    h.LineWidth = 1.8;
    h.YAxis.TickLabelInterpreter = 'none';
    
%     % if latex is used
%     if contains(Ylabels, '\begin{tabular}')
%         h.YAxis.TickLabelInterpreter = 'latex';
%     end
    if opts.strata % add another axis for N
        sepdist = 0.1;
        h(2) = axes('YLim', h(1).YLim, 'FontWeight', h(1).FontWeight,...
            'LineWidth', h(1).LineWidth, 'XLim', h(1).XLim, 'Box', 'on',...
            'XTickLabel', [], 'Position', h(1).Position, ...
            'YTick', h(1).YTick, 'Color', 'none', 'TickLength', [0 0], ...
            'Parent', fig);
        h(1).YAxis.TickValues = h(1).YAxis(1).TickValues + sepdist;
        h(2).YAxis.TickValues = h(2).YAxis(1).TickValues - sepdist;
        h(2).YTickLabel = N;
    end
    % h.TickLength = [0 0];

    ylabel(h(1), Qtitle,'FontWeight','Bold','FontSize',11,'FontName','Arial');
    title(h(1), YTEXT, 'Interpreter', 'none')

    if ~opts.binary % Beta
        xlabel(h(1), '\fontsize{12} \bf \beta', 'Interpreter', 'tex',...
            'FontWeight','Bold','FontSize',11,'FontName','Arial');
        line(h(1), [lineBegin,lineBegin],h(1).YLim,'LineStyle','--','LineWidth',.7,'Color',[128 128 128]./256)
    else % OR
        xlabel(h(1), 'OR', 'FontWeight','Bold','FontSize',11,'FontName','Arial');
        line(h(1), [lineBegin,lineBegin],h(1).YLim,'LineStyle','--','LineWidth',.7,'Color',[128 128 128]./256)
    end
    hold(h, 'on')
    for ii = 1:numel(eff)
        line(h(1), ci(ii, :),[ii,ii],'LineStyle','-','LineWidth',1.3,'Color','k')
        line(h(1), [ci(ii, 1),ci(ii, 1)],[ii-vertH,ii+vertH],'LineStyle','-','LineWidth',1.3,'Color','k')
        line(h(1), [ci(ii, 2),ci(ii, 2)],[ii-vertH,ii+vertH],'LineStyle','-','LineWidth',1.3,'Color','k')
        plot(h(1), eff(ii),ii,'Marker','O','MarkerEdgeColor','k',...
            'MarkerFaceColor',faceColor,'MarkerSize',markerSize)
    end
    axis(h, 'square')
    h(1).Position(1) = 0.97 - h(1).Position(3); % move axes to the right border
    if numel(h) > 1
        h(2).Position = h(1).Position;
    end
    hold(h, 'off')

    tmpout = opts.output + "." + matlab.lang.makeValidName(tab.SNP(i));
    if opts.saveplot
        exportgraphics(fig, ...
            regexprep(tmpout + "." + opts.format, '\.{2,}', '.'),...
            'Resolution', opts.resolution)
    end
    
    if opts.savefig
        savefig(fig, regexprep(tmpout + ".fig", '\.{2,}', '.'))
    end
    close(gcf)
end
end

%% ------------------------------------------------------------------------
function [covar_matrix, opts] = getcovar(get_covars, eid, opts)
f_names = fieldnames(get_covars);
covar_matrix = nan(numel(eid), numel(f_names));

exeid = cell(numel(f_names), 1); % exclusion criteria for covariates
for i = 1:numel(f_names)
    if isfield(get_covars.(f_names{i}), 'exeid')
        exeid{i} = get_covars.(f_names{i}).exeid;
    end
    [f1, f2] = ismember(eid, get_covars.(f_names{i}).eid);
    f2(f2<1) = [];
    if get_covars.(f_names{i}).numericFlag
        covar_matrix(f1, i) = double(get_covars.(f_names{i}).rawUKB(f2));
    else
        if ~strcmpi(get_covars.(f_names{i}).tag, 'sex') % ignore categories
            covar_matrix(:, i) = zeros(numel(eid), 1);
            covar_matrix(f1, i) = 1; 
            continue
        end

        % sex
        u_char = unique(get_covars.(f_names{i}).rawUKB(f2));
        if numel(u_char) < 2
            error('covariates cannot contain < 2 categories!')
        end
        covar_temp = zeros(numel(f2), 1);
        get_covars.(f_names{i}).rawUKB = get_covars.(f_names{i}).rawUKB(f2);
        get_covars.(f_names{i}).eid = get_covars.(f_names{i}).eid(f2);
        if isfield(get_covars.(f_names{i}), 'termMeaning')
            get_covars.(f_names{i}).termMeaning = get_covars.(f_names{i}).termMeaning(f2);
        end

        for j = 1:numel(u_char)
            idx = get_covars.(f_names{i}).rawUKB == u_char(j);
            covar_temp(idx) = j-1;

            % for stratification and descriptive stats
            opts.("sex_"+(j-1)) = string(j-1);
            if isfield(get_covars.(f_names{i}), 'termMeaning')
                opts.("sex_"+(j-1)) = unique(get_covars.(f_names{i}).termMeaning(idx));
            end
        end
        covar_matrix(f1, i) = covar_temp;
    end
end

exeid(cellfun(@isempty, exeid)) = [];
if ~isempty(exeid)
    exeid = unique(vertcat(exeid{:}));
    exeidx = ismember(eid, exeid);
    if any(exeidx)
        covar_matrix(exeidx, :) = nan;
        if opts.verbose
            fprintf('%d samples will be remove due to exclusion criteria.\n', sum(exeidx))
        end
    end
end

if ~opts.defaultcovar % don't need to include additional covariates
    return
end
% add extra covariates
ext_covar = load(fullfile(opts.phenopath,'UKB_CONFOUNDERS_NUMERIC.mat'));
ext_covar = ext_covar.UKB_CONFOUNDERS;


[f1, f2] = ismember(eid, ext_covar.pcaeid); f2(f2<1) = [];
covar_matrix(f1, numel(f_names)+1:numel(f_names)+10) = ext_covar.pcadata(f2, :);

if ~opts.wes % ignore this for WES
    [f1, f2] = ismember(eid, ext_covar.binarrayeid); f2(f2<1) = [];
    covar_matrix(f1, size(covar_matrix, 2) + 1) = ext_covar.binarraydata(f2);
end

if ~isempty(exeid) % not needed though 
    covar_matrix(exeidx, :) = nan;
end

end

%% ------------------------------------------------------------------------
function tabALL = appendTab(tabALL, tabS, ct, k, opts)
fis = ["A2FREQ", "p", "p_adj", "BETA", "SE", "CI", "p_int", "p_int_adj",...
    "p_var", "p_var_adj", "BETA_int", "SE_int", "BETA_var", "SE_var", ...
    "N", "p_flic", "p_adj_flic", "BETA_flic", "p_flac", "p_adj_flac", ...
    "BETA_flac"];
for i = 1:numel(fis)
    if isfield(tabS, fis(i))
        tabALL.(fis(i))(ct, k) = tabS.(fis(i));
    end
end

if opts.binary
    tabALL.N_case(ct, k) = tabS.N_case;
    tabALL.AF_case(ct, k) = tabS.AF_case;
    tabALL.AF_control(ct, k) = tabS.AF_control;
else
    [tabALL.N_case(ct, k), tabALL.AF_case(ct, k), ...
        tabALL.AF_control(ct, k)] = deal(NaN);
end
end

%% ------------------------------------------------------------------------
function des = makedes(des, snp, dv, ct, k, opts)
% generates overall/stratified by genotype stat.
% snp = gather(snp);
dv = gather(dv);
dv(snp < 0) = []; snp(snp < 0) = [];

% for dosage inputs, we follow similar to what PLINK does:
% https://www.cog-genomics.org/plink/2.0/input#dosage_import_settings
% [0, 0.1] --> 0, [0.9, 1.1] --> 1, [1.9, 2] --> 2
cusnp = unique(snp(~isnan(snp)));
hardcalled= numel(cusnp) <= 3 && all(ismember(cusnp, [0, 1, 2]));
if ~hardcalled && ~ismissing(des.A1(ct)) && ~opts.customgeno  % not need to hard-call PRS or custom inputs
    f0 = snp <= 0.1 & snp >= 0;
    f1 = snp <= 1.1 & snp >= 0.9;
    f2 = snp <= 2 & snp >= 1.9;
    snp(f0) = 0; snp(f1) = 1; snp(f2) = 2;
    snp(~f0 & ~f1 & ~f2) = nan;
end

usnp = unique(snp(~isnan(snp)));
if numel(usnp) > 4
    usnp = usnp(1:4); % # for future: should use percentiles for PRS instead.
end

if opts.binary 
    if any(ismember(lower(des.Pheno), ["sex", "sex, n(%)", "sex, male, n(%)", "sex, female, n(%)"]))
        des.Pheno(ct, 1) = "Sex, " + lower(opts.sex_1) + ", n(%)";
    elseif ~endsWith(des.Pheno(ct, 1), ", n(%)")
        des.Pheno(ct, 1) = des.Pheno(ct, 1) + ", n(%)";
    end
    des.n_case(ct, k)    = sum(dv);
    des.n_control(ct, k) = sum(dv == 0);
    des.n_perc(ct, k)    = round(des.n_case(ct, k).*100./numel(dv), 3, 'significant');
     for i = 1:numel(usnp)
        fix = snp == usnp(i);
        des.usnp(ct, k, i)        = string(usnp(i));
        des.n_case_s(ct, k, i)    = sum(dv(fix));
        des.n_control_s(ct, k, i) = sum(dv(fix) == 0);
        des.n_perc_s(ct, k, i)    = round(des.n_case_s(ct, k, i) .*100./numel(dv(fix)), 3, 'significant');
     end
     
else
    
    des.n(ct, k)      = numel(dv);
    des.mean(ct, k)   = round(mean(dv), 3, 'significant');
    des.sd(ct, k)     = round(std(dv), 3, 'significant');
    des.median(ct, k) = round(median(dv), 3, 'significant');
    des.iqr(ct, k)    = round(iqr(dv), 3, 'significant');
    for i = 1:numel(usnp)
        fix = snp == usnp(i);
        des.usnp(ct, k, i)     = string(usnp(i));
        des.n_s(ct ,k, i)      = sum(fix);
        des.mean_s(ct, k, i)   = round(mean(dv(fix)), 3, 'significant');
        des.sd_s(ct, k, i)     = round(std(dv(fix)), 3, 'significant');
        des.median_s(ct, k, i) = round(median(dv(fix)), 3, 'significant');
        des.iqr_s(ct, k, i)    = round(iqr(dv(fix)), 3, 'significant');
    end
    
end

end

%% ------------------------------------------------------------------------
function [tab1, tab2] = toCleanTabDES(des, opts)
% overall table
vars = ["Pheno", "SNP", "CHR", "POS", "A1", "A2"];
if opts.binary
    vars = [vars, "n_case", "n_control", "n_perc"];
else
    vars = [vars, "n", "mean", "sd", "median", "iqr"];
end
if opts.strata
    vars = [vars, "Strata"];
end

des = struct2table(des);
overall_idx = ismember(des.Properties.VariableNames, vars);
tab1 = des(:, overall_idx);

if opts.binary
    tab1.n = tab1.n_case + "(" + tab1.n_perc + ")";
    [tab1.("mean±SD"), tab1.("median(IQR)")] = deal(strings(numel(tab1.n), 1));
    rem_cols = ismember(tab1.Properties.VariableNames, ["n_case", "n_control", "n_perc"]);
    tab1(:, rem_cols) = [];
else
    tab1.("mean±SD") = tab1.mean + "±" + tab1.sd;
    tab1.("median(IQR)") = tab1.median + "(" + tab1.iqr + ")";
    rem_cols = ismember(tab1.Properties.VariableNames, ["mean", "sd", "median", "iqr"]);
    tab1(:, rem_cols) = [];
end

if strcmp(opts.model, 'add')
    genoRAW = [tab1.A1+tab1.A1, tab1.A1+tab1.A2, tab1.A2+tab1.A2];
elseif strcmp(opts.model, 'dom')
    genoRAW = [tab1.A1+tab1.A1, tab1.A2 + "(dom)"];
else
    genoRAW = [tab1.A1+tab1.A1, tab1.A2 + "(rec)"];
end
tab1 = standardizeMissing(tab1, "");

% stratify by genotype
tab2 = ({});
des = splitvars(des(:, ~overall_idx));
usnp = find(startsWith(des.Properties.VariableNames, 'usnp_'));
% des = varfun(@squeeze, des(:, ~overall_idx));
% des.Properties.VariableNames = replace(des.Properties.VariableNames, 'squeeze_', '');

for j = 1:height(des)
    tmp = des(j, :);
    usnp_in = usnp(~ismissing(tmp{:, usnp}));
    geno = tmp{:, usnp_in};
    if ~ismissing(tab1.A1(j)) && ~ismissing(tab1.A2(j)) % only supports biallelic autosomal for now
        geno = genoRAW(j, tiedrank(double(geno)));
    end
    t = table; t.Pheno = tab1.Pheno(j);
    t.SNP = tab1.SNP(j); t.CHR = tab1.CHR(j); t.POS = tab1.POS(j);
    t.A1 = tab1.A1(j); t.A2 = tab1.A2(j);
    for i = 1:numel(usnp_in)
        
        if opts.binary
            t.("n|"+geno(i)) = tmp.("n_case_s_"+usnp_in(i)) + "(" + tmp.("n_perc_s_"+usnp_in(i)) + ")";
            [t.("mean±SD|"+geno(i)), t.("median(IQR)|"+geno(i))] = deal(strings(numel(t.("n|"+geno(i))), 1));
        else
            t.("n|"+geno(i)) = tmp.("n_s_"+usnp_in(i));
            t.("mean±SD|"+geno(i)) = tmp.("mean_s_"+usnp_in(i)) + "±" + tmp.("sd_s_"+usnp_in(i));
            t.("median(IQR)|"+geno(i)) = tmp.("median_s_"+usnp_in(i)) + "(" + tmp.("iqr_s_"+usnp_in(i)) + ")";
        end
    end
    tab2{j, 1} = t;
end

end

%%
function pruneMore(opts)

for i = 1:numel(opts.cleanOutputcheck) % loop over sheets
    res = string(readcell(opts.output+".xlsx", 'Sheet', ...
        replace(truncateStr(opts.cleanOutputcheck(i)), ["?", "\", "/"], "_"), 'TextType', 'string'));
    hdr = res(1, :); res(1, :) = [];
    emptyrows = all(ismissing(res), 2); % empty rows that must be removed
    if ~any(emptyrows) % no need to munge the table
        continue
    end
    extrahdrsIdx = circshift(emptyrows, 1); % other headers to be merged with hdr
    extrahdrs = res(extrahdrsIdx, :);
    res(extrahdrsIdx | emptyrows, :) = []; % remove extra hdrs and empty rows
    
    if any(hdr ~= extrahdrs, 'all') % merge different headers into one
        extrahdrs(all(hdr == extrahdrs, 2), :) = [];
        didx = find(any(hdr ~= extrahdrs, 1));
        for j = 1:numel(didx)
            thisel = [hdr(didx(j));extrahdrs(:, didx(j))];
            thisel(ismissing(thisel)) = [];
            hdr(didx(j)) = join(unique(thisel, 'stable'), '|');
        end
    end
    
    if any(arrayfun(@(x)numel(x{1}), hdr) > 15) % avoid long header names
        % remove everything after | (this happens only in descriptive
        % stat sheet with multiple variants)

        idx = hdr == "β|OR.Firth" | hdr == "OR.Firth|β"; % keep it as it is!
        hdr(~idx) = regexprep(hdr(~idx), '(?=\|).*', '');
        if any(contains(hdr, "median(IQR)")) % descriptive stat (genotypes)
            
            if opts.model == "dom"
                hdrtag = [".WT", ".DOM"];
            elseif opts.model == "rec"
                hdrtag = [".WT", ".REC"];
            else
                hdrtag = [".WT", ".HET", ".HOM"];
            end
            idx = reshape(find(ismember(hdr, ["n", "median(IQR)", "mean±SD"])), [], numel(hdrtag));
            for j = 1:size(idx, 2)
                hdr(idx(:, j)) = hdr(idx(:, j)) + hdrtag(j);
            end
        end
    end

    if numel(unique(hdr)) ~= numel(hdr)
        % this happens because some variants are not bi-allelic (or when
        % continuous exposures are used instead of variants!).
        hdr = matlab.lang.makeUniqueStrings(hdr);
    end
    res = array2table(res, 'VariableNames', hdr);
    
    % check if contains double data
    rescheck = varfun(@double, res);
    rescheckIdx = find(all(~isnan(rescheck.Variables) | ismissing(res)));
    for j = 1:numel(rescheckIdx)
        res.(rescheckIdx(j)) = rescheck.(rescheckIdx(j));
    end
    writetable(res, opts.output+".xlsx", 'AutoFitWidth', true, ...
        'Sheet', replace(truncateStr(opts.cleanOutputcheck(i)), ["?", "\", "/"], "_"), 'WriteMode', 'overwritesheet')
end

end

%%
function [geno, wesinfo, qc] = prunewescalls(geno, wesinfo, qc, opts)
% only for WES calls from getWES function, to remove missing calls, and
% keep qced eids. The output is then fed into R for SKAT/SKAT-O/and Burden
% tests. Approximate coefficients are calculated based on collapsing
% variable, and should be treated only for a rough estimation of direction
% of association and not the effect size itself (refer to SAIGE-GENE paper
% for more details)

% check if there are unannotated variants missing from wesinfo table
if any(~ismember(geno.bim.snp, wesinfo.id))
    unannotatedIDs = setdiff(geno.bim.snp, wesinfo.id);
    hdrs = wesinfo.Properties.VariableNames;
    wesinfoUnannotated = splitvars(table(unannotatedIDs, ...
        strings(numel(unannotatedIDs), numel(hdrs) - 1)));
    wesinfoUnannotated.Properties.VariableNames = hdrs;
    wesinfo = [wesinfo; wesinfoUnannotated];
end

% match wesinfo to geno.bim.snp
[~, idx] = ismember(geno.bim.snp, wesinfo.id); idx(idx<1) = [];
wesinfo = wesinfo(idx, :);

[f_eid, f_rare] = ismember(qc, geno.fam);
f_rare(f_rare < 1) = [];
qc = qc(f_eid);
geno.bed = geno.bed(f_rare, :);
geno.fam = geno.fam(f_rare);
geno.bed = geno.bed(:, 1:2:end) + geno.bed(:, 2:2:end);

f_missing = geno.bed < 0; % missing idx
missingness = (sum(f_missing, 1) ./ size(f_missing, 1).*100).'; 
maf = zeros(size(geno.bed, 2), 1);
if opts.verbose; fprintf('\n'); end
for i = 1:numel(missingness)
    if (sum(geno.bed(:, i) == 0) < sum(geno.bed(:, i) == 2))
        geno.bed(:, i) = 2 - geno.bed(:, i);
        % swap a1/a2
        newa1 = geno.bim.a2(i);
        geno.bim.a2(i) = geno.bim.a1(i);
        geno.bim.a1(i) = newa1;
    end
    maf(i, 1) = mean(geno.bed(~f_missing(:, i), i))./2;

    if opts.verbose
        if missingness(i) > 2 % missigness > 2% to be removed
            fprintf('ID=%s\tMAF=%.3g\tmissigness=%.3g%%\n',...
                geno.bim.snp(i), maf(i), missingness(i))
        else
            fprintf('ID=%s\tMAF=%g\n', geno.bim.snp(i), maf(i))
        end
    end
end
geno.bed(f_missing) = 9; % Set missing genotypes to 9 (as in SKAT)
if opts.singlerare
    missingness = missingness == 100;
else
    missingness = missingness > (opts.missingness.*100);
end

f_zero = maf == 0;
if any(f_zero) && opts.verbose
    fprintf('%d variants were excluded with MAF = 0\n', sum(f_zero))
end

if ~opts.singlerare
    f_thresh = maf >= opts.raremaf;
    if any(f_thresh) && opts.verbose
        fprintf('%d variants were excluded with MAF >= %.2g\n', sum(f_thresh), opts.raremaf)
    end
end

if any(missingness)
    if opts.verbose
        if opts.singlerare
            fprintf('%d variants were excluded with 100%% missigness\n', sum(missingness))
        else
            fprintf('%d variants were excluded with missigness >%.2f%%\n', sum(missingness), opts.missingness.*100)
        end
    end
end

if opts.singlerare % for single rare variant analysis
    f_thresh = f_zero | missingness;
else
    f_thresh = f_zero | f_thresh | missingness;
end

geno.bed(:, f_thresh) = []; wesinfo(f_thresh, :) = []; % maf(f_thresh) = [];

bimfields = fieldnames(geno.bim);
for i = 1:numel(bimfields)
    geno.bim.(bimfields{i})(f_thresh) = [];
end

% remove missing (9)--> SKAT impuptes missing genotypes, if you want that
% comment followings
if opts.singlerare % should be treated differently
    geno.bed(geno.bed == 9) = -2;
else
%     f_9 = any(geno.bed == 9, 2);
%     fprintf('%d samples with missing genotype were removed.\n', sum(f_9))
%     geno.bed(f_9, :) = []; geno.fam(f_9) = []; qc(f_9) = [];
end

f_maf_0 = all(~geno.bed, 1);
if any(f_maf_0)
    if opts.verbose
        fprintf('%d variants were further removed with MAF 0.\n', sum(f_maf_0))
    end
    geno.bed(:, f_maf_0) = []; wesinfo(f_maf_0, :) = [];
    
    bimfields = fieldnames(geno.bim);
    for i = 1:numel(bimfields)
        geno.bim.(bimfields{i})(f_maf_0) = [];
    end
end

if opts.verbose
    fprintf('%d variants and %d samples passed QC.\n\n', ...
        size(geno.bed, 2), size(geno.bed, 1))
end

if ~opts.singlerare
    maf = (0);
    if ~isempty(geno.bim.snp)
        for i = 1:size(geno.bed, 2)
            tmp = geno.bed(:, i); tmp(tmp == 9) = [];
            maf(i, 1) = mean(tmp)./2;
        end
        wesinfo.maf = maf;
    end
end
end

%%
function [skat, res] = callSKAT(geno, dv, covar, qc, opts)
% returns results from SKAT R package (skat table), and from MATLAB on
% collapsing variable for both additive and dominant models (res table);
[f_geno, f_qc] = ismember(geno.fam, qc); f_qc(f_qc<1) = [];
geno.bed(~f_geno, :) = []; geno.fam(~f_geno) = [];
dv = dv(f_qc); 
clear qc
if ~isempty(covar)
    covar = covar(f_qc, :);
end

fnan = any(isnan(covar), 2);
dvClean = dv(~fnan);
if ~strcmp(opts.transform, "none") && ~opts.skiptransform
    dvClean = feval(opts.transform, gather(dvClean));
    dv = feval(opts.transform, gather(dv));
    if opts.gpu
        dvClean = gpuArray(dvClean);
        dv = gpuArray(dv);
    end
end

% R = corr(double(geno.bed), gather(dv));
% 
% if opts.debug == -1
%     geno.bed(:, R < 0) = [];
% elseif opts.debug == 1
%     geno.bed(:, R > 0) = [];
% end

% ax = drawGene(opts.wesinfo.gene(1), geno.bim, R);

% calculate freq. for this trait for both num markers (MATLAB) and decide
% on rare+common or rare variant analysis in SKAT
if any(geno.bed == 9, 'all') % to be imputed in SKAT
    checkmaf = (0);
    for i = 1:size(geno.bed, 2)
        tmp = geno.bed(:, i); tmp(tmp == 9) = [];
        checkmaf(i, 1) = mean(tmp)./2;
    end
    
else
    checkmaf = mean(geno.bed)./2;
end

% get effect size from collapsing score
score = sum(geno.bed > 0 & geno.bed ~= 9, 2); % additive score
if opts.gpu; score = gpuArray(score); end
score_dom = score; score_dom(score_dom > 1) = 1; % dominant

mdl = ({});
% A 2X2 cell:   unadj  adj
%           add (1,1) (1,2)
%           dom (2,1) (2,2)
if opts.binary 
    cols = ["P", "P-adj", "OR", "95% CI"];
    dblvars = ["P", "P-adj", "OR", "AF"];
    if ~opts.desonly
        mdl{1, 1} = compact(fitglm(score, dv, 'Distribution', 'binomial', 'Link', 'logit'));
        if isempty(covar)
            mdl{1, 2} = mdl{1, 1};
        else
            mdl{1, 2} = compact(fitglm([score(~fnan), covar(~fnan, :)], dvClean, 'Distribution', 'binomial', 'Link', 'logit'));
        end
    end
else
    cols = ["P", "P-adj", "β", "95% CI"];
    dblvars = ["P", "P-adj", "β", "AF"];
    if ~opts.desonly
        mdl{1, 1} = compact(fitlm(score, dv));
        if isempty(covar)
            mdl{1, 2} = mdl{1, 1};
        else
            mdl{1, 2} = compact(fitlm([score(~fnan), covar(~fnan, :)], dvClean));
        end
    end
end

% for dominant score
if any(score > 1) && ~opts.desonly
    if opts.binary
        mdl{2, 1} = compact(fitglm(score_dom, dv, 'Distribution', 'binomial', 'Link', 'logit'));
        if isempty(covar)
            mdl{2, 2} = mdl{2, 1};
        else
            mdl{2, 2} = compact(fitglm([score_dom(~fnan), covar(~fnan, :)], dvClean, 'Distribution', 'binomial', 'Link', 'logit'));
        end
    else
        mdl{2, 1} = compact(fitlm(score_dom, dv));
        if isempty(covar)
            mdl{2, 2} = mdl{2, 1};
        else
            mdl{2, 2} = compact(fitlm([score_dom(~fnan), covar(~fnan, :)], dvClean));
        end
    end
elseif ~opts.desonly
    mdl{2, 1} = mdl{1, 1};
    mdl{2, 2} = mdl{1, 2};
end

% remove nan values in covar matrix
dv = dvClean; clear dvClean
if ~isempty(covar)
    covar(fnan, :) = [];
    score_dom(fnan) = []; score(fnan) = [];
    geno.bed(fnan, :) = []; geno.fam(fnan) = [];
end

% format as a table
cols = [cols, "AF", "N", "N.case", "MAC.case", "MAC.control"];
res = cell(2, numel(cols));
for i = 1:2
    if ~opts.desonly
        res{i, 1} = mdl{i, 1}.Coefficients.pValue(2);
        res{i, 2} = mdl{i, 2}.Coefficients.pValue(2);
        coefCI = mdl{i, 2}.coefCI; coefCI = gather(coefCI(2, :));
        beta = mdl{i, 2}.Coefficients.Estimate(2);
    end
    
    if i == 1
        res{i, 5} = mean(score)./2;
    else
        res{i, 5} = mean(score_dom);
    end
    res{i, 6} = numel(dv);
    
    if opts.binary && ~opts.desonly
        beta = exp(beta); coefCI = exp(coefCI);
        res{i, 7} = sum(dv);
        if i == 1
            res{i, 8} = sum(score(dv == 1));
            res{i, 9} = sum(score(dv == 0));
        else
            res{i, 8} = sum(score_dom(dv == 1));
            res{i, 9} = sum(score_dom(dv == 0));
        end
    else
        res(i, 7:9) = {nan};
    end
    
    if ~opts.desonly
        res{i, 3} = beta;
        if opts.binary
            roundType = 'decimals';
        else
            roundType = 'significant';
        end
        res{i, 4} =  "[" + round(coefCI(1), opts.ciround, roundType) + ", " ...
            + round(coefCI(2), opts.ciround, roundType) + "]";
    end
end

res = cell2table(cellfun(@gather, res, 'uni', false), 'VariableNames', cols);
% round to 3 significant digits
if ~opts.desonly
    tab_double = varfun(@(x)round(x, opts.betaround, 'significant'), res,...
        'InputVariables', dblvars);
    for i = 1:numel(dblvars)
        res.(dblvars{i}) = tab_double.(i);
    end
end
res.Pheno = repmat(opts.westrait, 2, 1);
res.Model = ["Additive"; "Dominant"];
res = movevars(res, {'Pheno', 'Model'}, 'Before', 'P');
res.("N.marker") = sum(checkmaf > 0).*ones(size(res, 1), 1); % those wit a MAF ~ 0
clear tabALL_double


if opts.desonly
    res(:, 2:6) = [];
    skat = table;
    return
end
% prepare data for running in SKAT
if ~isempty(covar); covar = gather(covar); end
dv = gather(dv);

% modelNull = fitlm(covar, dv);
% v = optimvar('v',numel(dv),1,'Type','integer','LowerBound',0,'UpperBound',1);
% prb = optimproblem;
% G = v-covar*(((covar'*covar)\covar')*v);
% t = ((G)'*modelNull.Residuals.Raw);
% prb.Objective = t;
% optsLin = optimoptions('intlinprog','Display','final','MaxTime',160);
% prb.ObjectiveSense = 'maximize';
% sol = solve(prb,'options',optsLin);
% R = corr(double(geno.bed), sol.v);
% R = corr(double(geno.bed), dv);
% f = R < 0;
% 
% if opts.debug == -1
%     geno.bed(:, R < 0) = [];
% elseif opts.debug == 1
%     geno.bed(:, R > 0) = [];
% end

% check if rare+common variant analysis shoud be performed instead
if any(checkmaf >= 0.01)
    opts.rarecommon = true;
else
    opts.rarecommon = false;
end

% generate unique random string for parallel calls to gwasrunner
charSet = ['a':'z',upper('a':'z'),'0':'9'];
randStr = string(charSet(randi(numel(charSet),1, 10)));
maxloop = 1;
while exist("geno.skat."+randStr+".mat", 'file')
    randStr = string(charSet(randi(numel(charSet),1, 10)));
    maxloop = maxloop + 1;
    if maxloop > 5
        error('cannot generate unique random string for SKAT package :(')
    end
end

geno = geno.bed;
save("geno.skat."+randStr+".mat", 'geno')
if ~isempty(covar); save("covar.skat."+randStr+".mat", 'covar'); end
save("dv.skat."+randStr+".mat", 'dv')
clear dv geno

if opts.binary
    type = "D";
else
    type = "C";
end

currDir = strrep(pwd, '\', '/');
Rcode = string;
Rcode(1, 1) = "setwd('" + currDir + "')";
Rcode(2, 1) = "library(SKAT)";
Rcode(3, 1) = "geno <- R.matlab::readMat('geno.skat."+randStr+".mat')$geno";
Rcode(4, 1) = "dv <- R.matlab::readMat('dv.skat."+randStr+".mat')$dv";

if isempty(covar)
    Rcode(6, 1) = 'obj<-SKAT_Null_Model(dv ~ 1, out_type="' + type + '")';
else
    Rcode(5, 1) = "covar <- R.matlab::readMat('covar.skat."+randStr+".mat')$covar";
    Rcode(6, 1) = 'obj<-SKAT_Null_Model(dv ~ covar, out_type="' + type + '")';
end
Rcode(7, 1) = "rm(covar); rm(dv)";

if opts.rarecommon
    if opts.binary
        Rcode(8, 1) = 'out <- SKAT_CommonRare_Robust(geno, obj, method="SKAT", CommonRare_Cutoff = 0.01)';
        Rcode(9, 1) = 'outo <- SKAT_CommonRare_Robust(geno, obj, method="SKATO", CommonRare_Cutoff = 0.01)';
        Rcode(10, 1) = 'outb <- SKAT_CommonRare_Robust(geno, obj, method="Burden", CommonRare_Cutoff = 0.01)';
    else
        % method = "C" represents the combined sum test (default), and "A" represents the adaptive sum test
        Rcode(8, 1) = 'out <- SKAT_CommonRare(geno, obj, CommonRare_Cutoff = 0.01)';
        Rcode(9, 1) = 'outb <- SKAT_CommonRare(geno, obj, r.corr.rare=1, r.corr.common=1, CommonRare_Cutoff = 0.01)';
    end
else
    if opts.binary
        Rcode(8, 1) = 'out <-SKATBinary_Robust(geno, obj, kernel = "linear.weighted", weights.beta=c(1,25))';
        % SKAT-O: Combined Test of burden test and SKAT and a better type I error control
        Rcode(9, 1) = 'outo<-SKATBinary_Robust(geno, obj, method="SKATO")';
        Rcode(10, 1) = 'outb<-SKATBinary_Robust(geno, obj, method = "Burden")';
    else
        Rcode(8, 1) = 'out <-SKAT(geno, obj, kernel = "linear.weighted", weights.beta=c(1,25))';
        % SKAT-O: Combined Test of burden test and SKAT and a better type I error control
        Rcode(9, 1) = 'outo<-SKAT(geno, obj, method="SKATO")';
        Rcode(10, 1) = 'outb<-SKAT(geno, obj, r.corr=1)';
    end
end

if opts.rarecommon && ~opts.binary
    Rcode(11, 1) = 'res <- data.frame(Method = c("SKAT (1,25)", "Burden", "SKAT-O"), P = c(out$p.value, outb$p.value, NA), Q = c(out$Q, outb$Q, NA), N.marker = c(out$param$n.marker.test, outb$param$n.marker.test, NA))';
else
    Rcode(11, 1) = 'res <- data.frame(Method = c("SKAT (1,25)", "Burden", "SKAT-O"), P = c(out$p.value, outb$p.value, outo$p.value), Q = c(out$Q, outb$Q, NA), N.marker = c(out$param$n.marker.test, outb$param$n.marker.test, outo$param$n.marker.test))';
end
Rcode(12, 1) = "write.table(res, 'gwasrunner_skat_"+randStr+".txt', quote = F, col.names = T, row.names = F, sep = '\t')";
writematrix(Rcode, "gwasrunner_skat_"+randStr+".r", 'QuoteStrings', false, 'FileType','text')

if opts.verbose
    if opts.rarecommon
        fprintf('\tRunning SKAT and Burden tests (common + rare variant test)...')
    else
        fprintf('\tRunning SKAT, SKAT-O and Burden tests...')
    end
end

[timeM2R, failedM2R] = MATLAB2Rconnector("gwasrunner_skat_"+randStr+".r",...
    'delr', true, 'maxTime', opts.maxTime);

if ~exist("gwasrunner_skat_"+randStr+".txt", 'file') && opts.error
%     delete("geno.skat."+randStr+".mat")
%     delete("covar.skat."+randStr+".mat")
%     delete("dv.skat."+randStr+".mat")
    error('gwasrunner:callSKAT', 'something went wrong with %s!\n', ("gwasrunner_skat_"+randStr+".r"))
end

if ~failedM2R
    if opts.verbose; fprintf('\b\bDone.\n'); end
end

if opts.verbose
    fprintf('\tElapsed time in SKAT is %.2f seconds.\n', toc + timeM2R)
end

delete("geno.skat."+randStr+".mat")
delete("covar.skat."+randStr+".mat")
delete("dv.skat."+randStr+".mat")

if ~opts.error && ~exist("gwasrunner_skat_"+randStr+".txt", 'file')
    % don't interrupt gwasrunner flow if SKAT fails
    skat = table(["SKAT (1,25)"; "Burden"; "SKAT-O"], nan(3, 1),...
        nan(3, 1), nan(3, 1), ...
        'VariableNames', {'Method', 'P', 'Q', 'N.marker'});
else
    skat = readtable("gwasrunner_skat_"+randStr+".txt", 'ReadVariableNames', true,...
        'NumHeaderLines', 0, 'Delimiter', '\t', 'TextType', 'string',...
        'VariableNamingRule', 'preserve');
    delete("gwasrunner_skat_"+randStr+".txt");
    skat.P = round(skat.P, opts.betaround, 'significant');
end

skat.Pheno = repmat(opts.westrait, size(skat, 1), 1);
skat = movevars(skat, 'Pheno', 'Before', 'Method');

if opts.rarecommon
    skat.Method = ["SKAT (common&rare|cutoff=0.01)"; ...
        "Burden (common&rare|cutoff=0.01)"; ...
        "SKAT-O (common&rare|cutoff=0.01)"];
end
end

%%
function writeWESresults(skat, collapseres, opts)

%@22NOV2024: keep only subset of variants tested if method is REGENIE
if opts.regenie
    snplist = opts.regenied.snplistAll;
    snplist(cellfun(@isempty, snplist)) = [];
    snplist = vertcat(snplist{:});
    snplist.Pheno(snplist.Pheno == "") = missing;
    snplist = fillmissing(snplist, "previous", DataVariables="Pheno");
    snplist.mask = extractAfter(snplist.mask, "mask_");

    snplist_vars = arrayfun(@(x)x.split(","), snplist.variants, uni=false);
    snplistU = unique(vertcat(snplist_vars{:}));
    assert(all(ismember(snplistU, opts.wesinfo.id)), "some variants tested in REGENIE cannot be found in 'wesinfo'!")
    opts.wesinfo(~ismember(opts.wesinfo.id, snplistU), :) = [];
    wes_info = opts.wesinfo;

    % add extra columns to wesinfo to include snplist (variant.list is
    % deprecated now because of limit on maximum character allowed to store
    % in a cell in Excel)
    snplist.idx = cellfun(@(x) ismember(opts.wesinfo.id, x), snplist_vars, uni=false);
    upheno = unique(snplist.Pheno);
    for k = 1:numel(upheno)
        tmp = snplist(snplist.Pheno == upheno(k), :);
        mask_joined = strings(height(wes_info), height(tmp));
        for j = 1:height(tmp)
            mask_joined(tmp.idx{j}, j) = tmp.mask(j);
        end
        mask_joined = join(mask_joined, ",");
        mask_joined = regexprep(mask_joined, '(\,+)', ',');
        wes_info.(upheno(k)) = regexprep(mask_joined, '(^,$|,$|^,)', '');

    end

    opts.wesinfo = wes_info;
end

writetable(opts.wesinfo, opts.output+".xlsx", 'AutoFitWidth', true, ...
    'Sheet', opts.sheet + ".Info", 'WriteMode', 'overwritesheet')

skat(cellfun(@isempty, skat)) = [];
collapseres(cellfun(@isempty, collapseres)) = [];
if ~isempty(skat), skat = vertcat(skat{:}); end

if opts.regenie
    
    if isfield(opts.regenied.geno, 'bfile') % delete merged geno file
        if ~isfield(opts.regenied.geno, "anno_file")
            rmfiles = opts.regenied.geno.bfile;
            rmfiles = regexprep(rmfiles, ".bed$", "");
            rmfiles = rmfiles + [".bed", ".bim", ".fam",...
                ".pgen", ".pvar", ".psam", ".bgen", ".bgen.bgi", ".sample"];    
        else
            rmfiles = [opts.regenied.geno.anno_file, opts.regenied.geno.mask_def, opts.regenied.geno.set_list];
        end

        rmfiles(~isfile(rmfiles)) = [];
        if ~opts.debug %@20MARCH2025: keep junk files
            delete(rmfiles{:})
        end
    end
    
    if ~isempty(skat)
        [skat, leg] = pruneGeneBasedResults(skat);
        for i = 1:skat.dict.numEntries
            writetable(skat.("x"+i), opts.output + ".xlsx", ...
                Sheet=skat.dict("x"+i), AutoFitWidth=true, ...
                WriteMode="overwritesheet");
        end

        % write legend
        legkeys = leg.dict.keys;
        for i = 1:leg.dict.numEntries
            try
                writetable(leg.(legkeys(i)), opts.output + ".xlsx", ...
                    Sheet="Legend", AutoFitWidth=true, Range=leg.dict(legkeys(i)));
            catch
                writematrix(leg.(legkeys(i)), opts.output + ".xlsx", ...
                    Sheet="Legend", Range=leg.dict(legkeys(i)));
            end
        end
    end
    
    %@23NOV2024: deprecated, see above
    % snplist = opts.regenied.snplistAll;
    % snplist(cellfun(@isempty, snplist)) = [];
    % snplist = vertcat(snplist{:});
    % snplist = fillmissing(snplist, "previous", DataVariables="Pheno");
    % if ~isempty(snplist)
    %     writetable(snplist, opts.output + ".xlsx", ...
    %         'AutoFitWidth', true, 'Sheet', opts.sheet + ".variant.list", ...
    %         'WriteMode', 'overwritesheet')
    % end
    return
end

% from SKAT package
if ~isempty(skat)
    % reformat skat table for saving
    emptyrow = skat(1, :); emptyrow{:, :} = NaN;
    skat_out = cell({});
    for i = 1:size(skat, 1)/3
        skat_out{i, 1} = [skat(3*i-2:3*i, :); emptyrow];
    end
    skat_out = vertcat(skat_out{:});

    skat_out.Pheno(setdiff(1:size(skat_out,1),1:4:end)) = missing;
    writetable(skat_out, opts.output+".xlsx", 'AutoFitWidth', true, ...
        'Sheet', 'SKAT', 'WriteMode', 'overwritesheet')
end

if ~isempty(collapseres)
    % write collapsing variable results
    binaryflag = false(numel(collapseres), 1);
    for i = 1:numel(collapseres)
        if any(strcmp(collapseres{i}.Properties.VariableNames, 'OR'))
            binaryflag(i, 1) = true;
        end
    end

    if all(binaryflag) || all(~binaryflag)
        collapseres = vertcat(collapseres{:});
    else
        for i = 1:numel(collapseres)
            collapseres{i}.Properties.VariableNames{5} = 'β|OR';
        end
        collapseres = vertcat(collapseres{:});
    end

    emptyrow = collapseres(1, :); emptyrow{:, :} = NaN;
    collapseres_out = cell({});

    for i = 1:size(collapseres, 1)/2
        collapseres_out{i, 1} = [collapseres(2*i-1:2*i, :); emptyrow];
    end
    collapseres_out = vertcat(collapseres_out{:});
    collapseres_out.Pheno(setdiff(1:size(collapseres_out,1),1:3:end)) = missing;
    writetable(collapseres_out, opts.output+".xlsx", 'AutoFitWidth', true, ...
        'Sheet', 'Collapsed', 'WriteMode', 'overwritesheet')
end

end

%% ------------------------------------------------------------------------
function [res, leg] = pruneGeneBasedResults(tab)

% prunes gene-based tests from REGENIE 
% check ADD/DOM/REC test
adc = ["ADD", "DOM", "REC"];
adcIdx = arrayfun(@(x)contains(tab.A1, x), "." + adc, uni=false);
adc = adc(cellfun(@any, adcIdx));

tab.A1(tab.A1 == "Joint.GENE_P") = "Joint." + adc + "-GENE_P";
tests = extract(tab.A1, "." + adc + wildcardPattern + textBoundary("end"));
tab.Mask = strrep(tab.A1, tests, "");
pheno = unique(tab.Pheno);
mask = unique(tab.Mask);
% test = adc + ["-BURDEN-MINP", "-BURDEN-SBAT", "-BURDEN-ACAT", "",...
%     "-ACATO", "-ACATV", "-SKAT", "-SKATO",...
%     "-SKATO-ACAT", "-ACATO-FULL"];% for reformatting the table
test = unique(erase(tests, textBoundary("start") + "."))';

% maf groups (+singleton and Joint)
% maf = "0.01"; % maf cutoff to filter
maf = unique(extractAfter(mask, "."));
maf(ismissing(maf)) = [];
joint_mask = find(mask.lower == "joint", 1);
if ~isempty(joint_mask), maf = union(maf, mask(joint_mask)); end

mask(~endsWith(mask, maf)) = [];
selcols = string(union(setdiff(tab.Properties.VariableNames, ["A1", "P"], "stable"), "GO", "stable"));
sortedCols = ["Pheno", "ID", "CHR", "Mask", "BETA", "SE", test, "N"];
for j = 1:numel(maf)
    out = ({}); ct = 1;
    tmpout = tab;  
    tmpout(~endsWith(tmpout.Mask, maf(j)), :) = [];

    for i = 1:numel(pheno)
        tmp1 = tmpout;
        tmp1(tmp1.Pheno ~= pheno(i), :) = [];
        if isempty(tmp1), continue; end
        for k = 1:numel(mask)
            tmp = tmp1;
            tmp(tmp.Mask ~= mask(k), :) = [];
            if isempty(tmp), continue; end
            
            uID = unique(tmp.ID); 
            for l = 1:numel(uID) % loop over genes
                tmp2 = tmp(tmp.ID == uID(l), :);
                tmp2 = sortrows(tmp2, "SE");
                A1 = strrep(tmp2.A1, tmp2.Mask + ".", "");
                P = tmp2.P;
                tmp3 = nan(1, numel(test));
                [f1, f2] = ismember(test, A1);
                tmp3(1, f1) = P(f2(f1));
                
                cols = string(tmp2.Properties.VariableNames);
                if ~any(cols.lower.startsWith("go"))
                    selcols(selcols == "GO") = [];
                end
                tmp2 = tmp2(1, selcols);
                tmp2{:, test} = tmp3;
                out{ct} = tmp2; ct = ct + 1;
            end
        end 
    end
    
    out(cellfun(@isempty, out)) = [];
    out = vertcat(out{:});
    out = rmmissing(out, 2, "MinNumMissing", height(out));
    cols = string(out.Properties.VariableNames);
    beta_idx = contains(cols.lower, ["or|", "β|"]) | ismember(cols.lower, ["β", "or"]);
    scols = sortedCols;
    if any(beta_idx)
        scols(scols == "BETA") = cols(beta_idx);
    else
        scols(scols == "BETA") = [];
    end
    [f1,f2] = ismember(scols, cols);
    f3 = ~ismember(cols, scols);
    scols = [cols(f2(f1)), setdiff(cols(f3), scols, "stable")];

    res.("x"+j) = out(:, scols);
end

if iscolumn(maf), maf = maf'; end
res.dict = dictionary("x"+(1:numel(maf)), maf);

% add legend
% find different masks and tests
ucats = unique(tab.A1);
lastDot = regexp(ucats, '[.]','end');
leg.test = strings(numel(ucats), 1);
leg.mask = extractBefore(ucats, ".");
for i = 1:numel(lastDot)
    pos = lastDot{i}(end);
    tmp = ucats(i).char;
    leg.test(i) = string(tmp(pos+ 1:end));
end
leg.maf = replace(ucats, leg.mask + ".", "");
leg.maf = extractBefore(leg.maf, "." + leg.test);

fis = string(fieldnames(leg));
for i = 1:numel(fis)
    leg.(fis(i))(ismissing(leg.(fis(i))) | leg.(fis(i)) == "") = [];
    leg.(fis(i)) = unique(leg.(fis(i)));
end

leg.mask = erase(leg.mask, "mask_");
leg.mask(leg.mask.lower == "joint") = [];
tab = array2table(strings(numel(leg.mask), 5), "VariableNames", ...
    ["cadd", "revel", "polyphen", "sift", "clinvar"]);
cols = string(tab.Properties.VariableNames);
for i = 1:numel(cols)
    switch cols(i)
        case "cadd", tag = "[c](\d+)";
        case "revel", tag = "[r](\d+)";
        case "polyphen", tag = "[^\w]p(?<n>\w{1})([^\w]|$)";
        case "sift", tag = "[^\w]s(?<n>\w{1})([^\w]|$)";
        case "clinvar", tag = "(cvar)";
    end
    
    if any(ismember(cols(i), ["polyphen", "sift"]))
        idx = regexp(leg.mask, tag, 'names');
    else
        idx = regexp(leg.mask, tag, 'tokens');
    end

    if isempty(idx), continue; end % a custom mask
    rmidx = cellfun(@isempty, idx);
    idx(rmidx) = []; idx = vertcat(idx{:});
    if isstruct(idx), idx = struct2array(idx); end
    tab.(cols(i))(~rmidx) = string(idx);
    if cols(i) == "revel"
        tab.(cols(i)) = string(double(tab.(cols(i)))/10);
    elseif any(ismember(cols(i), ["polyphen", "sift"]))
        tab.(cols(i)) = replace(tab.(cols(i)), ["r", "s"], ["Relaxed", "Strict"]);
    elseif cols(i) == "clinvar"
        tab.(cols(i)) = replace(tab.(cols(i)), "cvar", "ClinVar");
    end
    
    if any(ismember(cols(i), ["cadd", "revel"]))
        tab.(cols(i))(~rmidx) = "≥" + tab.(cols(i))(~rmidx);
    elseif any(ismember(cols(i), ["polyphen", "sift"]))
        tab.(cols(i))(~rmidx) = "=" + tab.(cols(i))(~rmidx);
    end
end

tab = standardizeMissing(tab, "");
tab.clinvar(~ismissing(tab.clinvar)) = "";
cols = cols.upper;
cols = replace(cols, ["POLYPHEN", "CLINVAR"], ["PolyPhen", "ClinVar"]);
for i = 1:height(tab)
    rule = "";
    if contains(leg.mask(i).lower, "lof")
        rule = "LoF or ";
    end

    if contains(leg.mask(i).lower, "missense")
        rule = rule + "missnese or ";
    end
    
    if contains(leg.mask(i), "&")
        glue = " and ";
    else
        glue = " or ";
    end
    
    idx = ~ismissing(tab{i, :});
    if all(~idx)
        tab.term(i) = "UNKNOWN";
        continue
    end
    glued = join(cols(idx) + tab{i, idx}, glue);
    glued(ismissing(glued)) = [];

    if ~ismissing(glued)
        if any(contains(glued, glue))
            rule = rule + "(" + glued + ")";
        else
            rule = rule + glued;
        end
    else
        rule = regexprep(rule, " or $", "");
    end

    tab.term(i) = rule;
end
leg.mask = table(leg.mask, tab.term, 'VariableNames', {'Mask', 'Definition'});

leg.test = table(["ADD", "SKAT", "SKATO", "SKATO-ACAT", "ACATV", "ACATO", ...
    "ACATO-FULL", "ACAT", "SBAT", "GENE_P"]', ...
    ["Burden test"
    "Variance component test"
    "Omnibus test combining features of SKAT and Burden"
    "Same as SKATO but using Cauchy combination method to maximize power across SKATO models"
    "Test using Cauchy combination method to combine single-variant p-values"
    "Omnibus test combining features of ACATV, SKAT and Burden"
    "Same as ACATO but using the larger set of SKATO models used in the SKATO test"
    "Joint test: combines the p-values of the individual burden masks using the Cauchy combination method"
    "Joint test: Sparse Burden Association Test combines these burden masks in a joint model imposing constraints of same direction of effects"
    "uses ACAT to combine the p-values of the SKATO, ACATV, Burden and SBAT tests"], ...
    'VariableNames', {'Test', 'Desciption'});

sheetrange1 = height(leg.test) + 3;
leg.mask.Definition(leg.mask.Definition == "") = "NA";
sheetrange2 = sheetrange1 + height(leg.mask) + 2;
leg.maf = ["Alternate AF bins: ", join(leg.maf, ",")];

leg.dict = dictionary(["test", "mask", "maf"], "A" + [1, sheetrange1, sheetrange2]);

end

%% ------------------------------------------------------------------------
function mustBePheno(inpheno, pth)
% validation function for input traits
pheno = {dir(fullfile(pth, '*.mat')).name}.';
pheno = regexprep(pheno, '.mat$', '');
if ~all(ismember(inpheno, pheno))
    eid = 'Pheno:notFound';
    msg = 'All input traits must be present in phenoPath!';
    throwAsCaller(MException(eid,msg))
end

end

%% ------------------------------------------------------------------------
function sheetName = trimSheetName(sheetName, N)
% trims if sheet names contain > 31 characters

if nargin < 2, N = 31; end

sheetName = string(sheetName);
if sheetName.strlength > N
    sheetName = char(sheetName);
    sheetName = string(sheetName(1:N));
end

end

%%
function printLogo
logo = string;
logo(1, 1) = "---------------------------------------------------";
% logo(2, 1) = "    ⠰⠤⣤⣀⠀⠀⠀⠀⢀⣠⠤⠴⠦⢤⣀⠀⠀⠀⠀⠀⣀⡤⠴⠦⢤⣄⠀⠀⠀⠀⠀⣀⡤⠴⠄";
% logo(3, 1) = "    ⢸⠀⠀⡏⢳⣄⢀⡴⢫⡅⠀⢸⠀⠀⢨⠳⣄⠀⡴⠋⡇⠀⢸⡇⠀⢸⠓⣦⠀⣠⠞⢡⠀⢸⡇";
% logo(4, 1) = "    ⢸⠀⠀⡇⠀⢙⣾⠀⢸⡇⠀⢸⠀⠀⢸⠀ ⢘⣧⡁⠀⡇⠀⢸⡇⠀⢸⠀⠈⢿⡃⠀⢸⠀⢸⡇";
% logo(5, 1) = "    ⢸⠀⠀⢇⣰⠏⠘⢦⣸⠇⠀⢸⠀⠀⠸⣠⠞⠈⠳⣀⠇⠀⢸⡇⠀⠸⣠⠞⠀⠹⣄⠸⠀⢸⡇";
% logo(6, 1) = "    ⠨⠤⠖⠋⠁⠀⠀⠀⠉⠓⠦⠬⠤⠶⠋⠁⠀⠀⠀⠈⠓⠦⠬⠤⠴⠋⠁⠀⠀⠀⠈⠛⠲⠬⠅";
logo(7, 1) = "   _____      ___   ___                            ";
logo(8, 1) = "  / __\ \    / /_\ / __|_ _ _  _ _ _  _ _  ___ _ _ ";
logo(9, 1) = " | (_ |\ \/\/ / _ \\__ \ '_| || | ' \| ' \/ -_) '_|";
logo(10, 1) = "  \___| \_/\_/_/ \_\___/_|  \_,_|_||_|_||_\___|_|  ";
logo(11, 1) = "Version 1.2.1";
logo(12, 1) = "University of Gothenburg. Sahlgrenska Academy";
logo(13, 1) = string(datetime('now'));
logo(14, 1) = "---------------------------------------------------";
logo(ismissing(logo)) = [];
fprintf('%s\n', logo)
end

%% ------------------------------------------------------------------------
function pred_covars = readPredFiles(predfiles, chrom)
%@17NOV2024: loads columns from predfiles corresponding to chromosomes in
%'chrom'

chrom = natsort(string(chrom));
if ~iscolumn(chrom), chrom = chrom'; end
chrom_num = double(erase(chrom, textBoundary("start") + "chr"));
pred_covars = cell(numel(predfiles), 1);
for k = 1:numel(predfiles)
    df = readtable(predfiles(k), FileType="text", TextType="string",...
        NumHeaderLines=0, ReadVariableNames=false);
    df.Properties.VariableNames = ["pheno", "path"];
    if ispc
        df.path = makeWSLpath(df.path, true);
    end

    opts = detectImportOptions(df.path, FileType="text", Delimiter=" ", NumHeaderLines=1);

    % set DataLines
    breaks = [1; diff(chrom_num) > 1];
    group_start = find(breaks);
    group_end = [group_start(2:end) - 1; length(chrom_num)];
    opts.DataLines = [chrom_num(group_start), chrom_num(group_end)] + 1; % first header line

    data = readmatrix(df.path, opts, OutputType="double")';

    
    % read eids
    % opts = detectImportOptions(df.path, FileType="text", Delimiter=" ", TextType="string");
    opts.DataLines = [1, 1];
    eid = readmatrix(df.path, opts, OutputType="string").';
    eid = double(eid.extractBefore("_"));
    eid(1) = [];

    hdr = "chr" + data(1, :); data(1, :) = [];
    data = array2table([eid, data], VariableNames=["eid", string(hdr)]);
    data = data(:, ["eid";chrom]);
    clear eid

    pred_covars{k, 1} = df.pheno;
    pred_covars{k, 2} = data;
    clear data

end

end