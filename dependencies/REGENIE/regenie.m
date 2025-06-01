function cmd = regenie(p)
% for info see:
%       https://rgcgithub.github.io/regenie/options/
%       https://rgcgithub.github.io/regenie/recommendations/
%       https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1
% for gene-based tests:
%       https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/using-regenie-to-generate-variant-masks
% clc
% set input args
% 
% @03/15/2022: REGENIE 3.0 options were included for gene-/region-based tests 

arguments
    % on of bgen/bed/pgen is required
    p.bgen (1,1) {mustBeFile}
    p.bed (1,1) {mustBeFile}
    p.pgen (1,1) {mustBeFile}
    p.sample (1,1) {mustBeFile}
    p.ref_first (1,1) logical = false % (better to false for bed files, see PLINK2) Specify to use the first allele as the reference allele for BGEN or PLINK bed/bim/fam file input
    p.keep (1,1) {mustBeFile}
    p.remove (1,1) {mustBeFile}
    p.extract (1,1) {mustBeFile}
    p.exclude (1,1) {mustBeFile}
    p.phenoFile (1,1) {mustBeFile}
    p.phenoCol {mustBeText}
    p.covarFile (1,1) {mustBeFile}
    p.covarCol {mustBeText}
    p.covarColList {mustBeText}
    p.catCovarList {mustBeText}
    p.pred (1,1) {mustBeFile}
    p.interaction {mustBeTextScalar} % interaction term (GXE analysis)
    p.interaction_snp {mustBeTextScalar} % interaction term (GXG analysis)
    p.force_condtl (1,1) logical = false % to include the interacting SNP as a covariate in the marginal test 
    p.no_condtl (1,1) logical = false % to print out all the main effects from the interaction model 
    p.interaction_file {mustBeFile}
    p.interaction_file_sample {mustBeFile}
    p.rare_mac (1,1) double = 1000 % minor allele count (MAC) threshold below which to use HLM method for QTs [default is 1000]
    p.force_robust (1,1) logical = false % use robust SE instead of HLM for rare variant GxE test with quantitative traits
    p.print_prs (1,1) logical = false % flag to print whole genome predictions (i.e. PRS) without using LOCO scheme
    p.bt (1,1) logical = false % binary trait
    p.loocv (1,1) logical
	p.lowmem (1,1) logical = true % flag to reduce memory usage by writing level 0 predictions to disk (details below). This is very useful if the number of traits is large
    p.ignore_pred (1,1) logical = false % skip reading the file specified by --pred (corresponds to simple linear/logistic regression)
	p.firth (1,1) logical = false % specify to use Firth likelihood ratio test (LRT) as fallback for p-values less than threshold
	p.approx (1,1) logical = false % flag to use approximate Firth LRT for computational speedup (only works when option --firth is used)
    p.firth_se (1,1) logical = false
	p.spa (1,1) logical = true % specify to use Saddlepoint approximation as fallback for p-values less than threshold
    p.pThresh (1,1) double {mustBeInRange(p.pThresh, 1e-400, 1)} = 0.05 % P-value threshold below which to apply Firth/SPA correction [default is 0.05]
    p.test {mustBeMember(p.test, ["dominant", "recessive"])}
    p.chr (1,1) double % specify which chromosomes to test in step 2 (use for each chromosome to include)
	p.minMAC {mustBeGreaterThan(p.minMAC, 0)} = 5 % flag to specify the minimum minor allele count (MAC) when testing variants [default is 5]. Variants with lower MAC are ignored.
    p.minINFO {mustBeInRange(p.minINFO, 0, 1)} =  0.7
	p.threads (1,1) double = 64
    p.step {mustBeMember(p.step, [1 2])} = 1
    p.bsize (1,1) double =  2500
    p.verbose (1,1) logical = true 
    p.out {mustBeTextScalar} =  "regenie_out"
    p.condition_list {mustBeFile} % file with list of variants to condition on
    p.minCaseCount (1,1) double = 10 % flag to ignore BTs with low case counts [default is 10]
    % p.apply_rint (1,1) logical = true
    p.write_null_firth (1,1) logical = false
    p.use_null_firth {mustBeFile}
    
    p.af_cc (1,1) logical = true % to output A1FREQ in case/controls separately in the step 2 result file
    p.debug (1,1) logical = false

    % gene-based arguments
    % required
    p.anno_file	(1,1) {mustBeFile} % File with variant annotations for each set
    p.set_list	(1,1) {mustBeFile} % File listing variant sets
    p.mask_def (1,1) {mustBeFile} % File with mask definitions using the annotations defined in --anno-file
 
    % optional: multiple files can be specified for --extract-sets/--exclude-sets by using a comma-separated list.
    p.skat_params (1,2) double {mustBeNonnegative, mustBeNonNan} = [1,25] % a1,a2 values for the single variant weights computed from Beta(MAF,a1,a2) used in SKAT/ACAT-type tests [default is (1,25)]
    p.extract_sets (1,1) {mustBeFile} % Inclusion file that lists IDs of variant sets to keep
    p.exclude_sets (1,1) {mustBeFile} % Exclusion file that lists IDs of variant sets to remove
    p.extract_setlist {mustBeText, mustBeVector} % Comma-separated list of variant sets to keep
    p.exclude_setlist {mustBeText, mustBeVector} % Comma-separated list of variant sets to remove
    p.aaf_file (1,1) {mustBeFile} % File with variant AAF to use when building masks (instead of AAF estimated from sample)
    p.aaf_bins double {mustBeNonNan, mustBeNonnegative, mustBeVector} = [0.005, 0.01] % comma-separated list of AAF upper bounds to use when building masks [default is a single cutoff of 1%]
    p.build_mask {mustBeTextScalar, mustBeMember(p.build_mask, ["max", "sum", "comphet"])} = "max" % build masks using the maximum number of ALT alleles across sites ('max'; the default), or the sum of ALT alleles ('sum'), or thresholding the sum to 2 ('comphet')
    p.singleton_carrier (1,1) logical = true % to define singletons as variants with a single carrier in the sample (rather than alternative allele count=1)
    p.vc_tests {mustBeMember(p.vc_tests, ["skat", "skato", "skato-acat", "acatv", "acato", "acato-full"])} = ["skato", "acato-full", "skato-acat"]
    p.vc_maxAAF (1,1) double {mustBeInRange(p.vc_maxAAF, 1e-12, 1)} % = 0.01 % AAF upper bound to use for SKAT/ACAT-type tests [default is 100%]
    p.vc_MACthr (1,1) double {mustBeGreaterThan(p.vc_MACthr, 0)} = 10 % MAC threshold below which to collapse variants in SKAT/ACAT-type tests [default is 10]
    p.joint {mustBeMember(p.joint, ["minp", "acat", "sbat"])} = "acat" % comma-separated list of joint tests to apply on the generated burden masks
    p.rgc_gene_p (1,1) logical = true % to compute the GENE_P test
    p.skip_test (1,1) logical % to skip computing association tests after building masks and writing them to file
    p.mask_lovo {mustBeTextScalar} % to perform LOVO scheme e.g. "APOB,M1,0.01"
    p.mask_lodo (1,1) logical % to perform LODO scheme
    p.write_mask_snplist (1,1) logical = true % to write list of variants that went into each mask to file
    p.check_burden_files (1,1) logical = true % to check the concordance between annotation, set list and mask files
    p.strict_check_burden (1,1) logical % to exit early if the annotation, set list and mask definition files dont agree
    p.write_mask (1,1) logical = false % write mask to PLINK bed format (does not work when building masks with 'sum')
    
    p.weights_col

    %@16NOV2024: survival analysis v4.0 or above is needed
    p.t2e (1,1) logical = false
    p.phenoColList {mustBeText, mustBeVector} % Comma separated list of time names to include in the analysis
    p.eventColList {mustBeText, mustBeVector} % Comma separated list of columns in the phenotype file to include in the analysis that contain the events. These event columns should have 0=no event,1=event,NA=missing
    p.coxnofirth (1,1) logical
    p.coxscore_exact (1,1) logical

    % non-native arguments
    p.runcmd (1,1) logical = false
    p.regenie {mustBeTextScalar} = "regenie_v3_4_1" % executable regenie file
end

% Check files' paths and convert to wsl paths if needed -------------------
checkwslpath = ["bgen", "bed", "pgen", "keep", "remove", "exclude", ...
    "sample", "exclude", "phenoFile", "covarFile", "pred", "anno_file", ...
    "set_list", "mask_def", "extract_sets", "exclude_sets", "aaf_file", ...
    "condition_list", "extract", "interaction_file", "interaction_file_sample",...
    "use_null_firth"];
for ii = 1:numel(checkwslpath)
    if isfield(p, checkwslpath(ii)) && ispc
        p.(checkwslpath(ii)) = makeWSLpath(p.(checkwslpath(ii)));
        if ~ismember(checkwslpath(ii), ["bgen", "bed", "pgen", "sample", "interaction_file", "interaction_file_sample"])
            dos2unix(p.(checkwslpath(ii)));
        end
    end
end

if ispc, p.out = makeWSLpath(p.out); end

% reset gene-based arguments if unused ------------------------------------
if ~isfield(p, 'anno_file')
    rmfields = ["write_mask_snplist", "check_burden_files", ...
        "joint", "write_mask", "singleton_carrier", "vc_MACthr", ...
        "vc_maxAAF",  "vc_tests", "singleton_carrier", "build_mask",...
        "aaf_bins", "skat_params", "rgc_gene_p"];
    rmfields = intersect(rmfields, fieldnames(p));
    p = rmfield(p, rmfields);
end

% Check BGEN/BED/PGEN preferences -----------------------------------------
[bedFlag, pgenFlag, bgenFlag] = deal(false);
if isfield(p, 'bgen') && ~isempty(p.bgen) % ignore other pgen/bgen
    bgenFlag = true;
    p = rmfield(p, intersect(fieldnames(p), ["bed",  "pgen"]));
elseif isfield(p, 'pgen') && ~isempty(p.pgen)
    p.pgen = regexprep(p.pgen, '.pgen$', '');
    pgenFlag = true;
    p = rmfield(p, intersect(fieldnames(p), ["bed",  "bgen", "sample"]));
elseif isfield(p, 'bed') && ~isempty(p.bed)
    p.bed = regexprep(p.bed, '.bed$', '');
    bedFlag = true;
    p = rmfield(p, intersect(fieldnames(p), ["bgen",  "pgen", "sample"]));
end

if ~any([bedFlag, pgenFlag, bgenFlag]) % nothing is provided
    error('regenie:either bed,pgen or bgen is required!')
end

% DEPRECATED: user dependent
% askFlag = false;

% % only 1 data type must be fed into the software
% if sum([bedFlag, pgenFlag, bgenFlag]) > 1
%     error('only bed or pgen or bgen should be given as input!')
% elseif ~any([bedFlag, pgenFlag, bgenFlag]) % nothing is provided
%     askFlag = true;
%     inChoice = input('1-BED 2-PGEN 3-BGEN?[3] ');
%     if isempty(inChoice) || inChoice == 3
%         bgenFlag = true;
%         askText = 'Choose bgen file(s)';
%         askFilter = '*.bgen';
%     elseif inChoice == 1
%         bedFlag = true;
%         askText = 'Choose bed file(s)';
%         askFilter = '*.bed';
%     elseif inChoice == 2
%         pgenFlag = true;
%         askText = 'Choose pgen file(s)';
%         askFilter = '*.pgen';
%     else
%         fprintf('unknown input, set to BGEN!\n')
%         bgenFlag = true;
%     end
% end
% 
% if askFlag % not set in the input -----------------------------------------
%     [inFile, inFilePath] = uigetfile(askFilter, askText, 'MultiSelect', 'on');
%     if ~inFilePath
%         fprintf('No file has been selected!\n')
%         return
%     end
%     if bedFlag
%         p.bed = makeWSLpath(inFile, inFilePath);
%         p.bed = regexprep(p.bed, '.bed$', '');
%         inFiles = p.bed;
%     elseif pgenFlag
%         p.pgen = makeWSLpath(inFile, inFilePath);
%         p.pgen = regexprep(p.pgen, '.pgen$', '');
%         inFiles = p.pgen;
%         
%     elseif bgenFlag
%         p.bgen = makeWSLpath(inFile, inFilePath);
%         p.sample = regexprep(p.bgen, '.bgen$', '.sample');
%         p.with_bgi = true;
%         inFiles = p.bgen;
%         
%     else
%         error('unknown input file type!')
%     end
% end

%% Prepare REGENIE command ------------------------------------------------
regenieDir = fileparts(which('regenie.m'));
if ispc; regenieDir = makeWSLpath(regenieDir); end
runcmd = p.runcmd;
p = rmfield(p, 'runcmd');

cmd = makeregeniecmd(regenieDir, p);
if ispc; cmdpref = "wsl "; else; cmdpref = ""; end
if runcmd
    if p.verbose
        system(cmdpref + cmd);
    else
        [~, ~] = system(cmdpref + cmd);
    end
end

end % END

%% subfunctions ----------------------------------------------------------
function cmd = makeregeniecmd(regenieDir, p)
% cnt: file index

cmd = '"' + string(regenieDir) + "/" + p.regenie + '"';
p = rmfield(p, 'regenie');
fd = string(fieldnames(p));

for i = 1:numel(fd)
    if isempty(p.(fd(i)))
        continue
    end

    if islogical(p.(fd(i))) % FLAG
        if p.(fd(i))
            cmd = join([cmd, "--" + replace(fd(i), '_', '-')], " ");
        end
    else
        if numel(p.(fd(i))) > 1 % input files
            tmpval =  join(string(p.(fd(i))), ",");
        else
            tmpval = string(p.(fd(i)));
        end

        if contains(tmpval, ["/", "\"]) % ? path to file
            if fd(i) == "interaction_file", [~,~,ext] = fileparts(tmpval); end
            tmpval = '"' + tmpval + '"';
            if fd(i) == "interaction_file"
                tmpval = regexprep(ext, "^\.", "") + "," + tmpval;
            end
        end

        cmd = join([cmd, "--" + replace(fd(i), '_', '-'), tmpval], ' ');
        
    end
end
end % END
