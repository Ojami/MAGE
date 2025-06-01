function cmd = call_regenie(opts)
%@16APR2023: extract_setlist and exclude_setlist can be used now

%@12MARCH2024: 'force_robust' and 'rare_mac' options were added for
%               interaction analysis. These two options make it possible to
%               perform interaction using robust SE, HLM or mixture of
%               both. By default, 'force_robust' is false, meaning that for
%               variants > rare_mac, robust SE is used, while HLM is used
%               for rare variants (< rare_mac). If 'force_robust' is true,
%               then HLM is ignored (ergo 'rare_mac') and only robust SE
%               (HC3 sandwich estimator) is used. If one opts for only HLM,
%               then 'force_robust' should be false (default) and
%               'rare_mac' should be set to a fairly large value.
% @12JULY2024: 'annoChunk' option was added to split anno files into
%               separate chunks. By default, uses 'thread' to decide.
% @16NOV2024: '--t2e' to TTE analysis (>= v4.0) was added.

arguments
    % non-native arguments
    opts.trait {mustBeText, mustBeVector} = "" % if left empty, UKBtrait_toolbox will fetch traits
    opts.qc (1,1) double {mustBeMember(opts.qc, 0:6)} = 0 % see getQCEID function for more details
    opts.qcinc (:, 1) double {mustBeNumeric} = [] % custom eids to be included; in this case intersect of qc and qcinc will be used for analyses.
    opts.phenopath {mustBeFolder} % phenotype data should be here
    opts.transform {mustBeMember(opts.transform, ["none", "log", "log10", "irnt"])} = "irnt" % transformation function for continuous variables.
    opts.covar {mustBeText} = ["Age", "Sex", "BMI"]; % list of covariates, it automatically adds PC1-10 and array batch.
    opts.removenan (1,1) logical = true % remove samples with missing phenotype (only effective is number of query traits is 1).
    opts.agesexInt(1,1) logical = false; % to add age^2, age*sex and sex*age^2 to covariates (only if already present in covariate list).
    % opts.findphenocov (1,1) logical = true % for step 2, don't get phenotypes/covariates again, use the one generated in step 1
    opts.keepLogFiles (1,1) logical = true % to keep/remove log files from step 2
    opts.bgenhome {mustBeFolder} % for step 2, all bgen files must be present in this folder
    opts.wes (1, 1) logical = false % for WES analysis (array batch is removed from covariates). DON'T set to true in step 1 if directly genotyped files are used.
    opts.weshome {mustBeFolder}  % for step 2 gene-based tests, all bgen/bed files must be present in this folder
    opts.verbose (1,1) logical = true

    % to split the jobs over each part of anno_files (see geneWrapper doc).
    % If true (recommended), runs REGENIE jobs per each part of
    % 'anno_file'; this is only applicable when size of 'anno_file' > size
    % of geno files ('bed'/'bgen'/'pgen'). If false (not recommended), runs
    % all jobs at the same time. This is not recommended because multiple
    % calls to geno files renders the system in out of memory situation.
    opts.multi_anno (1,1) logical = true
    opts.regenie {mustBeTextScalar} = "regenie_v4_1" % executable regenie file
    opts.annoChunk (1,1) double % number of chunks to split
    

    % REGENIE arguments
    opts.phenoFile {mustBeFile} % if provided, skips processing trait (e.g. for step 2 when you already have generated the phenoFile). Function looks for bt.mat file in the same path (generated in step 1)
    opts.covarFile {mustBeFile} % 'phenoFile' overrides this option, i.e. use an already generated covarFile _only if_ 'phenoFile' argument is used, otherwise a new covarFile will be generated.
    opts.bed {mustBeFile}  % for step 1, must be a single bed file
    opts.bgen {mustBeFile}
    opts.pgen {mustBeFile}
    opts.sample {mustBeFile}
    opts.step {mustBeMember(opts.step, [1 2])} = 1 % REGENIE step 
    opts.pred {mustBeFile} % for step 2 
    opts.out {mustBeTextScalar} % file output prefix; by default pheno name is used for output
    opts.threads (1,1) double = 64
    opts.keep {mustBeFile}
    opts.remove {mustBeFile}
    opts.extract {mustBeFile}
    opts.exclude {mustBeFile}
    opts.minINFO (1,1) double {mustBeInRange(opts.minINFO, 0, 1)} = 0
    opts.firth (1,1) logical = false % specify to use Firth likelihood ratio test (LRT) as fallback for p-values less than threshold
	opts.approx (1,1) logical = false % flag to use approximate Firth LRT for computational speedup (only works when option --firth is used)
    opts.firth_se (1,1) logical = false
    opts.write_null_firth (1,1) logical = false
    opts.use_null_firth {mustBeFile}
    opts.print_prs (1,1) logical = false % flag to print whole genome predictions (i.e. PRS) without using LOCO scheme
    opts.bt (1,1) logical % binary trait
    opts.spa (1,1) logical % specify to use Saddlepoint approximation as fallback for p-values less than threshold
    opts.pThresh (1,1) double {mustBeInRange(opts.pThresh, 1e-400, 0.9999999999)} = 0.05 % P-value threshold below which to apply Firth/SPA correction [default is 0.05]
    opts.bsize (1,1) double % block size (if left empty, 2500 for step1 and 8000 fro step2 will be used)
    opts.rare_mac (1,1) double = 1000 % minor allele count (MAC) threshold below which to use HLM method for QTs [default is 1000]
    opts.minCaseCount (1,1) double = 10 % flag to ignore BTs with low case counts [default is 10]
    opts.minMAC (1,1) double

    % gene-based arguments
    % required
    opts.anno_file	{mustBeFile} % File with variant annotations for each set
    opts.set_list	{mustBeFile} % File listing variant sets
    opts.mask_def {mustBeFile} % File with mask definitions using the annotations defined in --anno-file

    % optional: multiple files can be specified for --extract-sets/--exclude-sets by using a comma-separated list.
    opts.skat_params (1,2) double {mustBeNonnegative, mustBeNonNan} = [1,25] % a1,a2 values for the single variant weights computed from Beta(MAF,a1,a2) used in SKAT/ACAT-type tests [default is (1,25)]
    opts.extract_sets {mustBeFile} % Inclusion file that lists IDs of variant sets to keep
    opts.exclude_sets {mustBeFile} % Exclusion file that lists IDs of variant sets to remove
    opts.extract_setlist {mustBeText, mustBeVector} % Comma-separated list of variant sets to keep
    opts.exclude_setlist {mustBeText, mustBeVector} % Comma-separated list of variant sets to remove
    opts.aaf_file {mustBeFile} % File with variant AAF to use when building masks (instead of AAF estimated from sample)
    opts.aaf_bins double {mustBeNonNan, mustBeNonnegative, mustBeVector} = 0.01 %e.g [0.05, 0.01] comma-separated list of AAF upper bounds to use when building masks [default is a single cutoff of 1%]
    opts.build_mask {mustBeTextScalar, mustBeMember(opts.build_mask, ["max", "sum", "comphet"])} = "max" % build masks using the maximum number of ALT alleles across sites ('max'; the default), or the sum of ALT alleles ('sum'), or thresholding the sum to 2 ('comphet')
    opts.singleton_carrier (1,1) logical = true % to define singletons as variants with a single carrier in the sample (rather than alternative allele count=1)
    opts.vc_tests {mustBeMember(opts.vc_tests, ["skat", "skato", "skato-acat", "acatv", "acato", "acato-full"])} = ["skato", "acato-full", "skato-acat"]
    opts.vc_maxAAF (1,1) double {mustBeInRange(opts.vc_maxAAF, 1e-12, 1)} = 0.01 % AAF upper bound to use for SKAT/ACAT-type tests [default is 100%]
    opts.vc_MACthr (1,1) double {mustBeGreaterThan(opts.vc_MACthr, 0)} = 10 % MAC threshold below which to collapse variants in SKAT/ACAT-type tests [default is 10]
    opts.joint {mustBeMember(opts.joint, ["minp", "acat", "sbat"])} = "acat" % comma-separated list of joint tests to apply on the generated burden masks
    opts.rgc_gene_p (1,1) logical = true % to compute the GENE_P test
    opts.skip_test (1,1) logical % to skip computing association tests after building masks and writing them to file
    opts.mask_lovo {mustBeTextScalar} % to perform LOVO scheme e.g. "APOB,M1,0.01"
    opts.mask_lodo (1,1) logical % to perform LODO scheme
    opts.write_mask_snplist (1,1) logical = true % to write list of variants that went into each mask to file
    opts.check_burden_files (1,1) logical = true % to check the concordance between annotation, set list and mask files
    opts.strict_check_burden (1,1) logical % to exit early if the annotation, set list and mask definition files dont agree
    opts.write_mask (1,1) logical = false % write mask to PLINK bed format (does not work when building masks with 'sum')
    opts.debug (1,1) logical = false
    opts.catCovarList {mustBeTextScalar} % comma separated list of categorical covariates. For now only supports when called from gwasrunner or 'covarFile' is used.

    % under development
    opts.condition_list {mustBeFile} % file with list of variants to condition on
    opts.interaction {mustBeTextScalar} = '' % interaction term (GXE analysis). 
    opts.interaction_snp {mustBeTextScalar} = ''% interaction term (GXG analysis).
    opts.force_condtl (1,1) logical = false % to include the interacting SNP as a covariate in the marginal test 
    opts.no_condtl (1,1) logical = false % to print out all the main effects from the interaction model
    opts.interaction_file {mustBeFile}
    opts.interaction_file_sample {mustBeFile}
    opts.loocv (1,1) logical
    opts.force_robust (1,1) logical = false % use robust SE instead of HLM for rare variant GxE test with quantitative traits

    %@01NOV2024
    opts.weights_col (1,1) double % custom weight column number for rare variant test

    %@16NOV2024: survival analysis v4.0 or above is needed
    opts.t2e (1,1) logical = false
    opts.phenoColList {mustBeText, mustBeVector} % Comma separated list of time names to include in the analysis
    opts.eventColList {mustBeText, mustBeVector} % Comma separated list of columns in the phenotype file to include in the analysis that contain the events. These event columns should have 0=no event,1=event,NA=missing
    opts.coxnofirth (1,1) logical
    opts.coxscore_exact (1,1) logical

    %@08MAY2025
    opts.cmdOnly (1,1) logical = false % doesn't run anything, just returns the executable ready commands
end


% prepare variant calls data ----------------------------------------------
if isfield(opts, 'anno_file') % gene-based tests
    opts.wes = true;
    if ~isfield(opts, "minMAC"), opts.minMAC = 5; end % default is 5
end

% check step specific inputs ==============================================
if opts.step == 1 % bed files is required for this step -------------------
    if ~isfield(opts, 'bed') % cannot left empty
        % try to find it in current wd
        opts.bed = getfilenames(pwd, "bed").bed;
        if ~isempty(opts.bed) && isscalar(opts.bed)
            opts.bed = fullfile(pwd, opts.bed);
            if opts.verbose, fprintf('found a bed file: %s\n', opts.bed); end
        else
            error('call_regenie: for step 1, bedfile input is required!')
        end
    end
    if ~isfile(opts.bed)
        error('call_regenie: %s not a valid file!', opts.bed)
    end
    % opts.bed = regexprep(opts.bed, '.bed$', '');
else % step 2 with bgen files ---------------------------------------------
    
    % check prediction file presence
    if ~isfield(opts, 'pred') % no prediction file
%         % search within this directory
%         opts.pred = getfilenames(pwd, "list").list;
%         if isempty(opts.pred) % ignore pred file
%             if opts.verbose, fprintf('couldn''t find pred file, skipping loco file!\n'); end
%             opts = rmfield(opts, 'pred');
            opts.ignore_pred = true;
%         else
%             opts.pred = strtrim(string(opts.pred));
%             if numel(opts.pred) > 1
%                 error('call_regenie: only 1 prediction (.list) file is allowed!')
%             end
%             if opts.verbose, fprintf('found prediction file: %s\n', opts.pred); end
%             opts.ignore_pred = false;
%         end
    end
    
    % decide between GWAS and gene-based analysis
    if isfield(opts, 'anno_file') % gene-based tests
        if opts.verbose, fprintf('gene-based tests will be performed\n\n'); end
        if ~isfield(opts, "bed")
            opts.bed = getfilenames(opts.weshome, "bed", 'fullpath', true).bed;
             if isempty(opts.bed) % no bed file was found in wes home, try bgen
                 if ~isfield(opts, 'bgen')
                     opts.bgen = getfilenames(opts.weshome, "bgen", 'fullpath', true).bgen;
                     if isempty(opts.bgen)
                        error('call_regenie: no bgen/bed was file found in %s', opts.weshome)
                     end
                 end
             end
        end
    else % GWAS
        if ~any(isfield(opts, ["bgen", "pgen"]))
            opts.bgen = getfilenames(opts.bgenhome, "bgen", 'fullpath', true).bgen;
            if isempty(opts.bgen) % bed files are there?
                if ~isfield(opts, 'bed')
                    opts.bed = getfilenames(opts.bgenhome, "bed", 'fullpath', true);
                    if isempty(opts.bed)
                        error('call_regenie: no bgen/bed files were found in %s', opts.bgenhome)
                    else
                        opts.bed = opts.bed.bed;
                    end
                end
                opts = rmfield(opts, 'bgen');
            end
        end
    end
end

if isfield(opts, 'bgen')
    opts.sample = regexprep(opts.bgen, '.bgen$', '.sample');
end

% get phenotype and covariates --------------------------------------------
if opts.verbose, fprintf('getting/writing phenotypes and covariates...\n'); end

if isfield(opts, 'phenoFile') % use already generated phenoFile (and/or covarFile)
    
    [~, traitTag, phenoExt] = fileparts(opts.phenoFile);
    traitTag = regexprep(traitTag, ".phenofile$", "", 'ignorecase');
    if ~isfield(opts, 'out'); opts.out = traitTag; end
    
    % check if binary/continuous flag exists
    btFile = regexprep(opts.phenoFile, phenoExt + "$", ".btFlag.mat");

    if isfile(btFile)
        opts.bt = load(btFile).bt;
    elseif ~isfield(opts, 'bt') % infer it from phenotype file
        checkbt = detectImportOptions(opts.phenoFile, ...
            'FileType', 'text', 'ReadVariableNames', ...
            true, 'NumHeaderLines', 0);
        checkbt.SelectedVariableNames(1:2) = []; % IID/FID
        checkbt = readtable(opts.phenoFile, checkbt);
        for l = 1:width(checkbt)
            if numel(unique(checkbt.(l))) > 2
                opts.bt(l) = false;
            else
                opts.bt(l) = true;
            end
        end
    end

else % generate phenoFile/covarFile
    [response, cov] = getPhenoCov('trait', opts.trait, 'qc', opts.qc,...
        'covar', opts.covar, 'phenopath', opts.phenopath, 'qcinc', opts.qcinc, ...
        'removenan', opts.removenan, 'transform', opts.transform, 'wes', opts.wes);
    response.Properties.VariableNames = matlab.lang.makeValidName(response.Properties.VariableNames);
    cov.Properties.VariableNames = matlab.lang.makeValidName(cov.Properties.VariableNames);
    
    assert(all(cov.eid == response.eid), "all eid in response/cov tables should be in the same order!")

    %@16NOV2024: survival analysis
    if opts.t2e
        tmp = load(opts.trait).UKB_STRUCT_ALL;
        idx = ~ismember(response.eid, tmp.eid); % to be removed
        response(idx, :) = [];
        cov(idx, :) = [];

        if isempty(response)
            error("no individual left for analysis!")
        end
        
        rcol = setdiff(colnames(response), "eid"); % trait column
        tcol = rcol + "_TIME"; % tt col
        response.(rcol)(:) = 0;
        response.(rcol)(ismember(response.eid, tmp.eid(tmp.censor == 1))) = 1; % cause-specific hazard
        
        % add time-to-event
        [~, idx] = ismember(response.eid, tmp.eid); 
        response.(tcol) = tmp.tt(idx);

        opts.phenoColList = tcol;
        opts.eventColList = rcol;
    end
    
    % add age/sex interaction to covariates
    if opts.agesexInt
        % age^2 + sex*age + sex*age2
        ageIdx = find(contains(lower(cov.Properties.VariableNames), 'age'));
        sexIdx = find(contains(lower(cov.Properties.VariableNames), {'sex', 'gender'}));
        if isscalar(ageIdx) 
        	ageTag = string(cov.Properties.VariableNames{ageIdx})+"_2";
            cov.(ageTag) = cov.(ageIdx).^2;
        end
        if isscalar(sexIdx)
            agesexTag = string(cov.Properties.VariableNames{sexIdx})+"_"+string(cov.Properties.VariableNames{ageIdx});
            cov.(agesexTag) = cov.(sexIdx).*cov.(ageIdx);
            agesexTag = string(cov.Properties.VariableNames{sexIdx})+"_"+string(cov.Properties.VariableNames{ageIdx})+"_2";
            cov.(agesexTag) = cov.(sexIdx).*(cov.(ageIdx).^2);
        end
    end

    % write response to pheno file
    [response.FID, response.IID] = deal(response.eid);
    response.eid = [];
    traitTag = string(response.Properties.VariableNames{1});
    response = movevars(response, {'FID', 'IID'}, 'Before', 1);
    
    %@15OCT2024: a bug was fixed for absolute path in 'out'
    if ~isfield(opts, 'out')
        opts.out = traitTag; 
    elseif fileparts(opts.out) == "" % only pheno tag in current directory
        opts.out = fullfile(pwd, opts.out);
    else % subfolder in current directory, get absolute path
        [rpath, rtag] = fileparts(opts.out);
        if ~isfolder(rpath), mkdir(rpath); end
        [~, info] = fileattrib(rpath);
        opts.out = string(info.Name);
        if rtag == "", rtag = traitTag; end
        opts.out = fullfile(opts.out, rtag);
    end

    writetable(response, opts.out + ".phenoFile.txt", 'Delimiter', '\t')
    opts.phenoFile = fullfile(opts.out + ".phenoFile.txt");

    opts.bt = response.Properties.UserData;
    if numel(unique(opts.bt)) > 1
        error('call_regenie:binaryCont', 'all phenotypes must be of the same type!')
    else
        opts.bt = unique(opts.bt);
    end
    bt = opts.bt;
    save(opts.out + ".phenoFile.btFlag.mat", 'bt')
    response(:, 3:end) = []; % keep it for covariates

    if ~isempty(cov)
        cov.eid = [];
        response = [response, cov];
        writetable(response, opts.out + ".covarFile.txt", 'Delimiter', '\t');
        clear response cov
        opts.covarFile = fullfile(opts.out + ".covarFile.txt");
    end
end

if opts.t2e % change 'out' argument
    opts.out = opts.out + "_TIME";
end

% -------------------------------------------------------------------------
if opts.verbose, fprintf('running REGENIE step %d\n', opts.step); end

% use 'spa' if binary traits
if ~isfield(opts, 'spa')
    if opts.bt
        opts.spa = true;
    else 
        opts.spa = false;
    end
end

if isfield(opts, 'firth') && opts.firth
    opts.spa = false;
end

% non-REGENIE arguments
rmfis = ["trait", "qc", "qcinc", "phenopath", "transform", "covar", ...
    "removenan", "agesexInt", "keepLogFiles", "cmdOnly", ...
    "bgenhome", "wes", "weshome", "multi_anno", "verbose", "annoChunk"];
newopts = opts;
newopts = rmfield(newopts, intersect(fieldnames(opts), rmfis));

% set some arguments to their default values
newopts.runcmd = false;

if isfield(opts, 'bed')
    sz = numel(opts.bed);
elseif isfield(opts, 'bgen')
    sz = numel(opts.bgen); 
else
    sz = numel(opts.pgen);
end

% only when there are multiple parts per each 'anno_file'
if isfield(opts, 'anno_file')
    if numel(opts.anno_file) <= sz
        opts.multi_anno = false;
    end
else
    opts.multi_anno = false; % only if 'anno_file' is used
end

% check size of annotation and set files for gene-based tests
if isfield(opts, 'anno_file') && numel(opts.anno_file) > sz && ~opts.multi_anno % multiple anno file per chromosome
    sz = numel(opts.anno_file); % should be same as set/mask size
    opts.multi_anno_all = true; % run all ".part" files at once
else
    opts.multi_anno_all = false;
    if isfield(opts, 'anno_file') && numel(opts.anno_file) < sz % only for a subset of chromosomes
        opts = subsetGenoFiles(opts); 
        sz = numel(opts.anno_file);
    end
end

opts.outpath = pwd; % where temp files should be saved (if it's called from another function)

% check if results should be written to a another dir than pwd
outdir = string(fileparts(opts.out));
if outdir ~= "" && ~isfolder(outdir)
    mkdir(outdir)
end

if opts.step == 1
    if isfield(opts, 'bsize')
        newopts.bsize = opts.bsize;
    else
        newopts.bsize = 2500;
    end
    newopts.runcmd = ~opts.cmdOnly; %@08MAY2025
    newopts = namedargs2cell(newopts);

    if opts.cmdOnly
        cmd = regenie(newopts{:});
        return
    else
        regenie(newopts{:});
    end
else % step 2 -------------------------------------------------------------
    if isfield(opts, 'bsize')
        newopts.bsize = opts.bsize;
    else
        newopts.bsize = 8e3;
    end

    % to run all parts ("part") of 'anno_file' at once ('multi_anno'
    % should be false in this case)
    if isfield(opts, 'anno_file') && (opts.multi_anno_all || opts.multi_anno) % multiple anno file per chromosome   
        % replicate bed/bgen files to have the same size as of 'anno_file'
        chrvec.str = string(regexp(opts.anno_file, "[.]chr(.*?)[.]", 'tokens'));
        [~, chrvec.count] = duplicates(chrvec.str);

        if isfield(opts, 'bed')
            opts.bed = repelem(opts.bed, chrvec.count);
        else % bgen (pgen should be implemented)
            opts.bgen = repelem(opts.bgen, chrvec.count);
            opts.sample = repelem(opts.sample, chrvec.count);
        end
    end

    if opts.multi_anno % get parts of each annotation file 
        parts = regexp(opts.anno_file, ".part\d+.", 'match');
        if all(cellfun(@isempty, parts))
            error('call_regenie:no "part" pattern was found in anno_file!')
        end
        parts = unique(string(parts));
        [cmd, extsetfiles] = deal(cell(numel(parts), 1));
        for j = 1:numel(parts)
            partopts = opts;
            part_idx = contains(partopts.anno_file, parts(j));
            if sum(part_idx) < sz % only for a subset of chromosomes
                sz = sum(part_idx);
            elseif sum(part_idx) > sz % not allowed!
                error('call_regenie: part size cannot be > geno file size')
            end

            partopts.anno_file = partopts.anno_file(part_idx);
            partopts.set_list = partopts.set_list(part_idx);
            if numel(partopts.mask_def) > 1
                partopts.mask_def = partopts.mask_def(part_idx);
            end
            if isfield(opts, 'bed')
                partopts.bed = partopts.bed(part_idx);
            else % bgen (pgen should be implemented)
                partopts.bgen = partopts.bgen(part_idx);
                partopts.sample = partopts.sample(part_idx);
            end
            [cmd{j}, extsetfiles{j}] = getREGENIEcmd(partopts, newopts, sz, parts(j));
        end
        cmd = vertcat(cmd{:});
        extsetfiles = vertcat(extsetfiles{:});
    else
        [cmd, extsetfiles] = getREGENIEcmd(opts, newopts, sz);
    end

    cmd = vertcat(cmd{:});

    %@08MAY2025
    if opts.cmdOnly, return; end

    [pth, outname, ext] = fileparts(opts.out);
    if pth ~= "", opts.outpath = pth; end
    outname = outname + ext;

    if ~isfield(opts, 'anno_file')
        if opts.verbose, t1 = tic; end
        for i = 1:numel(cmd)
            tmpout = extractBetween(cmd(i), "--out ", " ").erase('"') + ".log"; % check if log file exists (skip if so)
            tmpout = erase(tmpout, "");
            if ispc
                if isfile(makeWSLpath(tmpout, true)), continue; end
                if opts.verbose
                    system("wsl " + cmd(i));
                else
                    [~, ~] = system("wsl " + cmd(i));
                end
                % runbash(cmd(i), "regenie_temp", "conda", "regenie_env", ...
                %     "parallel", false, "verbose", opts.verbose)
            else
                if isfile(tmpout), continue; end
                if opts.verbose
                    system(cmd(i));
                else
                    [~, ~] = system(cmd(i));
                end
            end
        end
        if opts.verbose, toc(t1), end
    else

        if opts.multi_anno % memory-friendly jobs 
            seq = [1:opts.threads:numel(cmd), numel(cmd) + 1];
            if opts.verbose, t2 = tic; end
            for i = 1:numel(seq)-1
                runbash(cmd(seq(i):seq(i+1)-1), outname + ".call.regenie", ...
                    'parallel', true, 'dir', opts.outpath);
            end
            if opts.verbose, toc(t2), end
        else
            if opts.verbose, t2 = tic; end

             runbash(cmd, outname + ".call.regenie", ...
                 'parallel', true, 'dir', opts.outpath);
             if opts.verbose, toc(t2), end
        end
    end
    
end

% debug_names = getfilenames(opts.outpath).x_;
% disp(debug_names)

% Merge summary stats from different bgen files
if fileparts(opts.out) == "", opts.out = fullfile(pwd, opts.out); end
opts.mergedfile = opts.out + ".sumstat.step2.txt";
if opts.step == 2 % && sz > 1

    if opts.verbose, fprintf('final step: merging summary stats from different chromosomes...'); end
    summFiles = getfilenames(opts.outpath, ".regenie", fullpath=true).regenie;
    logFiles = getfilenames(opts.outpath, ".log", fullpath=true).log;
    
    if ~isempty(logFiles)
        logFiles(~startsWith(logFiles, opts.out)) = [];
        if opts.keepLogFiles
             mergeStatFiles(logFiles, opts.out + ".log");
        end
        logFiles = setdiff(logFiles, opts.out + ".log");
        delete(logFiles{:})
    end

    % check if REGENIE faced any error
    if isfile(opts.out + ".log")
        tmp_log = readlines(opts.out + ".log");
        tmp_log(tmp_log == "") = [];
        if tmp_log(end-3:end).lower.contains("error")
         error("REGENIE faced an error, see the log file!")
        end
    end

    if isfield(opts, 'anno_file') % gene-based tests
        n = "+3 ";
        snpListFiles = getfilenames(opts.outpath, "snplist", fullpath=true);

        % delete temp extract_set files
        if ~isempty(extsetfiles)
            extsetfiles = vertcat(extsetfiles{:});
            delete(extsetfiles{:})
        end

        if ~isempty(snpListFiles)
            snpListFiles = snpListFiles.snplist;
            snpListFiles(~startsWith(snpListFiles, opts.out)) = [];
            if ~isempty(snpListFiles)
                mergeStatFiles(snpListFiles, opts.out + ".snplist.txt");
                delete(snpListFiles{:});
            end
        end
        masksReportFiles = getfilenames(opts.outpath, "txt", fullpath=true).txt;
        if ~isempty(masksReportFiles)
            idx = startsWith(masksReportFiles, opts.out) & endsWith(masksReportFiles, "masks_report.txt");
            masksReportFiles(~idx) = [];
            if ~isempty(masksReportFiles)
                mergeStatFiles(masksReportFiles, opts.out + ".masks_report.txt");
                delete(masksReportFiles{:});
            end
        end
    else
        n = "+2 ";
    end

    mergeStatFiles(summFiles, opts.mergedfile, n)  
    delete(summFiles{:})
    if opts.verbose, fprintf('\b\b Done.\n'); end
end

% unstack each combination of TEST and ALLELE1 (mask/aaf combinations) and
% store to a new table 
if isfield(opts, 'anno_file') && isfile(opts.mergedfile)
    readRegenieGeneStat(opts.mergedfile)
end

% check if different tests are present in the summary stat file. If so,
% write the results of each test to a separate file (mainly for interaction
% GWAS)
if isfile(opts.mergedfile) && ~isfield(opts, 'anno_file') && opts.verbose
    checksum = bfilereader(opts.mergedfile, 'summary', 'only', 'header', true);
    if numel(unique(checksum.TEST)) > 1
        if isempty(gcp('nocreate')); parpool('threads'); end
        tests = unique(checksum.TEST);
        ds = tabularTextDatastore(opts.mergedfile, ...
            'NumHeaderLines', 0, 'TextType', 'string', ...
            'TreatAsMissing', 'NA', 'CommentStyle', '##');
        cols = ds.SelectedVariableNames;
        cols(ismember(cols, {'INFO', 'EXTRA'})) = [];
        ds.SelectedVariableNames = cols;
        ds = tall(ds);
    
        for i = 1:numel(tests)
            if opts.verbose, fprintf('processing/writing for %s (%d of %d files)\n', tests(i), i, numel(tests)); end
            idx = ds.TEST == tests(i);
            tab = gather(ds(idx, :));
            tab.TEST = [];
            name = replace(tests(i), ["_", "-"], ".") + ".txt";
            writetable(tab, name, 'Delimiter', '\t')
            clear tab idx
        end
    
        clear ds
    %     delete('REGENIE.summStat.step2.txt')
    end
end

end % END

%% subfunctions ===========================================================
function opts = createSetChunks(opts, n)
[~, name, ext] = fileparts(opts.set_list);
hm = fileparts(opts.out);
if hm == "", hm = pwd; end
tmpname = @(x) sprintf(name + ".tmp%d" + ext, x);
if isfile(fullfile(hm, tmpname(1))) % already exists!
    ct = 1;
    while true
        opts.extract_sets(ct, 1) = fullfile(hm, tmpname(ct));
        ct = ct + 1;
        if ~isfile(fullfile(hm, tmpname(ct))), break; end
    end
    return
end
genes = bfilereader(opts.set_list, 'extractCol', 1, 'verbose', 'off').(1);
genes(genes == "") = [];

%@16APR2023: extract_setlist and exclude_setlist can be used now
if isfield(opts, "extract_setlist")
    tgenes = opts.extract_setlist;
    tgenes = erase(split(tgenes, ","), '"');
    genes = intersect(tgenes, genes);
    if isempty(genes)
        error("regenie:sets in extract_setlist cannot be found in set files!")
    end
end

if isfield(opts, "exclude_setlist")
    tgenes = opts.exclude_setlist;
    tgenes = erase(split(tgenes, ","), '"');
    genes = setdiff(genes, tgenes);
    if isempty(genes)
        error("regenie:no set remained after applying exclude_setlist!")
    end
end

idx = unique(round(linspace(1, numel(genes) + 1, n + 1))); 
idx(end) = numel(genes) + 1;

for i = 1:numel(idx)-1
    opts.extract_sets(i, 1) = fullfile(hm, tmpname(i));

    writematrix(genes(idx(i):idx(i+1)-1,:), fullfile(hm, tmpname(i)),...
        'QuoteStrings', false, 'FileType', 'text');

end
end % END

%% ------------------------------------------------------------------------
function opts = subsetGenoFiles(opts)
% chrsubset = regexp(opts.anno_file, "[.]chr(.*?)[.]", 'tokens'); % sex chromosomes should be numbered sequentially
% chrsubset(cellfun(@isempty, chrsubset)) = [];
% if isempty(chrsubset) % files not numbers, peek inside!
chrsubset = nan(numel(opts.anno_file), 1);
for i = 1:numel(opts.anno_file)
    fid = fopen(opts.anno_file(i), 'r');
    tmp = string(fgetl(fid));
    fclose(fid);
    tmp = split(tmp);
    chrsubset(i, 1) = double(string(regexp(tmp(1), '^(\d+)?', 'tokens'))); % only works if variant ids are in chr[:,_,-]pos... format
end
% end

if isempty(chrsubset)
    error('call_regenie:unableToSubset', 'cannot assing annot_files to genotype files!\nfiles should be numbered properly!')
end

if isfield(opts, 'bed')
    bedfiles = opts.bed;
    [~, name, ext] = fileparts(opts.bed);
    genopatt = @(x)compose(strshare(name + ext, 'pattern', true), x);
    opts.bed = fullfile(opts.weshome, genopatt(chrsubset));

    if any(chrsubset == 23)
        bedX = bedfiles(contains(bedfiles, "_cX"));
        opts.bed(chrsubset == 23) = bedX;
    end

    if any(chrsubset == 24)
        bedX = bedfiles(contains(bedfiles, "_cY"));
        opts.bed(chrsubset == 24) = bedX;
    end
    
else % bgen. pgen should be implemented
    [~, name, ext] = fileparts(opts.bgen);
    genopatt = @(x)compose(strshare(name + ext, 'pattern', true), x);
    opts.bgen = fullfile(opts.weshome, genopatt(chrsubset));
    
    [~, name, ext] = fileparts(opts.sample);
    genopatt = @(x)compose(strshare(name + ext, 'pattern', true), x);
    opts.sample = fullfile(opts.weshome, genopatt(chrsubset));
end
end % END

%% ------------------------------------------------------------------------
function [cmd, extsetfiles] = getREGENIEcmd(opts, newopts, sz, parts)
if nargin < 4
    parts = [];
end
cmd = cell(sz, 1);
multiFields = ["anno_file", "set_list", "mask_def", "exclude",...
    "extract", "keep", "remove"]; % multi file options

% for set file chunk size to be run in parallel
if isfield(opts, 'anno_file')
    if isfield(opts, "annoChunk")
        chunknum = opts.annoChunk;
    else
        chunknum = round(opts.threads/sz);
    end
end
extsetfiles = cell(sz, 1); % temp extract set files 

for i = 1:sz % loop over chromosomes
    if isfield(opts, 'bgen')
        newopts.bgen = opts.bgen(i);
        newopts.sample = opts.sample(i);
    elseif isfield(opts, 'pgen')
        newopts.pgen = opts.pgen(i);
    else
        newopts.bed = opts.bed(i);
    end
    
    % per each chromosome
    for k = 1:numel(multiFields)
        if isfield(opts, multiFields(k)) && numel(opts.(multiFields(k))) > 1
            try
                newopts.(multiFields(k)) = opts.(multiFields(k))(i);
            catch
                fprintf("warning: field %s has fewer values than geno files, ignoring for this %s!\n", multiFields(k), opts.bgen(i))
                newopts = rmfield(newopts, multiFields(k));
            end
        end
    end
    
    % create chunk of set files for this chromosome 
    if isfield(opts, 'anno_file') && chunknum > 1
        newopts = createSetChunks(newopts, chunknum);
        extsetfiles{i} = newopts.extract_sets;
        for j = 1:numel(newopts.extract_sets)
            chunknewopts = newopts;
            chunknewopts.extract_sets = chunknewopts.extract_sets(j);
            if isempty(parts)
                chunknewopts.out = opts.out + "." + i + ".chunk" + j;
            else
                chunknewopts.out = opts.out + "." + i + parts + "chunk" + j;
            end
            
            %@16APR2023: extract_setlist and exclude_setlist are taken care
            %of in createSetChunks func. Should be removed because regenie
            %doesn't allow these with 'extract_sets' at the same time.
            rmfis = ["extract_setlist", "exclude_setlist"];
            for m = 1:numel(rmfis)
                if isfield(chunknewopts, rmfis(m))
                    chunknewopts = rmfield(chunknewopts, rmfis(m));
                end
            end
            loopopts = namedargs2cell(chunknewopts);
            cmd{i}(j, 1) = regenie(loopopts{:});
        end
    else
        newopts.out = opts.out + "." + i;
        loopopts = namedargs2cell(newopts); 
        cmd{i} = regenie(loopopts{:});
    end
    
end

extsetfiles(cellfun(@isempty, extsetfiles)) = [];
if isempty(extsetfiles); extsetfiles = []; end

end % END