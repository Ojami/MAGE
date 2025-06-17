function gwasCaller(opts)
% a wrapper for gwasrunner for gene-based/phenome-wide analysis.


% @07DEC2022: some bugs were fixed in gene_based_tests subfunction.
% @13MARCH2023: 'pheno_snps' can be a table with columns (only applicable
%               to UKBB data): file, bmi_adj, and mri. For all 4
%               combinations of bmi_adj and mri, different adjustments will
%               be applied. For instance, in case of bmi_adj and mri, are
%               false and true respectively for a trait, covariates to be
%               used are age at mri and sex.
% @26APR2023: 'mergeFolders' and 'mergeOutdir' have been added.
%             'mergeFolder' accepts a set of input folders each contains
%              summary stats from different calls to the function (e.g.
%              different adjustments, studies, etc.). Then
%              mergeDifferentGSSFolders subfunction merges all these
%              folders and writes them to 'mergeOutdir', which by default
%              is pwd/genes. Note that when this argument is used, the
%              function doesn't run anything and just merges the summary
%              stat files from previous runs.
% @22MAY2023: 'gtf37' and 'gtf38' options were added to replace Ensembl API
%              for finding gene positions (recommended). Previously, input
%              gene(s) symbols must have been the same between GRCh37/38 ,
%              and they were mapped later to GRCh37 for common variant
%              analysis. When using 'gtf37' and 'gtf38', input genes can be
%              a mix of GRCh38/37.
%              - a bug was fixed in mergeDifferentGSSFolders subfunction.
%              - 'dir' option was added: output dir for the results.
% @18JUNE2024: 'varset' option was added to include a set of desired rare
%              variant ids in the rare variant analysis only.
% 
% @03SEPT2024: 'qcinc' can be a cell of different sizes, with the same size
%               of pheno_genes or pred. This is useful, for instance when
%               one wants to analyze the same trait in different strata of
%               a covariate.
% @30SEP2024: 'interaction' argument was added (see gwasrunner). For the
%              moment, only applies to gene_based_tests subfunction.
% @26NOV2024: 'catCovarCell' option was added which accepts the same format
%              covarFile and phenoFile (i.e. cell of string vectors). This
%              is useful when categorical variables are different per each
%              pheno/covarfiles. This option is used only when
%              covarFile/phenoFile are used and in that case, overrides,
%              'catCovar' option. If there are multiple categorical
%              covariates for each pheno, they should be "," separated.

arguments
    opts.workers (1,1) double = 40
    
    opts.snps {mustBeA(opts.snps, "struct")} % list of snps to be analyzed (gwasrunner format).
    opts.genes {mustBeText, mustBeVector} % list of genes to be analyzed
    opts.varset {mustBeVector, mustBeText} % list of rare variants to be included in the analysis of rare variants

    % either a file for pheWas containing the list of traits, or a cell
    % array with each element corresponding to traits to be analyzed for
    % each gene in 'genes'
    opts.pheno_genes {mustBeVector} 
    opts.pheno_snps % same as above, but should be either a file, string vector, or table for 'snps' option. 
    opts.qc (1,1) double = 1
    opts.covar {mustBeText, mustBeVector} = ["Age", "Sex", "BMI"]; % list of covariates, it automatically adds PC1-10 and array batch.
    opts.qcinc {mustBeVector} % custom samples to be included
    opts.fitter {mustBeMember(opts.fitter, ["glm", "scoretest"])} = "scoretest"
    opts.win (1,1) double = 5e4 % region around each gene in kbp (only if 'genes' option is used)
    opts.bgenhome {mustBeFolder}  % sex.chr subfolder should contain sex chromosomes
    opts.weshome {mustBeFolder} 
    opts.maskdir {mustBeFolder}  % geneWrapper output
    opts.ignore_gene_based (1,1) logical = false % to ignore gene-based tests (REGENIE)
    opts.ignore_gene_scan (1,1) logical = false % to ignore gene scanning for common variant associations
    opts.ignore_regional_plots (1,1) logical = false % to ignore regional plots. This needs files to be present already in the target directories.
    opts.gene_scan_phewas (1,1) logical = true % perform PheWas on top hit per each trait from gene-scan analysis. 'pheno_snps' must be used in this case
    opts.report (1,1) logical = true % compile results into a report
    opts.report_only (1,1) logical = false % generate only the report and skip all other steps
    opts.report_name {mustBeTextScalar} = "report"
    opts.report_intro {mustBeText}
    opts.gtf37 {mustBeFile} = fullfile(fileparts(which("gwasCaller.m")), "GTF", "Homo_sapiens.GRCh37.87.gtf.gene.mat")
    opts.gtf38 {mustBeFile} = fullfile(fileparts(which("gwasCaller.m")), "GTF", "Homo_sapiens.GRCh38.107.gtf.gene.mat")
    opts.dir = pwd % creates a folder if it doesn't exist
    opts.infoscore (1,1) double = 0.7

    % merging different folders
    opts.mergeFolders {mustBeFolder}
    opts.mergeOutdir {mustBeTextScalar}

    % REGENIE inputs
    opts.pred {mustBeVector} % if left empty, typical linear/logisitc reg is performed
    opts.aaf_bins double {mustBeNonNan, mustBeNonnegative, mustBeVector} = [0.01, 0.005] % comma-separated list of AAF upper bounds to use when building masks [default is a single cutoff of 1%]
    opts.build_mask {mustBeTextScalar, mustBeMember(opts.build_mask, ["max", "sum", "comphet"])} = "max" % build masks using the maximum number of ALT alleles across sites ('max'; the default), or the sum of ALT alleles ('sum'), or thresholding the sum to 2 ('comphet')
    opts.skat_params (1,2) double {mustBeNonnegative, mustBeNonNan} = [1,25] % skat parameters
    opts.vc_maxAAF (1,1) double {mustBeInRange(opts.vc_maxAAF, 1e-12, 1)} = 0.01 % AAF upper bound to use for SKAT/ACAT-type tests [default is 100%]
    opts.phenoFile {mustBeVector} % the same size as covarFile and pred
    opts.covarFile {mustBeVector} % the same size as phenoFile and pred
    opts.condition {mustBeTextScalar}
    opts.write_mask (1,1) logical = false % write mask to PLINK bed format (does not work when building masks with 'sum')

    %@30SEP2024
    opts.interaction {mustBeTextScalar} = "" % adds SNP*trait interaction to the regression model. It's the name of MAT file and not the tag string.
    opts.catCovar {mustBeText, mustBeVector} % list of categorical covariates

    %@26NOV2024
    opts.catCovarCell {mustBeVector} % the same size as phenoFile and pred. overrides 'catCovar'

    %@20MARCH2025
    opts.debug (1,1) logical = false % for gwasrunner function to keep intermediate files
    
end

if isfield(opts, "mergeFolders")
    if isfield(opts, "mergeOutdir")
        mergeDifferentGSSFolders(opts.mergeFolders, opts.genes, outdir=opts.mergeOutdir)
    else
        mergeDifferentGSSFolders(opts.mergeFolders, opts.genes)
    end
    return
end

if isempty(gcp('nocreate')) && ~opts.report_only
    parpool("Processes", opts.workers);
    % parpool("Threads");
end

if opts.report_only, opts.report = true; end

if ~isfolder(opts.dir), mkdir(opts.dir); end

% get absolute path
[~, opts.dir] = fileattrib(opts.dir);
opts.dir = opts.dir.Name;

fixedFields = ["aaf_bins", "fitter", "build_mask", "skat_params", ...
    "vc_maxAAF", "pred", "qc", "qcinc", "covar", "condition", "write_mask"];

% check PheWide
if isfield(opts, 'snps')
    snps = opts.snps;
    pptime = tic;
    fid(1:numel(snps.snp), 1) = parallel.FevalFuture;
    for i = 1:numel(snps.snp)
        topts = opts;
        tmp.snp = snps.snp(i); tmp.chr = snps.chr(i);
        if isfield(snps, "tag"), tmp.tag = snps.tag(i); end
        topts.snps = tmp;
        % runPheWas(topts, fixedFields); %#DEBUG
        fid(i, 1) = parfeval(@runPheWas, 0, topts, fixedFields);
    end
    wait(fid)
    toc(pptime)    
end

% define genes and pertinent traits
if isfield(opts, 'genes')
    genes = opts.genes;
    pheno = opts.pheno_genes;

    if numel(pheno) ~= numel(genes)
        error("gwasCaller:wrongInputs", "pheno and genes must have the same lenght!")
    end
else
    return % nothing left to do!
end

qc = opts.qc; % unrelated Europeans (maximal unrelatedness)
win = opts.win; % region around each gene in kbp.

% load GTF files
if isfield(opts, "gtf37")
    opts.gtf37 = struct2table(load(opts.gtf37, "gene_name", "gene_id", "seqname", "start", "stop"));
    opts.gtf37 = renamevars(opts.gtf37, ["seqname", "stop"], ["chr", "end"]);
end

if isfield(opts, "gtf38")
    opts.gtf38 = struct2table(load(opts.gtf38, "gene_name", "gene_id", "seqname", "start", "stop"));
    opts.gtf38 = renamevars(opts.gtf38, ["seqname", "stop"], ["chr", "end"]);
end

for i = 1:numel(genes)
    if opts.report_only, continue; end
    if ~isfolder(fullfile(opts.dir, genes(i))), mkdir(fullfile(opts.dir, genes(i))); end

    fprintf('(%d of %d)-analysing gene %s\n', i, numel(genes), genes(i))

    if ~opts.ignore_gene_based || ~opts.ignore_gene_scan
        if isfield(opts, "gtf37")
            g = opts.gtf37(opts.gtf37.gene_name == genes(i), :);
            if isempty(g)
                g = opts.gtf38.gene_id(opts.gtf38.gene_name == genes(i));
                g = opts.gtf37(opts.gtf37.gene_id == g, :);
            else
                % match gene name to GRCh38
                try
                    genes(i) = opts.gtf38.gene_name(opts.gtf38.gene_id == g.gene_id);
                catch % different genes IDs between 37 and 38, keep the gene name
                end
            end
        else
            g = EnsemblREST(genes(i), "lookup", refGenome="37", geneSymbol=true);
        end
        s.snp = (g.start-win) + "-" + (g.end+win);

        if opts.bgenhome.contains("Imputed_TOPMED") % GRCh38
            g = opts.gtf38(opts.gtf38.gene_name == genes(i), :);
            s.snp = (g.start-win) + "-" + (g.end+win);
        end
        s.chr = string(g.chr);
    end
    
    if ~opts.ignore_gene_based
        % gene-based tests
        if isnan(double(string(g.chr)))
            weshome = fullfile(opts.weshome, "sex.chrom");
        else
            weshome = opts.weshome;
        end

        if ~isfile(fullfile(opts.dir, genes(i), genes(i)+".RVA.xlsx"))

            gopts = opts;
            rfis = setdiff(fieldnames(gopts), union(fixedFields, ...
                ["phenoFile", "covarFile", "maskdir", "interaction", ...
                "catCovar", "varset", "catCovarCell", "debug"]));
            gopts = rmfield(gopts, rfis);
            gopts.weshome = weshome; 
            gopts.output = fullfile(opts.dir, genes(i), genes(i)+".RVA");
            if isfield(opts, 'pred')
                gopts.pred = string(gopts.pred{i}); 
                gopts.qc = 0; 
            end

            if isfield(opts, 'varset')
                gopts.varset = opts.varset;
            end

            if isfield(opts, "phenoFile"), gopts.phenoFile = opts.phenoFile{i}; end
            if isfield(opts, "covarFile"), gopts.covarFile = opts.covarFile{i}; end

            %@26NOV2024: see if 'catCovarCell' argument is used along with
            %'phenoFile'
            if isfield(opts, "phenoFile") && isfield(opts, "catCovarCell")
                gopts.catCovarCell = opts.catCovarCell{i};
                assert(numel(gopts.catCovarCell) == numel(gopts.phenoFile), "catCovarCell and phenoFile size mismatch!")
                if isfield(gopts, "catCovar")
                    gopts = rmfield(gopts, "catCovar");
                end
            end

            %@03SEPT2024: 'qcinc' can be a cell array
            if isfield(opts, "qcinc") && iscell(opts.qcinc)
                gopts.qcinc = opts.qcinc{i};
            end

            gopts = namedargs2cell(gopts);
            gene_based_tests(thisgene, string(pheno{i}), gopts{:})
        end
    end

    if ~opts.ignore_gene_scan
        % common variant associations
        bgenhome = opts.bgenhome;
        infoscore = opts.infoscore;

        if ~isfile(fullfile(opts.dir, genes(i), genes(i)+".SVA.xlsx"))
            sopts = opts;
            rfis = setdiff(fieldnames(sopts), fixedFields);
            sopts = rmfield(sopts, rfis);
            sopts.output = fullfile(opts.dir, genes(i), genes(i)+".SVA");
            sopts.bgenhome = bgenhome;
            sopts.infoscore = infoscore;
            sopts.qc = qc;
            sopts.parallel=true;
            sopts.trait = string(pheno{i});
            if isfield(opts, 'pred')
                sopts.pred = string(opts.pred{i}); 
                sopts.qc = 0; 
                sopts.regenie = true;
            end

            if isfield(opts, "phenoFile"), sopts.phenoFile = opts.phenoFile{i}; end
            if isfield(opts, "covarFile"), sopts.covarFile = opts.covarFile{i}; end

            %@26NOV2024: see if 'catCovarCell' argument is used along with
            %'phenoFile'
            if isfield(opts, "phenoFile") && isfield(opts, "catCovarCell")
                sopts.catCovarCell = opts.catCovarCell{i};
                assert(numel(sopts.catCovarCell) == numel(sopts.phenoFile), "catCovarCell and phenoFile size mismatch!")
                if isfield(sopts, "catCovar")
                    sopts = rmfield(sopts, "catCovar");
                end
            end

            %@03SEPT2024: 'qcinc' can be a cell array
            if isfield(opts, "qcinc") && iscell(opts.qcinc)
                sopts.qcinc = opts.qcinc{i};
            end
            
            sopts = namedargs2cell(sopts);
            gwasrunner(s, sopts{:})

            % remove NaN columns
            gss = readtable(fullfile(opts.dir, genes(i), genes(i)+".SVA.xlsx"), ...
                "TextType", "string", "VariableNamingRule", "preserve");
            gss = rmmissing(gss, 2, "MinNumMissing", height(gss));
            writetable(gss, ...
                fullfile(opts.dir, genes(i), genes(i)+".SVA.xlsx"), ...
                "WriteMode", "replacefile");
        end

        if opts.gene_scan_phewas % do PheWas for top hits 
            gss = readtable(fullfile(opts.dir, genes(i), genes(i)+".SVA.xlsx"), ...
                "VariableNamingRule","preserve","TextType","string");
            gss = fillmissing(gss, "previous", "DataVariables", "Pheno");

            % only variants passing FDR BH cut-off: q <= 0.05
            gss = groupfilter(gss, 'Pheno', @(x) mafdr(x, "BH", true) <= 0.05, "P-adj");
            if ~isempty(gss)
                gss = groupfilter(gss, 'Pheno', @(x) x == min(x), "P-adj");
                % only 1 variant per each trait
                sumcheck = groupsummary(gss, "Pheno");
                sumcheck = sumcheck.Pheno(sumcheck.(2) > 1);
                if ~isempty(sumcheck)
                    gssTmp1 = gss;
                    gsstmp2 = cell(numel(sumcheck), 1);
                    for k = 1:numel(sumcheck)
                        idx = gss.Pheno == sumcheck(k);
                        gtmp = gss(idx, :);
                        [~, idx] = min(gtmp.("P-adj"));
                        gsstmp2{k, 1} = gtmp(idx, :);
                    end
                    gssTmp1(ismember(gssTmp1.Pheno, sumcheck), :) = [];
                    gss = [gssTmp1; vertcat(gsstmp2{:})];
                end

                for k = 1:height(gss)                
                    snp.snp = gss.SNP(k); snp.chr = string(gss.CHR(k));
                    snp.tag = genes(i) + "." + matlab.lang.makeValidName(snp.snp);
                    ropts = opts;
                    ropts.snps = snp;
                    ropts.bgenhome = bgenhome;
                    ropts.dir = fullfile(opts.dir, genes(i));
                    ropts.output = snp.tag + "." + matlab.lang.makeValidName(gss.Pheno(k));
                    runPheWas(ropts, fixedFields)
                end
            end
        end
    end % ignore_gene_scan

    if ~opts.ignore_regional_plots
        % draw regional plots
        call_genePlotter(fullfile(opts.dir, genes(i), genes(i)+".SVA.xlsx"), opts.bgenhome)
    end
end

% make report
if opts.report
    ropts = opts;
    fis = setdiff(fieldnames(ropts), ["snps", "report_name", "report_intro", "dir"]);
    ropts = rmfield(ropts, fis);
    ropts.genes = genes; ropts.traits = pheno;
    if isfield(ropts, "report_intro")
        ropts.intro = ropts.report_intro;
        ropts = rmfield(ropts, "report_intro");
    end
    ropts = namedargs2cell(ropts);
    makeReport(ropts{:})
end

end % END

%% subfunctions ===========================================================
function runPheWas(opts, fixedFields)

% gwasrunner options
if nargin > 1
    ropts = opts;
    fixedFields(fixedFields == "pred") = []; % reduce computational time
    rfis = setdiff(fieldnames(ropts), fixedFields);
    ropts = rmfield(ropts, rfis);
end

if ~istable(opts.pheno_snps) && isfile(opts.pheno_snps)
    pheno = load(opts.pheno_snps);
    fis = string(fieldnames(pheno));
    pheno = pheno.(fis(1));
else
    pheno = opts.pheno_snps;
end

if ~isfield(opts, 'dir'), opts.dir = pwd; end

% TODO: check for chr X
output = "PheWide";
snps = opts.snps;

if isfield(opts, "output")
    output = output + "." + opts.output;
elseif isfield(snps, "tag")
    output = output + "." + matlab.lang.makeValidName(snps.tag(1));
else
    output = output + "." + matlab.lang.makeValidName(snps.snp(1));
end

output = fullfile(opts.dir, output);
if istable(pheno)
    sc = ["bmiY_mriY", "bmiY_mriN", "bmiN_mriY", "bmiN_mriN"];
    fid(1:numel(sc), 1) = parallel.FevalFuture;
    for k = 1:numel(sc)
        thisc = split(sc(k), "_");

        curr_covars = opts.covar; % inout covariates
        curr_covars(startsWith(curr_covars.lower, ["age", "bmi"])) = []; % remove age and bmi
    
        if thisc(2).endsWith("Y") % MRI
            covars = ["AgeMRI", "BMI_2", "Sex"];
            idx = pheno.mri;
        else
            idx = ~pheno.mri;
            covars = ["Age", "BMI", "Sex"];
        end
    
        if thisc(1).endsWith("Y") 
            idx = idx & pheno.bmi_adj;
        else % don't adjust for BMI
            covars(ismember(covars, ["BMI_2", "BMI"])) = [];
            idx = idx & ~pheno.bmi_adj;
        end

        covars = union(covars, curr_covars); 
        
        if ~any(idx), continue; end
        
        checkfiles = getfilenames(opts.dir, "xlsx").xlsx;
        if ~isempty(checkfiles) && any(checkfiles.endsWith(sc(k) + ".xlsx"))
            continue
        end

        topts = ropts;
        topts.trait = pheno.file(idx);
        topts.output = output + "." + sc(k);
        topts.covar = covars;
        topts.gpu = false;
        topts.fitter = opts.fitter;
        topts.bgenhome = opts.bgenhome;
        topts.infoscore = opts.infoscore;
        topts = namedargs2cell(topts);
        fid(k, 1) = parfeval(@gwasrunner, 0, snps, topts{:});

    end

    if ~all(ismember({fid.State}, 'unavailable'))
        wait(fid)
    end
    
    % merge files
    checkfiles = getfilenames(opts.dir, "xlsx", "fullpath", true).xlsx;
    filepref = output + "." + sc;
    files = checkfiles(checkfiles.startsWith(filepref));
    df = fileDatastore(files, ...
        "ReadFcn", @(x)readtable(x, TextType="string", VariableNamingRule="preserve"));
    df = readall(df);
    for j = 1:numel(df)
        idx = find(ismember(colnames(df{j}), ["β|OR", "β", "OR"]));
        df{j}.Properties.VariableNames(idx) = "β|OR";
    end
    df = vertcat(df{:});
    df.P_Bonferroni = min(1, df.("P-adj").*height(df));
    df.P_FDR = mafdr(df.("P-adj"), "BH", true);
    df = rmmissing(df, 2, "MinNumMissing", height(df));
    df = sortrows(df, "P_FDR", "ascend");
    writetable(df, output + ".xlsx", "WriteMode", "replacefile");
    delete(files{:});

else
    if ~isfile(output + ".xlsx")
        topts = ropts;
        topts = namedargs2cell(topts);
        topts.output = output;
        topts.gpu =false;
        topts.parallel=false;
        gwasrunner(snps, topts{:})
    end

    % adjust for BH and Bonferroni
    gss = readtable(output + ".xlsx", "TextType", "string", "VariableNamingRule", "preserve");
    gss.P_Bonferroni = min(1, gss.("P-adj").*height(gss));
    gss.P_FDR = mafdr(gss.("P-adj"), "BH", true);
    gss = rmmissing(gss, 2, "MinNumMissing", height(gss));
    gss = sortrows(gss, "P-adj");
    writetable(gss, output + ".xlsx", "WriteMode", "replacefile");
end

end

%% ------------------------------------------------------------------------
function call_genePlotter(gss, bgenhome)

[out, gss, ext] = fileparts(gss);
if isempty(out) || string(out) == ""
    out = pwd;
end

if ~isfolder(fullfile(out, "regional"))
    mkdir(fullfile(out, "regional"))
end

gss = gss + ext;
gene = extractBefore(gss, ".");
gss = readtable(fullfile(out, gss),"VariableNamingRule","preserve","TextType","string");
gss = fillmissing(gss, "previous", "DataVariables", "Pheno");
gss(:, "95% CI") = [];
try gss.P = []; catch, end
try gss(:, ["N.case", "AF.case", "AF.control"]) = []; catch, end
betacol = find(ismember(gss.Properties.VariableNames, ["OR", "OR|β", "β|OR", "β"]));
betacol = string(gss.Properties.VariableNames(betacol));
gss = renamevars(gss, [betacol, "P-adj"], ["beta", "P"]);

gss.Pheno = matlab.lang.makeValidName(gss.Pheno);
pheno = unique(gss.Pheno);
for i = 1:numel(pheno)
    tab.(pheno(i)) = gss(gss.Pheno == pheno(i), :);
end

if any(string(gss.CHR(1)).lower == ["x", "23", "24", "y"])
    bgenhome = fullfile(bgenhome, "sex.chrom");
end

fis = string(fieldnames(tab));
for i = 1:numel(fis)
    genePlotter(tab.(fis(i)), 'ldmethod', "insample", ...
        "backup", true, 'parallel', true, ...
        'title', replace(fis(i), "_", "-") + " (" + gene + ")", ...
        'outdir', fullfile(out, "regional"), ...
        'save', true, 'savefig', true, 'resolution', 300, ...
        'significance', 'fdr', 'bgenhome', bgenhome);
end

end % END

%% ------------------------------------------------------------------------
function gene_based_tests(gene, traits, opts)

arguments
    gene
    traits {mustBeText}
    opts.qc (1,1) double
    opts.varset {mustBeVector, mustBeText} % list of rare variants to be included in the analysis of rare variants
    opts.qcinc {mustBeVector}
    opts.covar {mustBeText, mustBeVector} = ["Age", "Sex", "BMI"]; % list of covariates, it automatically adds PC1-10 and array batch.
    opts.weshome {mustBeFolder}
    opts.maskdir {mustBeFolder}
    opts.output {mustBeTextScalar} = ""
    opts.pred {mustBeVector, mustBeText} % if left empty, typical linear/logisitc reg is performed
    opts.aaf_bins double {mustBeNonNan, mustBeNonnegative, mustBeVector} = [0.01, 0.005] % comma-separated list of AAF upper bounds to use when building masks [default is a single cutoff of 1%]
    opts.build_mask {mustBeTextScalar, mustBeMember(opts.build_mask, ["max", "sum", "comphet"])} = "max" % build masks using the maximum number of ALT alleles across sites ('max'; the default), or the sum of ALT alleles ('sum'), or thresholding the sum to 2 ('comphet')
    opts.skat_params (1,2) double {mustBeNonnegative, mustBeNonNan} = [1,25] % skat parameters
    opts.vc_maxAAF (1,1) double {mustBeInRange(opts.vc_maxAAF, 1e-12, 1)} = 0.01 % AAF upper bound to use for SKAT/ACAT-type tests [default is 100%]
    opts.fitter {mustBeMember(opts.fitter, ["glm", "scoretest"])} = "scoretest"
    opts.phenoFile {mustBeFile} % the same size as covarFile and pred
    opts.covarFile {mustBeFile} % the same size as phenoFile and pred
    opts.condition {mustBeTextScalar}
    opts.interaction {mustBeTextScalar} = "" % adds SNP*trait interaction to the regression model. It's the name of MAT file and not the tag string.
    opts.catCovar {mustBeText, mustBeVector} % list of categorical covariates
    opts.catCovarCell {mustBeText, mustBeVector}
    opts.debug (1,1) logical = false
    opts.write_mask (1,1) logical = false % write mask to PLINK bed format (does not work when building masks with 'sum')
end

gene.gene = string(gene.gene);
if isfield(gene, "chr"), gene.chr = string(gene.chr); end

if string(opts.output) == ""
    opts.output = gene.gene;
    opts.dir = pwd;
else
    [opts.dir, name, ext] = fileparts(opts.output);
    opts.output = name + ext;
    if opts.dir == "" || isempty(opts.dir)
        opts.dir = pwd;
    end
    if ~isfolder(opts.dir)
        mkdir(opts.dir)
    end
end

% split traits into chunks to be run in parallel
%@03SEPT2024: qcinc can be a cell array
if isfield(opts, "qcinc") && iscell(opts.qcinc)
    idx = 1:(numel(traits)+1);
else
    idx = unique([1:3:numel(traits)+1, numel(traits)+1]);
    if numel(idx) < 3
        idx = unique(round(linspace(1, numel(traits)+1, 4)), "stable");
    end
end

% find anno/set/mask files
if isfield(opts, "maskdir")
    files = getfilenames(opts.maskdir, "txt").txt;
    ff.afiles = files(contains(files, ".annotation."));
    ff.sfiles = files(files.contains(".set."));
    ff.mfiles = files(files.contains(".mask."));
    ff.apatt = replace(strshare(ff.afiles, "pattern", true), "%d", "%s");
    ff.spatt = replace(strshare(ff.sfiles, "pattern", true), "%d", "%s");
    ff.mpatt = replace(strshare(ff.mfiles, "pattern", true), "%d", "%s");
    
    check.fi1 = ["anno_file", "set_list", "mask_def"];
    check.fi2 = ["a", "s", "m"] + "patt";
    check.fi3 = ["a", "s", "m"]  + "files";
    for i = 1:numel(check.fi1)

        opts.(check.fi1(i)) = compose(ff.(check.fi2(i)), gene.chr);
        if ~any(opts.(check.fi1(i)) == ff.(check.fi3(i))) % X/Y chromosomes?
            opts.(check.fi1(i)) = compose(ff.(check.fi2(i)), "X");
            if ~any(opts.(check.fi1(i)) == ff.(check.fi3(i)))
                opts.(check.fi1(i)) = compose(ff.(check.fi2(i)), "23");
                if ~any(opts.(check.fi1(i)) == ff.(check.fi3(i)))
                    opts.(check.fi1(i)) = compose(ff.(check.fi2(i)), "Y");
                    if ~any(opts.(check.fi1(i)) == ff.(check.fi3(i)))
                        opts.(check.fi1(i)) = compose(ff.(check.fi2(i)), "24");
                    end
                end
            end
        end

    end

    opts.anno_file = fullfile(opts.maskdir, opts.anno_file);
    opts.set_list = fullfile(opts.maskdir, opts.set_list);
    opts.mask_def = fullfile(opts.maskdir, opts.mask_def);
    opts = rmfield(opts, "maskdir");
end

tmp = getRandomName("g", 5); 
g.gene = gene.gene;

if isfield(opts, "varset")
    g.variants = opts.varset;
    opts = rmfield(opts, "varset");
end

% attach necessary files (note that gene_based_tests relies on parallel
% toolbox)
hd = fileparts(which('bfilereader.m'));
file = fullfile(hd, 'bFileReaderDep.class');
pobj = gcp("nocreate");
if isempty(pobj.AttachedFiles)
    addAttachedFiles(gcp, file);
end

% % decide on number of threads for each REGENIE job
pptime = tic;
fid(1:numel(idx)-1, 1) = parallel.FevalFuture;
for i = 1:numel(idx)-1
    if ~isfolder(fullfile(opts.dir, tmp + i)), mkdir(fullfile(opts.dir, tmp + i)); end
    ropts = opts;
    ropts.output = fullfile(opts.dir, tmp + i, tmp + i);
    ropts.trait = traits(idx(i):idx(i+1)-1);
    ropts.firth = true;
    ropts.regenie = true;
    if isfield(ropts, 'pred')
        ropts.pred = ropts.pred(idx(i):idx(i+1)-1);
    end

    if isfield(ropts, "phenoFile")
        ropts.phenoFile = ropts.phenoFile(idx(i):idx(i+1)-1);

        %@26NOV2024
        if isfield(ropts, "catCovarCell")
            ropts.catCovarCell = ropts.catCovarCell(idx(i):idx(i+1)-1);
        end
    end

    if isfield(ropts, "covarFile")
        ropts.covarFile = ropts.covarFile(idx(i):idx(i+1)-1);
    end

    %@03SEPT2024: qcinc can be a cell array
    if isfield(opts, "qcinc") && iscell(opts.qcinc)
        ropts.qcinc = ropts.qcinc{idx(i):idx(i+1)-1};
    end
    ropts.gpu = false;

    ropts = rmfield(ropts, "dir");
    ropts.phenopath = fullfile(fileparts(which("gwasrunner.m")), "UKBB.Toolbox", "UKB_PHENO");
    ropts = namedargs2cell(ropts);

    % gwasrunner(g, ropts{:}) % DEBUG
    fid(i, 1) = parfeval(@gwasrunner, 0, g, ropts{:});
end

if ~all(ismember({fid.State}, 'unavailable'))
    wait(fid)
end
toc(pptime)
delete(fid)

% merge output files
% first check what sheets overlap and create a map struct
data_key = cell(numel(idx)-1, 1);
for i = 1:numel(idx)-1
    file = fullfile(opts.dir, tmp + i, tmp + i + ".xlsx");
    if isfile(file)
        data_key{i, 1} = array2table([sheetnames(file), ...
            matlab.lang.makeValidName(sheetnames(file))], ...
            VariableNames=["s", "k"]);
    end
end
data_key = unique(vertcat(data_key{:}), 'rows', 'stable');

if isempty(data_key)
    disp("gene was not found or something else went wrong!")
    for i = 1:numel(idx)-1
        rmdir(fullfile(opts.dir, tmp + i), "s")
    end
    return
end

data = struct;
for i = 1:numel(idx)-1
    file = fullfile(opts.dir, tmp + i, tmp + i + ".xlsx");
    snames = sheetnames(file);
    [~, f2] = ismember(snames, data_key.s); f2(f2 < 1) = [];
    key = data_key(f2, :);
    legend_sheet_idx = key.s.lower.startsWith("legend");
    
    
    if any(legend_sheet_idx)
        data.(key.k(legend_sheet_idx)){i, 1} = readmatrix(file, ...
            Sheet=key.s(legend_sheet_idx), ...
            OutputType="string");
    end
    key(legend_sheet_idx, :) =  [];
    for j = 1:height(key)
        data.(key.k(j)){i, 1} = readtable(file, VariableNamingRule="preserve", Sheet=key.s(j));
    end    
    
    %@20MARCH2025: keep intermediate files
    if opts.debug
        junk_files = getfilenames(fullfile(opts.dir, tmp + i)).x_;
        junk_files(junk_files.endsWith([".xls", ".xlsx"])) = []; % already been read
        jpath1 = fullfile(opts.dir, tmp + i, junk_files);
        jpath2 = fullfile(opts.dir, junk_files);
        arrayfun(@(x,y)movefile(x,y),jpath1, jpath2);
    end
    rmdir(fullfile(opts.dir, tmp + i), "s")
    
end

fi = string(fieldnames(data));
for j = 1:numel(fi)
    if istable(data.(fi(j))), continue; end
    try
        data.(fi(j)) = vertcat(data.(fi(j)){:});
    catch % tables of different size
        tmp = data.(fi(j));
        tmpjoin = outerjoin(tmp{1}, tmp{2}, MergeKeys=true);
        for k = 3:numel(tmp)
            tmpjoin = outerjoin(tmpjoin, tmp{k}, MergeKeys=true);
        end
        data.(fi(j)) = tmpjoin;
    end
end

overall_sheet_idx = data_key.s.lower.startsWith("overall.info");
legend_sheet_idx = data_key.s.lower.startsWith("legend");
data_key = [data_key(overall_sheet_idx, :);
    data_key(~(overall_sheet_idx|legend_sheet_idx),:);
    data_key(legend_sheet_idx, :)];
legend_sheet_idx = data_key.s.lower.startsWith("legend");
lsheet = data_key.k(end);
data.(lsheet) = unique(data.(lsheet), "rows", "stable"); % Legend table
data.(lsheet) = array2table(data.(lsheet));
for i = 1:height(data_key)
    if legend_sheet_idx(i), writeFlag = false; else, writeFlag = true; end
    writetable(data.(data_key.k(i)), ...
        fullfile(opts.dir, opts.output + ".xlsx"), Sheet=data_key.s(i), ...
        WriteVariableNames=writeFlag, AutoFitWidth=true);
end


end % END

%% ------------------------------------------------------------------------
function makeReport(opts)

arguments 
    opts.snps {mustBeA(opts.snps, "struct")}
    opts.genes {mustBeText, mustBeVector}
    opts.traits {mustBeVector}
    opts.intro {mustBeVector} % introduction
    opts.theme {mustBeMember(opts.theme, ["flatly", "litera", "journal"])} = "journal"
    opts.fontsize (1,1) double = 12
    opts.resolution (1,1) double = 150
    opts.report_name {mustBeTextScalar} = "report"
    opts.dir {mustBeFolder} = pwd
end

rep.header = ["---", "title: Summary statistics report",...
"author: NA", 'date: "`r format(Sys.time(), ''%d %B, %Y'')`"', ...
"output:", "  html_document:", "    toc: true", "    toc_float: true",...
"    number_sections: true", "    theme:", "      bootswatch: " + opts.theme,...
"---", "", "```{=html}", '<style type="text/css">', "  body{", ...
"  font-size: " + opts.fontsize + "pt;", "}", "</style>", "```",...
"```{r include = FALSE}", "knitr::opts_chunk$set(echo=F, comment = NA)", "```"]';

if isfield(opts, "intro") 
    rep.intro = opts.intro;
else
    % default example
    rep.intro = ["", "# Summary"]';
end

rep.sep = ["*****"; "<br/><br/>"; ""]; % section break

% results section (add traits info) ---------------------------------------
[g, t, d] = deal(cell(numel(opts.genes), 1));
for i = 1:numel(opts.genes)
    tmp = string(opts.traits{i});
    for j = 1:numel(tmp)
        trt = load(tmp(j)).UKB_STRUCT_ALL;
        t{i}(j) = string(trt.tag);
        if isfield(trt, "desc") || (isfield(trt, "info") && isfield(trt.info, "desc"))
            d{i}(j) = string(trt.info.desc(1));
         else
             d{i}(j) = "";
        end

    end
    if isrow(t{i}), t{i} = t{i}'; end
    g{i} = repmat("", numel(t{i}), 1);
    g{i}(1) = opts.genes(i);
end
g = vertcat(g{:})'; t = vertcat(t{:})'; d = horzcat(d{:});

% remove dates from tags
rt = t;
idx = contains(t, regexpPattern('(\w+) (\d+)|(\w+)(\d+)'));
t(idx) = strtrim(regexprep(t(idx), '\(.*\)', ''));

rep.res = ["# Results", "All summary statistics (tables and figures)" + ...
    " have been provided. Traits analyzed per each gene have been" + ...
    " summarized below:", "", "```{r, warning=F}", "library(kableExtra)",... 
    "", "```"]';
rep.res(6) = "kbl(data.frame(Gene=" + rVec(g) + ", Trait=" + rVec(t) + ...
    ", Note=" + rVec(d) + "), " + ...
    "caption = 'Table 1: Traits analyzed per each gene')%>%"+...
    " kable_paper('striped', full_width = F)";

% add gene blocks ---------------------------------------------------------
tab = table(g', t', rt', d', 'VariableNames', {'g', 't', 'rt', 'd'});
tab.g(tab.g == "") = missing;
tab = fillmissing(tab, "previous","DataVariables", "g");

tmpFigs = cell(numel(opts.genes), 1); % png images (lower resol.) only for report and to be removed in th end
rep.genes = cell(numel(opts.genes), 1);
for i = 1:numel(opts.genes)
    hm = fullfile(opts.dir, opts.genes(i));
    tmp = tab(tab.g == opts.genes(i), :);
    tmp.mat = matlab.lang.makeValidName(tmp.rt);
    
    if isfolder(fullfile(hm, "regional"))
        % regional plots 
        close gcf force
        figs = getfilenames(fullfile(hm, "regional"), "fig", "fullpath", true).fig;
        geneR = cell(numel(figs), 1);
        tmpFigs2 = strings(numel(figs), 1);
        for j = 1:numel(figs)
    
            [~, name] = fileparts(figs(j));
            name = erase(name, [lower(opts.genes(i)), opts.genes(i)]);
    
            % what trait is this figure for?
            idx = find(arrayfun(@(x)startsWith(name, x), tmp.mat));
            if numel(idx) ~= 1
                idx = find(arrayfun(@(x)ismember(replace(name, "__", ""), x), tmp.mat));
                if numel(idx) ~= 1
                    continue % error("cannot find the trait!")
                end
            end
            
            if ~isfile(fullfile(hm, "regional", name + ".png"))
                f = openfig(figs(j), "visible");
                f.Children(2).Title.String = "";  
                exportgraphics(f, fullfile(hm, "regional", name + ".png"));
                close(f)
            end
    
            name = fullfile(hm, "regional", name + ".png");
            tmpFigs2(j, 1) = name;
            name = erase(name, opts.dir);
            name = regexprep(name, '(^\\|^//)', '');
            name = replace(name, filesep, [filesep, filesep]);
            geneR{j}(1, 1) = '```{r, out.width="100%", fig.align = "center", fig.cap = "' + ...
                tmp.t(idx) + '",fig.link="' + name + '"}';
            geneR{j}(2, 1) = "knitr::include_graphics('" + name + "', dpi = " + opts.resolution + ")";
            geneR{j}(3, 1) = "```";
            geneR{j}(4, 1) = "<br/><br/>";
        end
        
        tmpFigs{i, 1} = tmpFigs2; % to be deleted
        rep.genes{i} = [""; "### Regional plots"; vertcat(geneR{:})];
    end

    % phenome-wide and rare-variant analysis
    xfiles = getfilenames(hm, "xlsx", "fullpath", true).xlsx;
    
    rva = xfiles(contains(xfiles, ".RVA."));
    phe = xfiles(contains(xfiles, "PheWide."));

    if ~isempty(phe)
        %@13NOV2024: only show unique snps
        all_phe_snps = extractBetween(phe, opts.genes(i) + ".", ".");
        [~, u_phe_snps_idx] = unique(all_phe_snps);
        phe = phe(u_phe_snps_idx);

        geneRphe = cell(numel(phe), 1);
        for j = 1:numel(phe)
            snptag = extractBetween(phe(j), opts.genes(i)+".", ".");
            snptag = opts.genes(i) + " " + snptag;
            geneRphe{j}(1,1) = "tab" + j + ' = readxl::read_xlsx("' + ...
                replace(phe(j), filesep, [filesep, filesep]) + '")';
            geneRphe{j}(2,1) = "tab" + j + "=tab" + j + "[tab" + j + "$`P-adj` <= 0.05, ]";
            geneRphe{j}(3,1) =  "if (dim(tab" + j + ")[1] > 0) {";
            geneRphe{j}(4,1) = "DT::datatable(tab" + j + ", extensions" + ...
                "= c('FixedColumns','FixedHeader'),options = list(scrollX" + ...
                "= T,paging = T,fixedHeader=T, fixedColumns = " + ...
                "list(leftColumns = 2), pageLength = 5, " + ...
                "autoWidth = TRUE), rownames = T, caption = " + ...
                "htmltools::tags$caption(style = 'caption-side: " + ...
                "bottom; text-align: center;', htmltools::em('Phenome-wide summary statistics for " + ...
                snptag + "')))";
            geneRphe{j}(5,1) = "} else {cat('No P-value < 0.05')}";
        end
        geneRphe = [""; "### Phenome-Wide"; ""; "```{r}"; vertcat(geneRphe{:}); "```"; ""];
        rep.genes{i} = [rep.genes{i}; geneRphe];
    end

    if ~isempty(rva)
        geneRphe = cell(numel(rva), 1);
        for j = 1:numel(rva)

            % MAF < 1%
            geneRphe{j}(1,1) = "tab" + j + ' = readxl::read_xlsx("' + ...
                replace(rva(j), filesep, [filesep, filesep]) + ...
                '", sheet = "0.01")';
            geneRphe{j}(2,1) = "idx = apply(tab" + j + "[, startsWith(colnames(tab" + j + "), 'ADD')] <= 0.05, 1, function(x) any(x, na.rm = TRUE))";
            geneRphe{j}(3,1) = "tab" + j + "=tab" + j + "[idx, ]";
            geneRphe{j}(4,1) = "add_cols = which(startsWith(colnames(tab" + j + "), 'ADD'))";
            geneRphe{j}(5,1) = "tab" + j + "[add_cols] <- lapply(tab" + j + "[add_cols], function(x) ifelse(is.na(x), x, formatC(x, format = 'E', digits = 2)))";
            geneRphe{j}(6,1) =  "if (dim(tab" + j + ")[1] > 0) {";
            geneRphe{j}(7,1) = "colnames(tab" + j + ") <- gsub('^ADD-', '', colnames(tab" + j + "))";
            geneRphe{j}(8,1) = "tab" + j + " <- tab" + j + "[, !colnames(tab" + j + ") %in% c('ID', 'CHR')]";
            geneRphe{j}(9,1) = "tab" + j + "$Mask <- gsub('^mask_|\\.0\\.01$', '', tab" + j + "$Mask)";

            % change formats of beta/OR and SE
            geneRphe{j}(10,1) = "idx = grep('β|OR\\(', colnames(tab" + j + "))";
            geneRphe{j}(11,1) = "beta_se_col = c('SE', colnames(tab" + j + ")[idx[1]])";
            geneRphe{j}(12,1) = "tab" + j + "[beta_se_col] <- " + ...
                "lapply(tab" + j + "[beta_se_col], function(x) " + ...
                "ifelse(!is.na(x) & abs(x) < 0.01, formatC(x, format = 'E', digits = 2), sprintf('%.2f', x)))";

            geneRphe{j}(13,1) = "DT::datatable(tab" + j + ", extensions" + ...
                "= c('FixedColumns','FixedHeader'),options = list(scrollX" + ...
                "= T,paging = T,fixedHeader=T, fixedColumns = " + ...
                "list(leftColumns = 2), pageLength = 5, " + ...
                "autoWidth = TRUE), rownames = T, caption = " + ...
                "htmltools::tags$caption(style = 'caption-side: " + ...
                "bottom; text-align: center;', htmltools::em('Gene-based summary statistics (MAF < 1%) for " + ...
                opts.genes(i) + "')))";
            geneRphe{j}(14,1) = "} else {cat('No P-value < 0.05 for MAF < 1%')}";

            % Joint test
            geneRphe{j}(17,1) = "tab_" + j + ' = readxl::read_xlsx("' + ...
                replace(rva(j), filesep, [filesep, filesep]) + ...
                '", sheet = "Joint")';
            geneRphe{j}(18,1) = "idx = apply(tab_" + j + "[, startsWith(colnames(tab_" + j + "), 'ADD')] <= 0.05, 1, function(x) any(x, na.rm = TRUE))";
            geneRphe{j}(19,1) = "tab_" + j + "=tab_" + j + "[idx, ]";
            geneRphe{j}(20,1) = "add_cols = which(startsWith(colnames(tab_" + j + "), 'ADD'))";
            geneRphe{j}(21,1) = "tab_" + j + "[add_cols] <- lapply(tab_" + j + "[add_cols], function(x) ifelse(is.na(x), x, formatC(x, format = 'E', digits = 2)))";
            geneRphe{j}(22,1) =  "if (dim(tab_" + j + ")[1] > 0) {";
            geneRphe{j}(23,1) = "colnames(tab_" + j + ") <- gsub('^ADD-', '', colnames(tab_" + j + "))";
            geneRphe{j}(24,1) = "tab_" + j + " <- tab_" + j + "[, !colnames(tab_" + j + ") %in% c('ID', 'CHR', 'Mask')]";

            geneRphe{j}(25,1) = "DT::datatable(tab_" + j + ", extensions" + ...
                "= c('FixedColumns','FixedHeader'),options = list(scrollX" + ...
                "= T,paging = T,fixedHeader=T, fixedColumns = " + ...
                "list(leftColumns = 2), pageLength = 5, " + ...
                "autoWidth = TRUE), rownames = T, caption = " + ...
                "htmltools::tags$caption(style = 'caption-side: " + ...
                "bottom; text-align: center;', htmltools::em('Gene-based summary statistics (joint) for " + ...
                opts.genes(i) + "')))";

            geneRphe{j}(26,1) = "} else {cat('No P-value < 0.05 for joint test')}";
        end
        geneRphe = [""; "### Rare-variant analysis (RVA)"; ""; "```{r}"; vertcat(geneRphe{:}); "```"; ""];
        rep.genes{i} = [rep.genes{i}; geneRphe];
    end
    
    % add gene name heading
    rep.genes{i} = [""; "## " + opts.genes(i); ""; rep.genes{i}; ""; rep.sep];
end

% check Phenome-wide for single variants
if isfield(opts, 'snps')
    snps = opts.snps;
    xfiles = getfilenames(opts.dir, "xlsx").xlsx;
    xfiles(~startsWith(xfiles, "PheWide.")) = [];
    for j = 1:numel(snps.snp)
        idx = contains(xfiles, "." + snps.snp(j) + ".");
        file = xfiles(idx);
        if isempty(file), continue; end
        [~, name] = fileparts(file);
        name = regexprep(name, "^PheWide.", "");
        name = split(name, ".");
        geneName = name(1);
        name = join(name, " ");
        geneRphe = ({});
        geneRphe{j}(1,1) = "tab" + j + ' = readxl::read_xlsx("' + ...
            replace(file, filesep, [filesep, filesep]) + '")';
        geneRphe{j}(2,1) = "tab" + j + "=tab" + j + "[tab" + j + "$`P-adj` <= 0.05, ]";
        geneRphe{j}(3,1) =  "if (dim(tab" + j + ")[1] > 0) {";
        geneRphe{j}(4,1) = "DT::datatable(tab" + j + ", extensions" + ...
            "= c('FixedColumns','FixedHeader'),options = list(scrollX" + ...
            "= T,paging = T,fixedHeader=T, fixedColumns = " + ...
            "list(leftColumns = 2), pageLength = 5, " + ...
            "autoWidth = TRUE), rownames = T, caption = " + ...
            "htmltools::tags$caption(style = 'caption-side: " + ...
            "bottom; text-align: center;', htmltools::em('Phenome-wide summary statistics for " + ...
            name + "')))";
        geneRphe{j}(5,1) = "} else {cat('No P-value < 0.05')}";
        geneRphe = [""; "### Phenome-Wide"; ""; "```{r}"; vertcat(geneRphe{:}); "```"; ""];
        rep.genes{i+j} = [""; "## " + geneName; geneRphe; rep.sep];
    end
end

% compile final report
report = [rep.header; ""; rep.intro; ""; rep.sep; rep.res; vertcat(rep.genes{:})];
writematrix(report, fullfile(opts.dir, opts.report_name + ".Rmd"), "FileType", "text", "QuoteStrings", "none")

% workaround the R issue with render
rcode = "rmarkdown::render('" +  replace(fullfile(opts.dir, ...
    opts.report_name + ".Rmd"), filesep, [filesep,filesep]) + "')";
MATLAB2Rconnector(opts.report_name + ".r", "code", rcode);

tmpFigs = vertcat(tmpFigs{:});
try arrayfun(@delete, tmpFigs(tmpFigs ~= "")); catch, end
end % END

%% ------------------------------------------------------------------------
function c = rVec(c)
c = "c('" + join(c, "','") + "')";
end

%% ------------------------------------------------------------------------
function mergeDifferentGSSFolders(fos, genes, opts)

arguments
    fos {mustBeFolder}
    genes {mustBeVector, mustBeText}
    opts.outdir {mustBeTextScalar} % output folder
end

outdir = opts.outdir;
if outdir == "" || ismissing(outdir), outdir = pwd; end

for k1 = 1:numel(genes)
    for k2 = 1:numel(fos)
        xls.("g"+k1){k2,1} = getfilenames(fullfile(fos(k2), genes(k1)), "xlsx", fullpath=true).xlsx;
        regFolder = fullfile(fos(k2), genes(k1), "regional");
        if isfolder(regFolder)
            reg.("g"+k1){k2,1} = getfilenames(regFolder, fullpath=true).x_;
        else
            reg.("g"+k1){k2,1} = missing;
        end
    end
end

for k = 1:numel(genes)
    xls.("g"+k) = vertcat(xls.("g"+k){:});
    files = xls.("g"+k);
    mkdir(fullfile(outdir, genes(k)));
    if ~all(cellfun(@(x)any(ismissing(x)), reg.("g"+k)))
        mkdir(fullfile(outdir, genes(k), "regional"))
    end
    rva = files(files.endsWith(".RVA.xlsx"));
    sva = files(files.endsWith(".SVA.xlsx"));
    figs = vertcat(reg.("g"+k){:});
    figs(ismissing(figs)) = [];
    if ~isempty(figs)
        [~, names, ext] = fileparts(figs);
        names = names + ext;
        [~, ~, ~, idx] = duplicates(names);
        figs = figs(idx);
        names = names(idx);
        arrayfun(@(x,y)copyfile(x,y), figs, fullfile(outdir, genes(k), "regional", names))
    end

    % phewide files
    phefiles = files(files.contains(filesep + "PheWide."));
    if ~isempty(phefiles)
        [~, names, ext] = fileparts(phefiles);
        names = names + ext;
        arrayfun(@(x,y)copyfile(x,y), phefiles, fullfile(outdir, genes(k), names))
    end

    % merge sva
    if ~isempty(sva)
        svatab = cell(numel(sva), 1);
        for j = 1:numel(sva)
            svatab{j, 1} = readtable(sva(j), VariableNamingRule="preserve");
            if width(svatab{j}) < 15
               svatab{j}{:, ["AF.case", "AF.control", "N.case"]} = nan(height(svatab{j}), 3);
            end
            idx = ismember(colnames(svatab{j}), ["β|OR", "β", "OR"]);
            svatab{j}.Properties.VariableNames{idx} = 'β|OR';
        end
    
        writetable(vertcat(svatab{:}), ...
            fullfile(outdir, genes(k), genes(k) + ".SVA.xlsx"))
    end

    % merge rva
    if isempty(rva), continue; end
    data_key = cell(numel(rva), 1);
    for j = 1:numel(rva)
        data_key{j, 1} = array2table([sheetnames(rva(j)), ...
            matlab.lang.makeValidName(sheetnames(rva(j)))], ...
            VariableNames=["s", "k"]);
    end
    data_key = unique(vertcat(data_key{:}), 'rows', 'stable');
    
    str = struct;
    for i = 1:numel(rva)
        snames = sheetnames(rva(i));
        [~, idx] = ismember(snames, data_key.s); idx(idx < 1) = [];
        key = data_key(idx, :);

        legend_sheet_idx = key.s.lower.startsWith("legend");
        
        if any(legend_sheet_idx)
            str.(key.k(legend_sheet_idx)){i, 1} = readmatrix(rva(i), ...
                Sheet=key.s(legend_sheet_idx), ...
                OutputType="string");
        end
        key(legend_sheet_idx, :) =  [];
        for j = 1:height(key)
            str.(key.k(j)){i, 1} = readtable(rva(i), VariableNamingRule="preserve", Sheet=key.s(j));
        end
    end

    fi = string(fieldnames(str));
    for j = 1:numel(fi)
        if istable(str.(fi(j))), continue; end
        try
            str.(fi(j)) = vertcat(str.(fi(j)){:});
        catch % tables of different size
            tmp = str.(fi(j));
            tmpjoin = outerjoin(tmp{1}, tmp{2}, MergeKeys=true);
            for k2 = 3:numel(tmp)
                tmpjoin = outerjoin(tmpjoin, tmp{k2}, MergeKeys=true);
            end
            str.(fi(j)) = tmpjoin;
        end
    end

    overall_sheet_idx = data_key.s.lower.startsWith("overall.info");
    legend_sheet_idx = data_key.s.lower.startsWith("legend");
    data_key = [data_key(overall_sheet_idx, :);
        data_key(~(overall_sheet_idx|legend_sheet_idx),:);
        data_key(legend_sheet_idx, :)];
    legend_sheet_idx = data_key.s.lower.startsWith("legend");
    lsheet = data_key.k(end);
    str.(lsheet) = unique(str.(lsheet), "rows", "stable"); % Legend table
    str.(lsheet) = array2table(str.(lsheet));
    for i = 1:height(data_key)
        if legend_sheet_idx(i), writeFlag = false; else, writeFlag = true; end
        
        writetable(str.(data_key.k(i)), ...
            fullfile(outdir, genes(k), genes(k) + ".RVA.xlsx"), Sheet=data_key.s(i), ...
            WriteVariableNames=writeFlag, AutoFitWidth=true);
    end

end

end
