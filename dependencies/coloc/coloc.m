function coloc(gwas, index, opts)
% a wrapper for coloc R package using either eQTL Catalogue database or
% another summary stat file

% Notes:
% coloc.abf() assumes single causal variant
%   H0: neither trait has a genetic association in the region
%   H1: only trait 1 has a genetic association in the region (eQTL/2nd GWAS)
%   H2: only trait 2 has a genetic association in the region (GWAS)
%   H3: both traits are associated, but with different causal variants
%   H4: both traits are associated and share a single causal variant
% 
% Default priors: 
% p1 (prior probability that a SNP is associated with the
% eQTL) = 1e−4; p2 (prior probability that a SNP is associated with the
% GWAS trait) = 1e−4; and p12 (prior probability that a SNP is
% associated with both GWAS trait and eQTL) = 1e−5
% 
% fine-mapped eQTL: http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/credible_sets/
% 
% The error "The estimated prior variance is unreasonably large" in susie,
% means mismatch between LD and summary stat (in this case eQTLs): https://github.com/stephenslab/susieR/issues/151
% 
% @19/08/2022: some bugs were fixed.

% TODO: GTEX API for LD matrix:
% g = webread("https://gtexportal.org/rest/v1/reference/gene?geneId=ENSG00000035681&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38&pageSize=250", options);
% x = webread("https://gtexportal.org/rest/v1/dataset/ld?format=json&gencodeId=ENSG00000035681.7&datasetId=gtex_v8", options);
% x = string(horzcat(x.ld{:}))';


arguments
    gwas {mustBeFile} % GWAS summary stat file
    index % index variants (clumped, COJO or any other filters): can be table, string or a file

    % table for the other GWAS summary stats. If provided, ignores eQTL
    % arguments. Note that at the moment it's assumed that both GWAS are
    % from the same study/cohort (e.g. UK Biobank), therefore doesn't
    % check for ambiguous/palindromic snps (and therefore, no strand
    % flipping is done). For cross-platforms, follow this:
    % https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html
    opts.gss {mustBeA(opts.gss, "table")}
    opts.workers (1,1) double = 50 % for parpool
    opts.output {mustBeTextScalar, mustBeNonzeroLengthText}

    % for eQTL datasets
    opts.tissue {mustBeText} = "" %["adipose_subcutaneous", "liver"] % tissues for colocalization (case insensitive); if left empty all tissues will be tested (not recommended)
    opts.geneid {mustBeText, mustBeVector} = "" % consider only these genes
    opts.win (1,1) double = 500 % in kb. Window around the index variants and genes
    opts.liftover (1,1) logical = true % liftovers GWAS summary stats from GRCh37 -> GRCh38
    opts.quant {mustBeMember(opts.quant, ["ge", "exon", "microarray", "tx", "txrevise"])} = "ge" % quantification method. default: RNA-Seq (ge). https://www.ebi.ac.uk/eqtl/api-docs/
    
    % study to consider. 'tissue' overrides this option
    opts.study {mustBeMember(opts.study, ["all", "Alasoo_2018", "BLUEPRINT", ...
    	"BrainSeq", "Braineac2", "CAP", "CEDAR", "CommonMind", "FUSION", ...
    	"Fairfax_2012", "Fairfax_2014", "GENCORD", "GEUVADIS", "GTEx", ...
    	"HipSci", "Kasela_2017", "Lepik_2017", "Naranbhai_2015", ...
        "Nedelec_2016", "Peng_2018", "PhLiPS", "Quach_2016", "ROSMAP", ...
    	"Schmiedel_2018", "Schwartzentruber_2018", "Steinberg_2020", ...
    	"TwinsUK", "Young_2019", "iPSCORE", "van_de_Bunt_2015", "GTEx_V8"])} = "all"

    opts.prior (1, 3) double = [1e-4 1e-4 1e-5] % priors for coloc.abf: in order: P1 P2 P12 
    opts.verbose (1, 1) logical = true
    opts.sdY (1,1) double = 1 % SD of dependent variable. If NaN, will be estimated from MAF and sample size. If IRNT transformation was done, set to 1.
    opts.type (1,1) {mustBeMember(opts.type, ["quant", "cc"])} = "quant" % quant: continuous, cc: binary response.
    opts.method {mustBeMember(opts.method, ["abf", "susie"])} = "abf" % colocalization method, if susie fails, will use abf
    opts.fdrcutoff (1,1) double = 0.1 % FDR cutoff of eQTL dataset
    opts.egenedir {mustBeFolder} % for GTEx to choose list of eGenes around a GWAS index SNP based on the 'fdrcutoff'.

    % for susie
    opts.lddir {mustBeFolder} % needed for susie (from PolyFun paper) only if ldmethod is "external"
    opts.ldmethod {mustBeMember(opts.ldmethod, ["internal", "external"])} = "internal" % "internal": calculate LD matrix, "external": use precalculated LD matrix 
    opts.ldsample {mustBeNumeric, mustBeVector} = getQCEID(3, false) % samples for LD matrix (only if ldmethod is internal)
    opts.bgenhome {mustBeTextScalar} % bgen files for LD matrix (only if ldmethod is internal)
    opts.coverage (1,1) double {mustBeInRange(opts.coverage, 1e-6, 0.9999)} = 0.9 % default is 0.95 -> 95% credible set
    opts.cachedir {mustBeFolder} % tabix cache local files. If not provided, tabix tries to fetch data from server (not recommended in using in parallel).
end

if isempty(gcp('nocreate'))
    parpool("Processes", opts.workers);
else
    opts.workers = min(opts.workers, gcp('nocreate').NumWorkers);
end

if isfield(opts, 'gss'); opts.eqtl = false; else; opts.eqtl = true; end
if ~isfield(opts, 'output')
    [~, opts.output] = fileparts(gwas);
    opts.output = string(opts.output) + ".coloc";
end

% read only summary stat of variants around the index variants
if opts.verbose, fprintf('reading gwas summary stat...\n'); end
[gwastab, ~, index] = readGWASfile(gwas, 'index', index, 'parallel', true, ...
    'light', false, 'full', true, 'n', true, 'win', opts.win, "maxIndexSize", 40);

% 'chr' between index and gwastabs may be have different data type
if isnumeric(gwastab.chr) && ~isnumeric(index.chr)
    index.chr = double(index.chr);
elseif ~isnumeric(gwastab.chr) && isnumeric(index.chr)
    index.chr = string(index.chr);
end

if ~any(ismember(index.Properties.VariableNames, {'pos', 'chr'}))
    [~, idx] = ismember(index.snp, gwastab.snp);
    index.pos = gwastab.pos(idx);
    index.chr = gwastab.chr(idx);
end

% keep only bi-allelic variants
if opts.verbose, fprintf('prunning gwas summary stat...\n'); end
snpcol = ismember(gwastab.Properties.VariableNames, {'annotate'});
gwastab.Properties.VariableNames(snpcol) = {'snp'};
gwastab.id = gwastab.chr + ":" + gwastab.pos;
dups = duplicates(gwastab.id);
gwastab(ismember(gwastab.id, dups), :) = [];

% do liftOver using rtracklayer R pacakge
if opts.liftover
    if opts.verbose, fprintf('performing liftOver...\n'); end
    gwastab = liftover(gwastab, 'useEnsemble', true);
    % remove those with failed mapping
    if any(isnan(gwastab.pos))
        if opts.verbose
            fprintf('\tremoved %d snp due to failed lift over.\n', sum(isnan(gwastab.pos)))
        end
        gwastab(isnan(gwastab.pos), :) = [];
    end

    % after lift over, some variants may have the same position (can be due
    % to non-mapped e.g. rs7247683 in GRCh38). These problematic variants
    % should be removed also 
    gwastab.id = gwastab.chr + ":" + gwastab.pos;
    dups = duplicates(gwastab.id);
    gwastab(ismember(gwastab.id, dups), :) = [];
end

% if opts.method ~= "susie"
%     gwastab = modifymaf(gwastab, 'afreq', opts.verbose);
% end

% loop over each index variant and extract a subset of GWAS summary stat
if opts.liftover
    poscol = "pos_old";
else
    poscol = "pos";
end

if opts.eqtl
    % get tissues and studies from eQTL catalogue
    tryMe = true;
    cnt = 1;
    while tryMe
        try
            tabix_paths = readtable("https://raw.githubusercontent.com/eQTL" + ...
                "-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths" + ...
                ".tsv", 'FileType', 'text', 'TextType', 'string', ...
                'VariableNamingRule', 'preserve', 'NumHeaderLines', 0);
            
            % @11APR2023: some changes have been made to eQTL Catalogue git
            % repo
            if any(colnames(tabix_paths) == "study_label")
                tabix_paths.study = tabix_paths.study_label;
                tabix_paths.qtl_group = tabix_paths.sample_group;

                % ftp paths have been changed. To avoid inconsistencies with
                % already downloaded GTEx files in the cache_dir, I modify ftp
                % paths, but this needs to be properly taken care of for future
                % use with other stduies
                if opts.verbose
                    warning("FTP paths have been changed manually")
                end
    
                [ftp_pth, ftp_name, ftp_ext1] = fileparts(tabix_paths.ftp_path);
                ftp_ext2 = extract(ftp_name, "." + alphanumericsPattern + ".tsv");
                ftp_name = tabix_paths.study + "_" + tabix_paths.quant_method + "_" + tabix_paths.qtl_group;
                tabix_paths.ftp_path_raw = tabix_paths.ftp_path;
                tabix_paths.ftp_path = ftp_pth + "/" + ftp_name + ftp_ext2 + ftp_ext1;
            end

            tryMe = false;
        catch exception
            pause(1) % wait for 1 sec!
            if (cnt > 5)
                throw(exception)
            end
            cnt = cnt + 1;
        end
    end
    
    % in eQTL catalogue release 4: GTEx shows GTEX V.8 and GTEx_V8 is original
    % data
%     gtex8 = readtable("https://raw.githubusercontent.com/eQTL" + ...
%         "-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported" + ...
%         ".tsv", 'FileType', 'text', 'TextType', 'string', ...
%         'VariableNamingRule', 'preserve', 'NumHeaderLines', 0);
%     gtex8.sample_size = nan(height(gtex8), 1);
%     tabix_paths = [tabix_paths; gtex8];
    
    opts.tissue = lower(opts.tissue); 
    
    tabix_paths(tabix_paths.quant_method ~= opts.quant, :) = []; % only keep this quantification method
    if opts.tissue ~= "" 
        idx = contains(lower(tabix_paths.tissue_label), opts.tissue) |...
            contains(lower(tabix_paths.qtl_group), opts.tissue);
        if any(colnames(tabix_paths) == "tissue_ontology_term")
            idx = idx | contains(lower(tabix_paths.tissue_ontology_term), opts.tissue);
        end
        tabix_paths = tabix_paths(idx, :);
    end
    if opts.study ~= "all"
       tabix_paths(~ismember(tabix_paths.study, opts.study), :) = []; 
    end

    % check eGenes
    tabix_paths = fetchEgeneFiles(tabix_paths, opts);
end

res = cell(numel(index.snp), 1);

% for temp files
if fileparts(opts.output) == ""
    opts.dir = pwd;
else
    opts.dir = fileparts(opts.output);
end

if opts.method == "abf" && ~isempty(gcp('nocreate')) && ~opts.eqtl % use parfor for GWAS files (non-eQTL)
    parfor i = 1:numel(index.snp)
        optsIn = opts;
        indexIn = index;
        gwastabIn = gwastab;

        if optsIn.verbose
            fprintf('processing region %d of %d\n', i, numel(indexIn.snp))
        end
        start = indexIn.pos(i) - optsIn.win*1e3;
        stop = indexIn.pos(i) + optsIn.win*1e3;
        if start <= 0; start = 1; end
        idx = gwastabIn.(poscol) >= start & gwastabIn.(poscol) <= stop & gwastabIn.chr == indexIn.chr(i);
        gwas = gwastabIn(idx, :);
        if isempty(gwas)
            fprintf('\tskipped (no summary stat was found within the region)\n')
            fprintf('-------------------------------------------------\n\n')
            continue
        end
        region = gwas.chr(1) + ":" + min(gwas.pos) + "-" + max(gwas.pos);
        
        gssidx = optsIn.gss.pos >= start & optsIn.gss.pos <= stop & optsIn.gss.chr == indexIn.chr(i);
        gss = optsIn.gss(gssidx, :);
        vrbs = optsIn.verbose; optsIn.verbose = false;
        res{i} = doColoc(gwas, gss, optsIn);
        optsIn.verbose = vrbs;
        if ~isempty(res{i})
            res{i}.region = repmat(region, height(res{i}), 1);
            res{i}.chr = repmat(indexIn.chr(i), height(res{i}), 1);
            res{i}.pos = repmat(indexIn.pos(i), height(res{i}), 1);
            if any(strcmp("snp", indexIn.Properties.VariableNames))
                res{i}.index = repmat(indexIn.snp(i), height(res{i}), 1);
            else
                try
                    [~, minpIdx] = min(gss.p);
                    res{i}.index = repmat(gss.snp(minpIdx), height(res{i}), 1);
                    res{i}.chr = repmat(gss.chr(minpIdx), height(res{i}), 1);
                    res{i}.pos = repmat(gss.pos(minpIdx), height(res{i}), 1);
                catch % snp/p missing?
                end
            end
        end
        
    end 

else
    for i = 1:numel(index.snp)
        if opts.verbose
            fprintf('processing region %d of %d\n', i, numel(index.snp))
        end
        start = index.pos(i) - opts.win*1e3;
        stop = index.pos(i) + opts.win*1e3;
        if start <= 0; start = 1; end
        idx = gwastab.(poscol) >= start & gwastab.(poscol) <= stop & gwastab.chr == index.chr(i);
        
        % when "egenedir" is present, multiple regions are needed depending on the eGenes around each GWAS index
        if opts.eqtl 

            % is there are some datasets with missing eGenes, use gnomAD API to fetch genes around this GWAS index
            if any(~cellfun(@istable, tabix_paths.egene)) 
                regiongenes = gnomad('start', min(gwastab.pos(idx)), ...
                    'stop', max(gwastab.pos(idx)), 'chrom', ...
                    string(index.chr(i)), 'compact', true, ...
                    'dataset', 'gnomad_r3', 'reference_genome', 'GRCh38');
                regiongenes = struct2table(regiongenes.data.region.genes);
                regiongenes = renamevars(regiongenes, ["gene_id", "symbol", "stop"], ["id", "name", "end"]);
                regiongenes = convertvars(regiongenes, ["name", "id"], @string);
            end
            
            [rstart, rstop] = deal(nan(height(tabix_paths), 1));
            tabix_paths.regiongenes = deal(cell(height(tabix_paths), 1));
            tabix_paths.region = strings(height(tabix_paths), 1); % regions around eQTL genes (eGenes)
            for j = 1:height(tabix_paths)
                if istable(tabix_paths.egene{j})
                    regiongenesTmp = tabix_paths.egene{j}(...
                        tabix_paths.egene{j}.chr == string(index.chr(i)) ...
                        & max(tabix_paths.egene{j}.start, min(gwastab.pos(idx))) <= min(tabix_paths.egene{j}.end, max(gwastab.pos(idx))), :); 
                else % single pvalue cut-off: use gnomAD API genes
                    regiongenesTmp = regiongenes;        
                end
                
                % check custom genes to keep
                if all(opts.geneid ~= "") 
                    rmidx = ~ismember(regiongenesTmp.id, opts.geneid);
                    if ~any(rmidx)
                        rmidx = ~ismember(regiongenesTmp.name, opts.geneid);
                    end
                    regiongenesTmp(rmidx, :) = [];
                end

                if isempty(regiongenesTmp) % no gene is left for this dataset due to either 'geneid' option or no FDR corrected eGenes within the region
                    continue
                end

                rstart(j) = min(regiongenesTmp.start) - opts.win*1e3;
                rstop(j) = max(regiongenesTmp.end) + opts.win*1e3;
                if rstart(j) <= 0; rstart(j) = 1; end
                tabix_paths.region(j) = string(index.chr(i)) + ":" + rstart(j) + "-" + rstop(j);
                tabix_paths.regiongenes{j} = regiongenesTmp(:, ["id", "name"]);
            end

            % remove datasets with empty 'regiongenes' (see above
            % condition)
            rmidx = cellfun(@isempty, tabix_paths.regiongenes);
            rstart(rmidx) = []; rstop(rmidx) = [];
            tabix_paths(rmidx, :) = [];
            
            if isempty(rstart) % no gene has remained for all datasets within this GWAS index region!
                idx = [];
            else
                % now expand the region around GWAS index to cover genes within
                % 'region' variable
                start = min(rstart);
                stop = max(rstop);
                if start <= 0; start = 1; end
                idx = gwastab.pos >= start & gwastab.pos <= stop & gwastab.chr == index.chr(i);
            end

        end
        
        gwas = gwastab(idx, :);
        if isempty(gwas)
            fprintf('\tskipped (no summary stat was found within the region)\n')
            fprintf('-------------------------------------------------\n\n')
            continue
        end
        region = gwas.chr(1) + ":" + min(gwas.pos) + "-" + max(gwas.pos);

        % check if backup files exist (in case of crashes)
        if isfile(fullfile(opts.dir, "chr" + replace(region, [":", "-"], "_") + ".mat")) && opts.eqtl
            res{i} = load(fullfile(opts.dir, "chr" + replace(region, [":", "-"], "_") + ".mat")).eqtl;
            continue
        end
        
        if opts.method == "susie"
            gwas = getLDmat(gwas, opts);
        end
        
        if ~opts.eqtl % 2nd GWAS summary stat
            gssidx = opts.gss.pos >= start & opts.gss.pos <= stop & opts.gss.chr == index.chr(i);
            gss = opts.gss(gssidx, :);
            res{i} = doColoc(gwas, gss, opts);
            if ~isempty(res{i})
                res{i}.region = repmat(region, height(res{i}), 1);
                res{i}.chr = repmat(index.chr(i), height(res{i}), 1);
                res{i}.pos = repmat(index.pos(i), height(res{i}), 1);
                if any(strcmp("snp", index.Properties.VariableNames))
                    res{i}.index = repmat(index.snp(i), height(res{i}), 1);
                else
                    try
                        [~, minpIdx] = min(gss.p);
                        res{i}.index = repmat(gss.snp(minpIdx), height(res{i}), 1);
                        res{i}.chr = repmat(gss.chr(minpIdx), height(res{i}), 1);
                        res{i}.pos = repmat(gss.pos(minpIdx), height(res{i}), 1);
                    catch % snp/p missing?
                    end
                end
            end
        else
            % fetch eQTL data for this subset
            eqtl = tabix_paths;
            [eqtl.coloc, tabix_paths] = fetcheqtlCatalogue(tabix_paths, opts);
            eqtl.regiongenes = tabix_paths.regiongenes;
            eqtl.egene = [];
        
            % remove those eQTL datasets with no significant SNPs within this region
            emptyEqtlDatasets = cellfun(@isempty, eqtl.coloc); 
            eqtl(emptyEqtlDatasets, :) = []; 
            
            eqtl.coloc = doColoc(gwas, eqtl.coloc, opts);
            eqtl.region = repmat(region, height(eqtl), 1);
            if ~isempty(eqtl.coloc)
                eqtl.region = repmat(region, height(eqtl), 1);
                eqtl.chr = repmat(index.chr(i), height(eqtl), 1);
                eqtl.pos = repmat(index.pos(i), height(eqtl), 1);
                if any(strcmp("snp", index.Properties.VariableNames))
                    eqtl.index = repmat(index.snp(i), height(eqtl), 1);
                else
                    eqtl.index = repmat("", height(eqtl), 1);
                end
            end
            res{i} = eqtl;
            save(fullfile(opts.dir, "chr" + replace(region, [":", "-"], "_") + ".mat"), 'eqtl')
            
        end
    
        if opts.verbose
            fprintf('-------------------------------------------------\n\n')
        end
        
    end 
end

res(cellfun(@isempty, res)) = [];

if opts.eqtl
    resp = prunecolocResult(res);
else
    res = vertcat(res{:});
    resp = res(res.("PP.H4.abf") >= 0.8, :);
end

save(opts.output + ".mat", 'res', 'resp')

end % END

%% subfunctions ===========================================================
function [coloc, pths] = fetcheqtlCatalogue(pths, opts, verbose)
if nargin < 4
    verbose = true;
end

% prepare some inputs
coloc = cell(height(pths), 1);
ftp_path = pths.ftp_path;
region = pths.region;
study = pths.study;

% check if eGenes are missing for some datasets, use single p-value cutoff to define eGenes
singleCutoffIdx = ~cellfun(@istable, pths.egene);
if any(singleCutoffIdx) 
    singleCutoff = pths.egene{find(singleCutoffIdx, 1)};
else
    singleCutoff = 1e-6; % to avoid parfor error 
end

% genes within the GWAS index region (either eGenes or all genes in the
% region; singleCutoff is used for the latter)
regiongenes = pths.regiongenes;

if ~isfield(opts, 'cachedir')
    cachedir = pwd;
else
    cachedir = opts.cachedir;
end
cacheExists = isfolder(fullfile(cachedir, study));
if sum(~cacheExists) > 10
    N = max(round(opts.workers/4), 1);
else
    N = opts.workers;
end

% step 1: find eGenes: cis-eGenes with at least one FDR-corrected
% significant variant

% From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5142122/: eGenes are
% genes that have at least one significant variant (p value ≤ 10−5 when
% corrected for multiple hypothesis) associated with the expression of that
% gene. In another word, they are genes with at least one significant eQTL.
% In GTEx project, eGenes are based on FDR corrected p-values per each
% tissue (i.e. using all associations with each gene in each tissue file,
% which is absent from eQTL-catalogue).
% 
% As also mentioned here
% (https://www.nature.com/articles/s41467-021-24051-6#Sec7), a single
% p-value cut-off within the region around the GWAS locus can be used to
% find eGenes (molecular traits). So, this seems reasonable to employ FDR
% corrected p-value (q-value) "around the GWAS locus region", instead of
% taking all SNPs into account. The caveat is that this means the region
% (+/- bp around GWAS locus) affects the q-value, and hence the number of
% genes to be tested. The more accurate approach would be to first define
% all eGenes (https://www.nature.com/articles/nature24277 and
% https://www.medrxiv.org/content/10.1101/2020.12.26.20248491v3.full).
% However since it's computationally intensive this is cutrrently what we
% use:
%   1-check if 'egenedir' exists, and use 'fdrcutoff' to define list of
%   eGenes for genes around the GWAS index variant
% 
%   2- if 'egenedir' not present a single cut-off of <= 1E-6 within the
%   region will be used (for non-GTEx databases mainly), as used here:
%   https://www.nature.com/articles/s41467-021-24051-6

parfor (k = 1:numel(ftp_path), N)
% for k = 1:numel(ftp_path)
    if verbose
        fprintf('\tfetching data %d of %d\n', k, numel(ftp_path))
    end
    
    if cacheExists(k)
        try
            tmp = tabix(ftp_path(k), region(k), 'cachedir', fullfile(cachedir, study(k)), 'offline', 'lazy');
        catch
            fprintf('\t%d-problematic region: %s\ncache file:%s\n', k, string(region(k)), fullfile(cachedir, study(k)))
            tmp = tabix(ftp_path(k), region(k), 'cachedir', fullfile(cachedir, study(k)), 'offline', 'lazy');
        end
    else
        tmp = tabix(ftp_path(k), region(k), 'offline', 'lazy');
    end

    if isempty(tmp)
        continue
    end

    % remove rsid duplicates
    rsid = tmp.rsid; tmp.rsid = [];
    [tmp, idx] = unique(tmp, 'rows');
    tmp.rsid = rsid(idx);

    if singleCutoffIdx(k)
        % use single p-value cutoff (default = 1e-6) to find eGenes
        tmp = groupfilter(tmp, 'gene_id', @(x) any(x <= singleCutoff), 'pvalue');
        
        % modify 'regiongenes', and keep only single p-value based eGenes
        rmidx = ~ismember(regiongenes{k}.id, tmp.gene_id);
        regiongenes{k}(rmidx, :) = [];

        % see if there are genes in 'tmp' table not present in
        % 'regiongenes', this may happen due to different GENCODE versions
        % between gnomAD (regiongenes come from gnomAD API) and eQTL
        % datasets
        missingGenes = unique(setdiff(tmp.gene_id, regiongenes{k}.id));
        if ~isempty(missingGenes)
            missingGenes = table(missingGenes, ...
                strings(numel(missingGenes)), 'VariableNames', {'id', 'name'});
            regiongenes{k} = [regiongenes{k}; missingGenes];
        end

    else % dataset has already the list of eGenes
        id = regiongenes{k}.id;
        tmp(~ismember(tmp.gene_id, id), :) = [];
    end
    coloc{k} = tmp;
end

pths.regiongenes = regiongenes;

end

%% ------------------------------------------------------------------------
function res = doColoc(gwasall, alleqtl, opts)
res = [];
gwasall.id = gwasall.chr + ":" + gwasall.pos;

if ~opts.eqtl
    if opts.verbose, fprintf('\trunning coloc for both traits\n'); end
    alleqtl.id = alleqtl.chr + ":" + alleqtl.pos;
    alleqtl = unique(alleqtl, 'rows');
    dups = duplicates(alleqtl.id);
    alleqtl(ismember(alleqtl.id, dups), :) = [];
    if isempty(alleqtl); return; end

    % joint summary stats between eQTL and GWAS
    [f1,f2] = ismember(gwasall.id, alleqtl.id); f2(f2<1) = [];
    gwasall = gwasall(f1, :); alleqtl = alleqtl(f2, :); 
    
    % find flipped and different variants between eQTL and GWAS
    flippedIdx = find(alleqtl.allele1 ~= gwasall.allele1);
    
    % check if their other allele is the same, then flip the beta sign
    % in GWAS
    toFlipIdx = flippedIdx(alleqtl.allele0(flippedIdx) == gwasall.allele1(flippedIdx));
    newA1 = gwasall.allele0(toFlipIdx); 
    gwasall.allele0(toFlipIdx) = gwasall.allele1(toFlipIdx);
    gwasall.allele1(toFlipIdx) = newA1;
    gwasall.beta(toFlipIdx) = -1.*gwasall.beta(toFlipIdx);
    gwasall.afreq(toFlipIdx) = 1 - gwasall.afreq(toFlipIdx);
    
    % remove those failed to compare
    rmIdx = (alleqtl.allele1 ~= gwasall.allele1) | (alleqtl.allele0 ~= gwasall.allele0);
    alleqtl(rmIdx, :) = []; gwasall(rmIdx, :) = [];
    if isempty(alleqtl), return; end

    if opts.method == "susie"
        gwasall.Properties.UserData = gwasall.Properties.UserData(f1, f1);
        % flip R sign
        gwasall.Properties.UserData(toFlipIdx, :) = ...
            -gwasall.Properties.UserData(toFlipIdx, :); 
        gwasall.Properties.UserData(:, toFlipIdx) = ...
            -gwasall.Properties.UserData(:, toFlipIdx); 
        gwasall.Properties.UserData(rmIdx, :) = [];
        gwasall.Properties.UserData(:, rmIdx) = [];
    end
    
    % prepare for coloc package
    try gwasN = gwasall.n(1); catch; gwasN = nan; end
    if isnan(opts.sdY)
        gwasall = gwasall(:, ["id", "afreq", "beta", "se"]); % keep only required columns
    else
        gwasall = gwasall(:, ["id", "beta", "se"]); % keep only required columns
    end
    gwasall.se = gwasall.se.^2;
    
    if ~isnumeric(alleqtl.afreq); alleqtl.afreq = double(string(alleqtl.afreq)); end
    if ~isnumeric(alleqtl.se); alleqtl.se = double(string(alleqtl.se)); end
    if ~isnumeric(alleqtl.beta); alleqtl.beta = double(string(alleqtl.beta)); end

    try eqtlN = alleqtl.n(1); catch; eqtlN = nan; end
    alleqtl = alleqtl(:, ["id", "afreq", "beta", "se"]);
    alleqtl.se = alleqtl.se.^2;
    nanidx = isnan(alleqtl.se);
    if any(nanidx)
        alleqtl(nanidx, :) = [];
        gwasall(nanidx, :) = [];
    end
    res = callColoc(gwasall, alleqtl, opts, gwasN, eqtlN, "NA");
    res = convertvars(res, [1, 4:8], "double");
    return
end

N = opts.workers;
if opts.method == "susie"
    N = round(N/2); % to avoid out of memory issues
end
res = cell(numel(alleqtl), 1);
for i = 1:numel(alleqtl) % loop over studies/tissue combinations
    if opts.verbose
        fprintf('\trunning coloc for eqtl dataset %d of %d (%d snps)\n', i, numel(alleqtl), height(alleqtl{i}))
    end

    alleqtl{i}.id = alleqtl{i}.chromosome + ":" + alleqtl{i}.position;
    genes = unique(alleqtl{i}.gene_id);
    out = cell(numel(genes), 1);
    [eqtlN, gwasN] = deal(zeros(numel(genes), 1));

    [eqtl, gwas] = deal(cell(numel(genes), 1));

    for j = 1:numel(genes) % loop over genes

        eqtl{j} = alleqtl{i}(alleqtl{i}.gene_id == genes(j), :); % for this gene

        % remove multi-allelic variants from eQTL data
        dups = duplicates(eqtl{j}.id);
        eqtl{j}(ismember(eqtl{j}.id, dups), :) = [];

        if isempty(eqtl{j}); continue; end

        % joint summary stats between eQTL and GWAS
        [f1,f2] = ismember(gwasall.id, eqtl{j}.id); f2(f2<1) = [];
        gwas{j} = gwasall(f1, :); eqtl{j} = eqtl{j}(f2, :); 
        
        % find flipped and different variants between eQTL and GWAS
        flippedIdx = find(eqtl{j}.alt ~= gwas{j}.allele1);
        
        % check if their other allele is the same, then flip the beta sign
        % in GWAS
        toFlipIdx = flippedIdx(eqtl{j}.ref(flippedIdx) == gwas{j}.allele1(flippedIdx));
        newA1 = gwas{j}.allele0(toFlipIdx); 
        gwas{j}.allele0(toFlipIdx) = gwas{j}.allele1(toFlipIdx);
        gwas{j}.allele1(toFlipIdx) = newA1;
        gwas{j}.beta(toFlipIdx) = -1.*gwas{j}.beta(toFlipIdx);
        gwas{j}.afreq(toFlipIdx) = 1 - gwas{j}.afreq(toFlipIdx);
        
        % remove those failed to compare
        rmIdx = (eqtl{j}.alt ~= gwas{j}.allele1) | (eqtl{j}.ref ~= gwas{j}.allele0);
        eqtl{j}(rmIdx, :) = []; gwas{j}(rmIdx, :) = [];

        if opts.method == "susie"
            gwas{j}.Properties.UserData = gwasall.Properties.UserData(f1, f1);
            % flip R sign
            gwas{j}.Properties.UserData(toFlipIdx, :) = ...
                -gwas{j}.Properties.UserData(toFlipIdx, :); 
            gwas{j}.Properties.UserData(:, toFlipIdx) = ...
                -gwas{j}.Properties.UserData(:, toFlipIdx); 
            gwas{j}.Properties.UserData(rmIdx, :) = [];
            gwas{j}.Properties.UserData(:, rmIdx) = [];
        end
        
        % prepare for coloc package
        try
            gwasN(j) = gwas{j}.n(1);
        catch
            gwasN(j) = nan;
        end
        if isnan(opts.sdY)
            gwas{j} = gwas{j}(:, ["id", "afreq", "beta", "se"]); % keep only required columns
        else
            gwas{j} = gwas{j}(:, ["id", "beta", "se"]); % keep only required columns
        end
        gwas{j}.se = gwas{j}.se.^2;
        
        if ~isnumeric(eqtl{j}.an); eqtl{j}.an = double(string(eqtl{j}.an)); end
        if ~isnumeric(eqtl{j}.maf); eqtl{j}.maf = double(string(eqtl{j}.maf)); end
        if ~isnumeric(eqtl{j}.se); eqtl{j}.se = double(string(eqtl{j}.se)); end
        if ~isnumeric(eqtl{j}.beta); eqtl{j}.beta = double(string(eqtl{j}.beta)); end

        eqtlN(j) = round(mean(eqtl{j}.an, 'omitnan'))/2; % allele num / 2 -> N
        eqtl{j} = eqtl{j}(:, ["id", "maf", "beta", "se"]);
        eqtl{j}.se = eqtl{j}.se.^2;
        nanidx = isnan(eqtl{j}.se);
        if any(nanidx)
            eqtl{j}(nanidx, :) = [];
            gwas{j}(nanidx, :) = [];
        end
    end

    emptyIdx = cellfun(@isempty, eqtl);
    eqtl(emptyIdx) = []; gwas(emptyIdx) = []; genes(emptyIdx) = [];
    
    verbose = opts.verbose;
    parfevalOnAll(@warning, 0, 'off', 'all');
    if verbose
        fprintf('\t\trunning for %d genes\n', numel(genes))
    end

    parfor (j = 1:numel(genes), N)
%     for j = 1:numel(genes)
        out{j} = callColoc(gwas{j}, eqtl{j}, opts, gwasN(j), eqtlN(j), genes(j));
    end
    parfevalOnAll(@warning, 0, 'on', 'all');
    
    res{i} = vertcat(out{:});
    res{i} = convertvars(res{i}, [1, 4:8], "double");
    res{i} = sortrows(res{i}, 'PP.H4.abf', 'descend');
end
end

%% ------------------------------------------------------------------------
function out = callColoc(gwas, eqtl, opts, gwasN, eqtlN, gene)

name = getRandomName("colocR");
if isnan(opts.sdY); gmaf = gwas.afreq; end
gbeta = gwas.beta;
gse = gwas.se; 

if opts.method == "susie"
    LD = gwas.Properties.UserData;
end
clear gwas

id = cellstr(eqtl.id);
if opts.eqtl
    emaf = eqtl.maf;
else
    emaf = eqtl.afreq;
end
ebeta = eqtl.beta;
ese = eqtl.se; 
clear eqtl

if isnan(opts.sdY)
    save(name + ".mat", 'gmaf', 'gbeta', 'gse', 'id', 'emaf', 'ebeta', 'ese');
else
    save(name + ".mat", 'gbeta', 'gse', 'id', 'emaf', 'ebeta', 'ese');
end

if opts.method == "susie"
    save(name + ".mat", 'LD', '-append') 
end

rcode(1) = "library(coloc)";
rcode(2) = "data <- R.matlab::readMat('" + name + ".mat" + "')";

if ~isnan(opts.sdY)
    rcode(3) = "gwas <-  list(beta = as.vector(data$gbeta), " + ...
        "varbeta = as.vector(data$gse), sdY = " + opts.sdY + ", " + ...
        "type = '" + opts.type +"', snp = unlist(data$id, use.names = F))";
else
    rcode(3) = "gwas <-  list(beta = as.vector(data$gbeta), " + ...
        "varbeta = as.vector(data$gse), N = " + gwasN + ", MAF = " + ...
        "as.vector(data$gmaf), type = '" + opts.type + "', snp = unlist(data$id, use.names = F))";
end

if isnan(eqtlN) % assume SD == 1
    rcode(4) = "eqtl <-  list(beta = as.vector(data$ebeta), varbeta = " + ...
        "as.vector(data$ese), sdY = 1, MAF = as.vector(data$emaf)," + ...
        " type = 'quant', snp = unlist(data$id, use.names = F))";

else
    rcode(4) = "eqtl <-  list(beta = as.vector(data$ebeta), varbeta = " + ...
        "as.vector(data$ese), N = " + eqtlN + ", MAF = as.vector(data$emaf)," + ...
        " type = 'quant', snp = unlist(data$id, use.names = F))";
end

if opts.method == "susie"
    rcode(5) = "data$LD = data$LD + t(data$LD)";
    rcode(6) = "diag(data$LD) <- 1";
    rcode(7) = "colnames(data$LD) <- eqtl$snp";
    rcode(8) = "rownames(data$LD) <- eqtl$snp";
    rcode(9) = "gwas$LD = data$LD";
    rcode(10) = "eqtl$LD = data$LD";
end

rcode(11) = "rm(data)";

if opts.method == "susie"
    % maybe coverage 0.9 is more useful for coloc, also r_tol for semi-definitive check in LD matrix
    % Note: instead of maxit (default is 100), one could set
    % repeat_until_convergence to true, however, not sure how intensive the
    % problem can get. 
    rcode(12) = "res<- tryCatch({";
    rcode(13) = "gwas1 = coloc::runsusie(gwas, r_tol = 1e-4, coverage = " + opts.coverage + ", maxit = 1000)"; 
    rcode(14) = "eqtl1 = coloc::runsusie(eqtl, r_tol = 1e-4, coverage = " + opts.coverage + ", maxit = 1000)";
    rcode(15) = "coloc_res = coloc:::coloc.susie(eqtl1, gwas1" + ...
        ",p1 = " + opts.prior(1) + ", p2 = " + ...
        opts.prior(2) + ", p12 = " + opts.prior(3) + ")";
    rcode(16) = "coloc_res$summary}, error = function(e) {return(NULL)})";

    rcode(17) = "if (is.null(res)){";
    rcode(18) = "coloc_res = coloc::coloc.abf(dataset1 = eqtl, " + ...
        "dataset2 = gwas,p1 = " + opts.prior(1) + ", p2 = " + ...
        opts.prior(2) + ", p12 = " + opts.prior(3) + ")";
    rcode(19) = "res <- coloc_res$summary";
    rcode(20) = "res <- c(res, c(method = 'abf'))";
    rcode(21) = "} else { res$method = 'susie'";
    rcode(22) =  "res = t(res)}";

else
    rcode(12) = "coloc_res = coloc::coloc.abf(dataset1 = eqtl, " + ...
        "dataset2 = gwas,p1 = " + opts.prior(1) + ", p2 = " + ...
        opts.prior(2) + ", p12 = " + opts.prior(3) + ")";
    rcode(13) = "res <- coloc_res$summary";
    rcode(14) = "res <- c(res, c(method = 'abf'))";
end

rcode(23) = "write.table(as.data.frame(res), " + ...
    "row.names =  T, col.names = F, file = '" + name + ".txt', " + ...
    "quote = F, sep = '\t')";
writematrix(rcode', name + ".r", 'FileType', 'text', 'QuoteStrings', false)
MATLAB2Rconnector(name + ".r", 'delr', true);
out = readmatrix(name + ".txt", 'NumHeaderLines', 0, 'FileType', 'text', 'OutputType', 'string');

delete(name + ".txt");
delete(name + ".mat")

if isfile(name + ".r"); delete(name + ".r"); end
if isfile(name + ".r.Rout"); delete(name + ".r.Rout"); end
hdr = out(:, 1);
out = splitvars(table(out(:, 2:end).'));
% out = [varfun(@double, out, 'InputVariables', 1:width(out)-1), out(:, end)];
out.Properties.VariableNames = hdr;
if out.method == "abf"
    [out.hit1, out.hit2] = deal(strings(height(out), 1));
    out = movevars(out, {'hit1', 'hit2'}, 'After', 1);
    [out.idx1, out.idx2] = deal(nan(height(out), 1));
    out = movevars(out, {'idx1', 'idx2'}, 'Before', 'method');
end
out.gene_id = repmat(gene, height(out), 1);

end

%% ------------------------------------------------------------------------
function subset = getLDmat(gwas, opts)

if opts.verbose
    fprintf('\treading LD mat ...\n')
end

if opts.ldmethod == "internal"
    subset = gwas;
    subset.Properties.UserData = getbulkgeno(subset.snp, ...
        string(subset.chr), 'parallel', true, 'datatype', 'single', ...
        'verbose', false, 'chunk', round(height(subset)/70),...
        'home', opts.bgenhome, 'eid', opts.ldsample);
    fi = fieldnames(subset.Properties.UserData);
    subset.Properties.UserData = subset.Properties.UserData.(fi{1});
    
    [~, idx] = ismember(subset.Properties.UserData.snp, subset.snp);
    idx(idx < 1) = [];
    subset = subset(idx, :);

    % only keep samples used in GWAS (not needed anymore, 'eid' option was
    % added to getbulkgeno function)
    % subset.Properties.UserData.bed(~ismember(subset.Properties.UserData.eid, opts.ldsample), :) = [];

    % it's assumed that these are the bgen files used for GWAS, so we don't
    % compare the alternative allele to confirm they're the same variants
    % (especially because gwasall is already pruned).
    idx = subset.Properties.UserData.a2 ~= subset.allele1; % flip R sign
    if ~all(subset.Properties.UserData.a1(idx) == subset.allele1(idx))
        error('genotype data do not correspond to GWAS summary stat!')
    end

    subset.Properties.UserData.bed = normalize(subset.Properties.UserData.bed, 'zscore');
    subset.Properties.UserData = ...
        tril(subset.Properties.UserData.bed'*subset.Properties.UserData.bed./size(subset.Properties.UserData.bed, 1));
%     subset.Properties.UserData = tril(corr(subset.Properties.UserData.bed));
    
    subset.Properties.UserData(logical(eye(size(subset.Properties.UserData, 1)))) = 0.5;
    subset.Properties.UserData(idx, :) = -subset.Properties.UserData(idx, :);
    subset.Properties.UserData(:, idx) = -subset.Properties.UserData(:, idx);

    return
end

% external: uses precomputed LD structure from PolyFun paper
% https://stackoverflow.com/questions/4383571/importing-files-from-different-folder
% P = py.sys.path;
% P.insert(int32(1), opts.lddir);

if opts.liftover
    pos = "pos_old";
else
    pos = "pos";
end

gwas.id = gwas.chr + ":" + gwas.(pos) + ":" + gwas.allele0 + ":" + gwas.allele1;

% find LD files needed to be read
ldfiles = string({dir(fullfile(opts.lddir, '*.npz')).name}.');
indexfiles = regexprep(ldfiles, '.npz$', '.gz');
ldregion = split(regexprep(ldfiles, '.npz$', ''), '_');
ldregion(:, 1) = regexprep(ldregion(:, 1), '^chr', '');
ldregion = array2table(double(ldregion), 'VariableNames', {'chr', 'start', 'end'});
ldregion.index = fullfile(opts.lddir, indexfiles);
ldregion.files = fullfile(opts.lddir, ldfiles);
ldregion = sortrows(ldregion, 'start', 'ascend');
ldregion = sortrows(ldregion, 'chr', 'ascend');
ldregion(~ismember(ldregion.chr, gwas.chr), :) = [];

% fire-up python
warning('off', 'all')
terminate(pyenv)
spar = py.importlib.import_module("scipy.sparse");
np = py.importlib.import_module("numpy");

for i = 1:height(ldregion) % check/read LD files one by one 

    idx = gwas.chr == ldregion.chr(i) & ...
        gwas.(pos) >= ldregion.start(i) & ...
        gwas.(pos) <= ldregion.end(i);

    if ~any(idx); continue; end

    subset = gwas(idx, :);
    index = bfilereader(ldregion.index(i), parallel=true, header=true, verbose="off");
    index.id = index.chromosome + ":" + index.position + ":" + index.allele2 + ":"  + index.allele1;

    % intersect variants in subset gwas and varinats in the LD index file
    missIdx = ~ismember(subset.id, index.id);
    if any(missIdx) % then flip the alleles and recheck: in this case R sign must be reversed
        subset.id(missIdx) = subset.chr(missIdx) + ":" + ...
            subset.(pos)(missIdx) + ":" + subset.allele1(missIdx) +...
            ":" + subset.allele0(missIdx);
    end
    [idx1, idx2] = ismember(subset.id, index.id);
    subset.flipR = false(height(subset), 1);
    subset.flipR(idx1 & missIdx) = true;
    subset(~idx1, :) = [];
    idx2(idx2 < 1) = [];

    R = pyrun("R = spar.load_npz(npz).toarray()", "R", spar = spar, npz = ldregion.files(i));
    R = pyrun("R = R[np.meshgrid(tuple(n1),tuple(n2), sparse=True, indexing='ij')]",...
        "R", n1 = np.uint(idx2-1), n2 = np.uint(idx2-1), R = R, np = np);
    R = R.single;

    if any(subset.flipR)
        R(subset.flipR, subset.flipR) = -1.*R(subset.flipR, subset.flipR);
    end
    subset.id = [];
    subset.Properties.UserData = R;

    if sum(~idx) > 2 % there are still variants in GWAS that couldn't be found in single LD mat file
        error('fix here')
    else
        break
    end
end

if opts.verbose
    fprintf('\b\b Done.\n')
end

terminate(pyenv)
warning on
end

%% ------------------------------------------------------------------------
function out = prunecolocResult(res)
out = cell(numel(res), 1);
for i = 1:numel(res)
    coloc = cellfun(@(x) groupfilter(x, 'method', @(x) x >= 0.8, 'PP.H4.abf'), res{i}.coloc, 'uni', false);
    emptyIdx = cellfun(@isempty, coloc);
    
    if all(emptyIdx) % try a less stringent cutoff
        coloc = cellfun(@(x) groupfilter(x, 'method', @(x) x >= 0.5, 'PP.H4.abf'), res{i}.coloc, 'uni', false);
        emptyIdx = cellfun(@isempty, coloc);
    end
    
    coloc(emptyIdx) = [];
    res{i}(emptyIdx, :) = [];
    res{i}.ftp_path = []; res{i}.coloc = [];
    if isempty(coloc); continue; end

    % add gene names
    gene_id = vertcat(coloc{:});
    gene_id = unique(gene_id.gene_id);
    
    tab = table(gene_id);
    [tab.gene_name, tab.biotype] = deal(strings(height(tab), 1));
    tryMe = true;
    cnt = 1;
    while tryMe
        try
            genes = EnsemblREST(tab.gene_id, 'lookup', 'refGenome', '38');
            tryMe = false;
        catch exception
            pause(1) % wait for 1 sec!
            if (cnt > 6)
                throw(exception)
            end
            cnt = cnt + 1;
        end
    end
    
    [f1, f2] = ismember(tab.gene_id, genes.id); f2(f2<1) = [];
    tab.gene_name(f1) = genes.name(f2);
    tab.biotype(f1) = genes.biotype(f2);
    
    for j = 1:numel(coloc)
        [f1, f2] = ismember(coloc{j}.gene_id, tab.gene_id); f2(f2<1) = [];
        [coloc{j}.gene_name, coloc{j}.biotype] = deal(strings(height(coloc{j}), 1));
        coloc{j}.gene_name(f1) = tab.gene_name(f2);
        coloc{j}.biotype(f1) = tab.biotype(f2);

        idx = coloc{j}.biotype == "protein_coding"; % keep "protein_coding" only if available
        if any(idx)
            coloc{j}(~idx, :) = [];
        end

        tmp = res{i}(j, :);
        tmp = repmat(tmp, height(coloc{j}), 1);
        coloc{j} = [coloc{j}, tmp];
        
%         if all(coloc{j}.hit1 == "")
%             coloc{j}.hit1 = []; coloc{j}.hit2 = [];
%         end

%         if all(isnan(double(coloc{j}.idx1))) || all(coloc{j}.idx1 == "")
%             coloc{j}.idx1 = []; coloc{j}.idx2 = [];
%         end
    end

    out{i} = vertcat(coloc{:});
end

out = vertcat(out{:});

end

%% ------------------------------------------------------------------------
function tabix_paths = fetchEgeneFiles(tabix_paths, opts)
if opts.eqtl
    if ~isfield(opts, 'egenedir') % use a single cut-off 1E-6 for eQTL variants around the GWAS index
        tabix_paths.egene = cell(height(tabix_paths), 1);
        tabix_paths.egene = deal({1e-6});
       return
    end

    egenedir = struct2table(dir(opts.egenedir));
    egenedir = egenedir.name(egenedir.isdir);
    egenedir = string(egenedir);
    egenedir(ismember(egenedir, [".", ".."])) = [];
    studies = unique(tabix_paths.study);
    tabix_paths.egene = cell(height(tabix_paths), 1);
    exCols = ["gene_id", "gene_chr", "gene_name", "gene_start", "gene_end"]; % columns to be read from eGene files (GTEx format)

    if opts.verbose
        fprintf('reading eGene files...\n')
    end

    for i = 1:numel(studies)
        idx = find(tabix_paths.study == studies(i));

        if opts.verbose
            fprintf("\t(%d of %d)-study: %s", i, numel(studies), studies(i))
            pbar = progressGen(numel(idx));
        end

        if any(egenedir == studies(i)) % eGene files are present for this study
            efiles = getfilenames(fullfile(opts.egenedir, studies(i)), "gz", "fullpath", false).gz;
            for j = 1:numel(idx)
                if opts.verbose
                    progressGen(pbar, j);
                end
            
                idx2 = find(startsWith(efiles.lower, tabix_paths.qtl_group(idx(j)).lower));
                if isempty(idx2) % different naming?
                    idx2 = find(contains(efiles.lower, tabix_paths.qtl_group(idx(j)).lower));
                    if isempty(idx2) && strcmpi(tabix_paths.qtl_group(idx(j)).lower, "esophagus_gej") % manually set this dataset for GTEx
                        idx2 = find(startsWith(efiles.lower, "esophagus_gastroesophageal_junction"));
                    end

                    if isempty(idx2) && strcmpi(tabix_paths.qtl_group(idx(j)).lower, "lcl") % manually set this dataset for GTEx
                        idx2 = find(startsWith(efiles.lower, "cells_ebv-transformed_lymphocytes"));
                    end
                end
                if ~isempty(idx2)
                    if numel(idx2) > 1
                        if opts.verbose
                            fprintf('warning: found %d eGene files for %s, will pick the first!\n', numel(idx2), tabix_paths.qtl_group(idx(j)))
                        end
                        idx2 = idx2(1);
                    end
                    tabix_paths.egene{idx(j)} = bfilereader(fullfile(opts.egenedir, studies(i), efiles(idx2)), ...
                        "header", true, "filter", opts.fdrcutoff, ...
                        "filterCol", "qval", "operator", "<=", ...
                        "extractCol", exCols, "verbose", "off", "parallel", true);
                    tabix_paths.egene{idx(j)}.gene_chr = regexprep(tabix_paths.egene{idx(j)}.gene_chr, "^chr", "");
                    tabix_paths.egene{idx(j)}.gene_id = extractBefore(tabix_paths.egene{idx(j)}.gene_id, "."); % Gencode format

                    tabix_paths.egene{idx(j)} = renamevars(tabix_paths.egene{idx(j)}, exCols, erase(exCols, "gene_"));

                    % check if start is always < end (make strands the
                    % same)
                    if ~all(tabix_paths.egene{idx(j)}.start < tabix_paths.egene{idx(j)}.end)
                        tabix_paths.egene{idx(j)}.s = min(tabix_paths.egene{idx(j)}.start, tabix_paths.egene{idx(j)}.end);
                        tabix_paths.egene{idx(j)}.e = max(tabix_paths.egene{idx(j)}.start, tabix_paths.egene{idx(j)}.end);
                        tabix_paths.egene{idx(j)}(:, ["start", "end"]) = [];
                        tabix_paths.egene{idx(j)} = renamevars(tabix_paths.egene{idx(j)}, ["s", "e"], ["start", "end"]);
                    end
                else
                    tabix_paths.egene(idx(j)) = {1e-6};
                end
            end
        else % use default 1E-6 cut-off for cis-eQTL variants within the GWAS index region
            tabix_paths.egene(idx) = deal({1e-6});
        end
    end
end

if opts.verbose, fprintf('\n'); end

end 


% function meta = geteQTLmeta
% api = "https://www.ebi.ac.uk/eqtl/api/";
% studies = webread(api + "studies");
% studies = string({studies.x_embedded.studies.study_accession}.');
% 
% tissues = webread(api + "tissues?size=1000");
% tissues = tissues.x_embedded.tissues;
% tissues = struct2table(tissues, 'AsArray', true);
% tissues.tissue = string(tissues.tissue);
% for i = 1:numel(tissues.tissue)
%     tmp = webread(api + "tissues/"+ tissues.tissue(i) + "/studies?size=1000");
%     tmp = string({tmp.x_embedded.studies.study_accession}.');
% end
% end