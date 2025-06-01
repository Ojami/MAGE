function  [geno, wesinfo] = getWES(opts)
% getWES is mainly designed to be integrated in gwasrunner.m function for
% fetching WES UKBB variants. OUTPUTs
%   geno: a struct with fields bed/bim/fam containing wes calls
%   wesinfo: a table containing annotations for query gene/variants
% 
% Oveis Jamialahmadi, GU, March 2021.
% 
% [EDITED]: variants field can be set to "ALL" to get all variants on the
%           gene. In this case gnomAD API is used to get gene boundaries.
%           26 Apr 2021.
% [EDITED]: supports parallele computing for multiple genes per time.
%           27 Apr 2021
% [EDITED]: 'genesymbol' flag was added. true (default): gene
%           names/symbols are provided; false: Ensembl (ENSG) identifiers
%           in input query. 30 July 2021.
% [EDITED]: 'genePos' option was added. A 1X2 double vector containing gene
%           boundaries [start, end] for getting all variants on gene.
%           Effective if only 'variants' set to "ALL". 2 August 2021.
% [EDITED]: 'loftee' functionality was added to return only LoF variants
%           with a LOFTEE annotation of HC. 4 August 2021.
% [EDITED]: 'exclude' falg was added to remove varaints from analysis. 26
%           April 2022.
% [EDITED]: 'regenie' flag was added (default:false). (1), only returns
%           variant IDs to be sent to REGENIE for gene-based tests, and (2)
%           for multi-gene mode, also subsets genotype call files. 27 April
%           2022.
% [EDITED]: now supports multiple genes/variant sets and can merge calls
%           over multiple genes (multi-gene mode). If this latter option is
%           desired, 'gene' should be in format of "gene1,gene2,gene3" (,
%           or ; separate). 'variants' in this case should still be in the
%           same vector. variants in different cell elements will be
%           treated separately. 29 April 2022.
%           Example: 
%                   variants1={["v1","v2"], ["v3","v4","v5"]};
%                   variants2=["v6","v7"];
%           variants1 has two variant sets, each will be processed
%           distinctively, and it's assumed both sets belong to different
%           genes/regions. The same is true for variants2. But if in each
%           set, variants belong to different chromosomes (e.g. in
%           variants2 if "v6" and "v7" belong to different chromosomes),
%           then they merged together. So, variants within each set are
%           merged (either belong to the same or different chromosomes).
%           Note that 'gene' overrides 'variants' option. Therefore, if
%           only variants are desired, don't use 'gene' argument. If used,
%           first all annotated variants on the query gene are fetched,
%           then those matching 'variants' will be kept. Hence, only
%           variants present on the query 'gene' wil be filtered in the
%           end.
% [EDITED]: % 'ignore_anno' flag was added and if ture, ignores reading
%           annotation of variants in 'annfile'. Important note:
%           wesinfo table in this case cannot be used for annotation (as
%           used in gwasrunner 'vep' flag). 20 December 2022.
% [EDITED]: 'parallel' flag was added to read bed files using bedreader in
%           parallel (use only if there is a large number of variants). 
%           12 February 2024.

arguments
    opts.gene {mustBeText, mustBeVector} % if left empty, variant must be supplied
    opts.home {mustBeNonzeroLengthText} % home directory of WES data (BED/BIM/FAM files)
    opts.bed {mustBeTextScalar} = "" % [DEPRECATED] Q1 had only 1 BED file; but, Q4 has different BED files per each chromosome. So this field can be left empty, and function search for appropriate bed file based on chromosome
    opts.plinkdir {mustBeTextScalar} = "" % home directory of PLINK, if left empty opts.home will be used (PLINK1.9 must be in opts.home directory).
    opts.annhome {mustBeNonzeroLengthText} % home directory of WES annotated files (see runVEP.m function)
    opts.annfiles {mustBeNonzeroLengthText} % can be more than one file
    opts.lof (1,1) logical = false % only subset LoF variants
    opts.loftee (1,1) logical = false % only subset those LoF with a true LOFTEE flag.
    opts.variants {mustBeVector} % if left empty, gene must be supplied. For consistency, variants must be varIDs provived by UKBB in WES bim files.
    opts.genePos (:, 2) double % gene boundary. Works only if 'variants' set to "ALL" (all variants on the gene)
    opts.loader {mustBeTextScalar, mustBeMember(opts.loader, ["matlab", "bfilereader"])} = "matlab"
    opts.genesymbol (1,1) logical = true; % true: gene symbol, false: Ensembl (ENSG) identifier.
    opts.verbose (1,1) logical = true % prints related info
    opts.exclude {mustBeText, mustBeVector}
    opts.regenie (1,1) double {mustBeMember(opts.regenie, [0, 1, 2, 3])} = 0 % (1) only reads variant IDs, (2) also extract genotype files for multi-gene mode, (3) only return bed file names (if outputs of geneWrapper func are used)
    opts.outdir {mustBeFolder} = pwd % name of output directory for bed files
    opts.set_list {mustBeText} = "" % set-list file when 'regenie' is 3 (for regenie gene-based analysis)
    opts.anno_file {mustBeText} = "" % anno-file file when 'regenie' is 3 (for regenie gene-based analysis)
    opts.mask_def {mustBeText} = "" % mask file when 'regenie' is 3 (for regenie gene-based analysis)
    opts.ignore_anno (1,1) logical = false % ignores reading annotation of variants in 'annfile'. 
    opts.createBED (1,1) logical = false % to create PLINK files for the query gene/variant sets. Effective only when 'regenie' is used. 
    opts.parallel (1,1) logical = false % for bedreader
end

if ~ispc
    opts.loader = "matlab"; % current bfilereader class is slow on Linux
end

% check inputs: either 'gene' or 'variant' must be specified
if isfield(opts, 'gene')
    opts.geno = string(opts.gene);
    eidx = (opts.gene == "") | ismissing(opts.gene);
    opts.gene(eidx) = [];
else
    opts.gene = [];
end

if isfield(opts, 'variants') && ~isempty(opts.variants)
    % reformat variants (to diffrentiate between variant sets)
    if ~isempty(opts.variants) && ~any(cellfun(@any, cellfun(@(x) strcmp(x, ""), opts.variants, 'uni', false)))
        checkvars = opts.variants;
        if isstring(checkvars), checkvars = cellstr(checkvars); end
        if iscellstr(checkvars), checkvars = {checkvars}; end
        opts.variants = checkvars;
    end

    % e.g. opts.variants = {["snp1", "snp2"], ["snp3", "snp4", "snp5"]};
    if ~iscell(opts.variants) 
        eidx = (opts.variants == "") | ismissing(opts.variants);
        opts.variants(eidx) = [];
        if ~isempty(opts.variants)
            opts.variants = {opts.variants};
        else
            opts.variants = {};
        end
    end
else
    opts.variants = {};
end

if isempty(opts.gene) 
    if isempty(opts.variants)
        error("getWES:wrongInputs", "both variants and gene cannot be left empty!")
    elseif any(cellfun(@(x)any(strcmpi(x, "all")), opts.variants))
        error("getWES:wrongInputs", "when variants set to 'all', gene name should be present")
    end
end

if any(opts.set_list ~= "")
    if any(~isfile(opts.set_list))
        error('getWES:wrongInputs', "set_list must be a valid file!")
    end
    opts.regenie = 3;
    [geno, wesinfo] = readAnnoSetFiles(opts);
    return
end

if ismissing(opts.plinkdir) || opts.plinkdir == "" && opts.regenie ~= 3
    opts.plinkdir = opts.home; % PLINK1.9 or PLINK2 must be there!
end

% variant inclusion/exclusion criteria
if ~isfield(opts, 'exclude') || any(opts.exclude == ""); opts.exclude = []; end

% loop over genes
[geno, wesinfo] = deal({});

if isempty(opts.gene)
    opts.gene = repmat("", numel(opts.variants), 1);
end

for i = 1:numel(opts.gene)
    if opts.verbose
        fprintf('(%d of %d)-reading genotype calls for query %s\n', i, numel(opts.gene), opts.gene(i))
    end

    inopts = opts;
    inopts.gene = opts.gene(i);
    inopts.gene = split(inopts.gene, [",", ";"]); % multi-gene mode: genes to be merged

    % check gene/variants query
    if isempty(inopts.gene)
        inopts.gene = "";
    end
    
    inopts.getWholeGene = false;
    if isempty(inopts.variants)
        inopts.variants = "";
    else
        inopts.variants = string(opts.variants{i});
        if all(inopts.variants.lower == "all") 
            inopts.getWholeGene = true; % fetches all variants on the gene.
            inopts.variants = "";
        end
    end
    
    if any(inopts.gene == "") && (all(inopts.variants == "") || all(inopts.variants == "ALL"))
        error("both variants and gene cannot be left empty!")
    end

    if numel(inopts.gene) > 1 && (all(inopts.gene ~= "") || ~isempty(inopts.gene))
        if inopts.regenie
            inopts.regenie = 2;
        end
    elseif opts.createBED && inopts.regenie
        inopts.regenie = 2;
    end

    % check if the variants in set i belong to the same chromosome, if not,
    % they should be merged also
    if all(inopts.variants ~= "")
        chr = extractBefore(inopts.variants, regexpPattern(":|_|-"));
        chrU = unique(chr);
        if numel(chrU) > 1 && inopts.regenie
            inopts.regenie = 2;
            inopts.gene = repmat("", numel(chrU), 1); % fake gene names
            varlist = inopts.variants;
            inopts.variants = cell(numel(chrU), 1);
            for k = 1:numel(chrU)
                inopts.variants{k} = varlist(chr == chrU(k));
            end
        elseif opts.createBED && inopts.regenie
            inopts.regenie = 2; % create a bed file even for single genes
        end

    end

    
    [genot, wesinfot] = deal(cell(numel(inopts.gene), 1));
    for j = 1:numel(inopts.gene)
        inoptsj = inopts;
        inoptsj.gene = inopts.gene(j);
        if inopts.gene(j) == "" && numel(inopts.gene) > 1 % merge variant set
            inoptsj.variants = inopts.variants{j};
        end

        if ~isfield(inoptsj, 'genePos'), inoptsj.genePos = nan(1, 2); end
        [genot{j}, wesinfot{j}] = getWESsingle(inoptsj);
    end

    [geno{i}, wesinfo{i}] = geneCallMerger(genot, wesinfot, inopts);
    clear genot wesinfot
    
    if opts.verbose; fprintf('\n'); end
end

end % END

%% subfunctions ===========================================================
function [callData, genedata] = getWESsingle(opts)

if opts.regenie ~= 1 % for 'regenie' single-gene only read bim file
    if ~isfile(fullfile(opts.plinkdir, 'plink.exe')) % should be also implemented for linux
        if ~isfile(fullfile(opts.plinkdir, 'plink2.exe'))
            error('couldn''t find plink within %s!', opts.plinkdir)
        else
            opts.plinkdir = fullfile(opts.plinkdir, 'plink2.exe');
        end
    else
        opts.plinkdir = fullfile(opts.plinkdir, 'plink.exe');
    end
end

if opts.getWholeGene % gets gene boundaries
    try % try Ensembl REST
        restdata = EnsemblREST(opts.gene, 'lookup', 'geneSymbol', opts.genesymbol, 'refGenome', '38');
        geneInfo.pos = [restdata.start, restdata.end];
        geneInfo.chr = restdata.chr;
    catch % try gnomAD
        if opts.genesymbol
            gnomadterm = 'gene_symbol';
        else
            gnomadterm = 'gene_id';
        end
        gnomadData = gnomad(gnomadterm, opts.gene, 'dataset', 'gnomad_r3', 'reference_genome', 'GRCh38');
        geneInfo.pos = [gnomadData.data.gene.start, gnomadData.data.gene.stop];
        geneInfo.chr = string(gnomadData.data.gene.chrom);
    end

    if all(~isnan(opts.genePos))
        geneInfo.pos = opts.genePos;
    end
end

[bed, bim, fam, genedata] = deal({});
if opts.verbose
    fprintf('fetching variant IDs from annotation files...\n')
end
for ii = 1:numel(opts.annfiles)
    [~, wesfile] = fileparts(fullfile(opts.annhome, opts.annfiles(ii)));

    if opts.ignore_anno
        queryinfo = [];
    else
        extractCols = ["Uploaded_variation", "Location", "Gene", ...
            "Consequence", "Existing_variation", "SYMBOL", "Amino_acids"];
        wes_ann_file = fullfile(opts.annhome, wesfile + ".txt");

        if opts.gene ~= ""
            if strcmp(opts.loader, "bfilereader")
                patt = "(?<!-)\b" + opts.gene + "(?!-)\b";
                if opts.genesymbol
                    pattCol = "SYMBOL";
                else
                    pattCol = "Gene";
                end
                   
                queryinfo = bfilereader(wes_ann_file, ...
                    'header', true, 'pattern', patt, 'patternCol', pattCol,...
                    'verbose', 'off', 'parallel', true, ...
                    'extractCol', extractCols);

                if ~isempty(queryinfo)
                    if opts.loftee
                        fprintf('WARNING:LOFTEE not implemented')
                        opts.loftee = false;
                    else
                        queryCols = {'id', 'pos', 'ENSG', 'consequence',...
                            'rsid', 'gene', 'aa.change'};
                    end
                    queryinfo.Properties.VariableNames = queryCols;
                end
            else
                % note: my benchmark showed readtable is faster than load for
                % network drive fetch, though it's a bit slower for local
                % access. Therefore, I decided to go with network mode
                fileOpts = detectImportOptions(fullfile(opts.annhome,wesfile + ".txt"),...
                    'PreserveVariableNames', true, 'TextType', 'string', ...
                    'ReadVariableNames', true);
                if opts.loftee
                    fprintf('WARNING:LOFTEE not implemented')
                    opts.loftee = false;
                end
                
                if opts.genesymbol
                    fileOpts.SelectedVariableNames = {'SYMBOL'};
                else
                    fileOpts.SelectedVariableNames = {'Gene'};
                end
                wes_anndata = readtable(fullfile(opts.annhome,wesfile + ".txt"), fileOpts);
                % wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)),'gene'); % see note above
                idx1 = ismember(lower(wes_anndata.(1)), lower(opts.gene));
                idx2 = contains(wes_anndata.(1), ";"+lower(opts.gene)+";", 'IgnoreCase',true);
                idx3 = startsWith(wes_anndata.(1), lower(opts.gene)+";", 'IgnoreCase', true);
                idx4 = endsWith(wes_anndata.(1), ";"+lower(opts.gene), 'IgnoreCase', true);
                query_idx = idx1 | idx2 | idx3 | idx4;
                clear idx1 idx2 idx3 idx4
            end
        elseif ~isempty(opts.variants)
            if strcmp(opts.loader, "bfilereader")
                queryinfo = bfilereader(wes_ann_file, ...
                        'header', true, 'pattern', opts.variants, ...
                        'patternCol', "Uploaded_variation",...
                        'verbose', 'off', 'parallel', true, ...
                        'extractCol', extractCols);
                if ~isempty(queryinfo)
                    if opts.loftee
                        fprintf('WARNING:LOFTEE not implemented')
                        opts.loftee = false;
                    else
                        queryCols = {'id', 'pos', 'ENSG', 'consequence',...
                            'rsid', 'gene', 'aa.change'};
                    end
                    queryinfo.Properties.VariableNames = queryCols;
                end
            else
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'Uploaded_variation');
                query_idx = ismember(wes_anndata.Uploaded_variation, opts.variants);
            end
        else
            error('No input for WES data!')
        end

        if ~strcmp(opts.loader, "bfilereader") && exist(fullfile(opts.annhome,wesfile + '.txt'), 'file') % Text file
            % Large MAT files are not efficient to load. Instead we read an equivalent txt file 
            if any(query_idx)

                query_idx = find(query_idx) + 1; % 1 for header
                fileOpts = detectImportOptions(fullfile(opts.annhome,wesfile + '.txt'),...
                    'PreserveVariableNames', true, 'TextType', 'string', ...
                    'ReadVariableNames', true);
                fileOpts.SelectedVariableNames = extractCols;
                if opts.loftee
                    fileOpts.SelectedVariableNames = ...
                        {fileOpts.VariableNames, 'LOFTEE'};
                end
                linePosDiff = find(diff(query_idx) > 1);
                C = (0);
                if ~isempty(linePosDiff)
                    for i = 1:numel(linePosDiff)
                        if i == 1
                            C(i, 1:2) = [query_idx(1), query_idx(linePosDiff(i))];
                        else
                            C(i, 1:2) = [query_idx(linePosDiff(i - 1) + 1),...
                                query_idx(linePosDiff(i))];
                        end
                    end
                    C(i + 1, 1:2) = [query_idx(linePosDiff(i) + 1), query_idx(end)];
                else
                    C = [query_idx(1), query_idx(end)];
                end
                fileOpts.DataLines = C;
                queryinfo = readtable(fullfile(opts.annhome,wesfile + '.txt'), fileOpts);
                queryinfo.Properties.VariableNames = {'id', 'pos', 'ENSG',...
                        'consequence', 'rsid', 'gene', 'aa.change'};
                
                % check if text file matches annotation MAT file
                if opts.genesymbol
                    geneCol = "gene";
                else
                    geneCol = "ENSG";
                end
                if ~all(contains(lower(queryinfo.(geneCol)), lower(opts.gene)))
                    error('mismatched annotation text/MAT files! Try to generate a new text file!')
                end
            else
                queryinfo = [];
            end

        elseif ~strcmp(opts.loader, "bfilereader") % MAT file
            if any(query_idx)
                queryinfo = strings(sum(query_idx), 7);
                queryinfo(:, 4) = wes_anndata.gene(query_idx);
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'geneid');
                queryinfo(:, 5) = wes_anndata.geneid(query_idx);
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'id');
                queryinfo(:, 1) = wes_anndata.id(query_idx);
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'dbsnp_id');
                queryinfo(:, 2) = wes_anndata.dbsnp_id(query_idx);
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'pos');
                queryinfo(:, 3) = wes_anndata.pos(query_idx);
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'effect');
                queryinfo(:, 6) = wes_anndata.effect(query_idx);
                wes_anndata = load(fullfile(opts.annhome, opts.annfiles(ii)), 'aa_change');
                queryinfo(:, 7) = wes_anndata.aa_change(query_idx);
            else
                queryinfo = [];
            end
        end
    end
    
    if isempty(queryinfo) && ~opts.ignore_anno
        % comment since there may be unannotated variants in the query list
        if any(opts.variants ~= "")
            fprintf('WARNING: input variants are not annotated!\n')
        else % nothing found for this gene
            fprintf('oops! query was not found in WES Q4!\n')
            [callData.bed, callData.bim, callData.fam] = deal([]);
            genedata = table;
            return
        end
    end
    
    if opts.loftee
        queryCols = {'id', 'rsid', 'pos', 'gene', 'ENSG', 'consequence',...
            'aa.change', 'LOFTEE'};
    else
        queryCols = {'id', 'rsid', 'pos', 'gene', 'ENSG', 'consequence',...
            'aa.change'};
    end

    if ~istable(queryinfo) && ~isempty(queryinfo)
        queryinfo = array2table(queryinfo, 'VariableNames', queryCols);
    elseif isempty(queryinfo)
        queryinfo = table("", "", "", "", "", "", "",...
            'VariableNames', queryCols);
    end
    
    if opts.lof
        % LOF filters: https://www.nature.com/articles/s41586-020-2853-0
        % A more stringent definition does not contain start_/stop_lost as
        % the one used in LOFTEE
        % In UKBB WES 450 K paper, start-/stop-lost are also in pLoF
        % variants: https://www.nature.com/articles/s41586-021-04103-z
        disp("@06NOV2024: start-/stop-lost are now considered in LoF: https://www.nature.com/articles/s41586-021-04103-z")
        lof_filters = ["stop_gained"; "splice_donor_variant"; ...
            "frameshift_variant"; "splice_acceptor_variant"; ...
            "start_lost"; "stop_lost"];
    else
        lof_filters = ["stop_gained"; "start_lost"; "splice_donor_variant"; ...
            "frameshift_variant"; "splice_acceptor_variant"; "stop_lost";...
            "transcript_ablation"; "missense_variant"; ...
            "protein_altering_variant"; "inframe_insertion"; "inframe_deletion"];
    end
    
    % check if queryinfo is empty. if empty, it means query variant are
    % unannotated
    if size(queryinfo, 1) ~= 1 && queryinfo.id(1) ~= ""
        lof_idx = contains(queryinfo.consequence, lof_filters);
        queryinfo = queryinfo(lof_idx, :);
        
        if opts.loftee % only applies to LoF variants
            if opts.lof
                queryinfo(queryinfo.LOFTEE ~= "true", :) = [];
            else 
                lof_filters = ["stop_gained"; "splice_donor_variant"; ...
                    "frameshift_variant"; "splice_acceptor_variant"];
                lof_idx = contains(queryinfo.consequence, lof_filters);
                loftee_idx = queryinfo.LOFTEE ~= "true";
                queryinfo(loftee_idx & lof_idx, :) = [];
            end
        end
    end
        
    % check if there are unannotated variants in opts.variants
    if any(opts.variants ~= "")
        leftoverIDs = setdiff(opts.variants, queryinfo.id);
        warning('off', 'MATLAB:table:RowsAddedExistingVars');
        queryinfo.id(end+1:end+numel(leftoverIDs)) = leftoverIDs;
        warning('on', 'MATLAB:table:RowsAddedExistingVars');
        queryinfo = fillmissing(queryinfo, 'constant', "-", 'DataVariables', @isstring);

        % keep only custom set of variants in 'variants' field
        queryinfo(~ismember(queryinfo.id, opts.variants), :) = [];
    end
    
    % remove empty rows, this happens if query variants are unannotated.
    queryinfo(queryinfo.id == "", :) = [];

    % apply variant inclusion/exclusion criteria
    if ~isempty(opts.exclude) && ~isempty(queryinfo)
        idx = ismember(queryinfo.id, opts.exclude);
        queryinfo(idx, :) = [];
        if any(idx) && opts.verbose
            fprintf('%d variant were removed due to exclusion criteria\n', sum(idx))
        end
    end
    
    if opts.verbose
        if ~isempty(queryinfo)
            summarycons = groupsummary(queryinfo, 'consequence');
            summarycons.Properties.VariableNames{2} = 'Count';
            summarycons = sortrows(summarycons, 'Count', 'descend');
            summarycons.(2) = string(summarycons.(2));
            summarycons(end+1, :) = {'<strong>Total</strong>', ...
                "<strong>" + sum(double(summarycons.Count)) + "</strong>"};
            fprintf('summary of consequences for query gene %s:\n', join(unique(queryinfo.gene),'|'))
            summarycons = convertvars(summarycons, 1:2, 'categorical');
            disp(summarycons)
        end
    end
    
    % it's quite unlikely that a gene only contains unannotated variants.
    % in this case queryinfo is empty; so to avoid error/skippig variant
    % call read, getWholeGene must be false. Because the only possible
    % (though unlikely) scenario that queryinfo is empty is getWholeGene to
    % be true
    if ~isempty(queryinfo) || opts.getWholeGene 
        
        % generate random strings for parallel computing
        randStr = getRandomName("t", 10);

        % Write IDs to a file for plink
        if ~isempty(queryinfo)
            queryIDs = cellstr(queryinfo.id);
            if ~opts.getWholeGene
                writecell(queryIDs, fullfile(opts.outdir, "tmpID."+randStr+".txt"))
            end
        else
            queryIDs = {''};
        end
       
        if isempty(opts.bed(ii)) || opts.bed(ii) == ""
            % Chromosome number can be found from query_gene_id_str, and
            % therefore opts.bed can be set.
            if opts.getWholeGene && isempty(queryinfo)
                chr = string(geneInfo.chr);
            else
                chr = string(regexp(queryIDs(1, 1), '^(.+?)[^:]*',...
                    'match'));
            end

            if double(chr) > 22 || any(chr.lower == ["x", "y"]) % sex chr
                if isfolder(fullfile(opts.home, "sex.chrom/"))
                    opts.sexchr = true;
                    opts.home = fullfile(opts.home, "sex.chrom/");
                    chr = replace(chr, ["23","24"], ["X", "Y"]);
                end
            else
                opts.sexchr = false;
            end
            
            existingBEDfiles = getfilenames(opts.home, "bed").bed;
            if isempty(existingBEDfiles)
                error('no bed file found within opts.home!')
            end
            if numel(existingBEDfiles) < 2 % only found 1 BED file in WES_homeidr
                opts.bed(ii) = existingBEDfiles{1};
            else % multiple BED files 
                if opts.sexchr
                    chrPos = strshare(existingBEDfiles, 'pattern', true, 'patternstr', "(?<=_c)(.*?)(?=_)");
                    chrPos = replace(chrPos, "%d", "%s");
                    opts.bed(ii) = string(sprintf(chrPos, chr));
                else
                    chrPos = strshare(existingBEDfiles, 'pattern', true);
                    opts.bed(ii) = string(sprintf(chrPos, double(chr)));
                end        
                
            end
            opts.bed(ii) = string(regexprep(opts.bed(ii), '.bed$', ''));
        end
        
        if opts.getWholeGene % fetch all variants available on the query gene
            tbim = bimreader(fullfile(opts.home, opts.bed(ii) + ".bim"), ...
                'cols', ["snp", "pos"]);
            tbim = tbim.snp(tbim.pos>=geneInfo.pos(1) & tbim.pos<=geneInfo.pos(2));

            % apply variant inclusion/exclusion criteria
            if ~isempty(opts.exclude)
                idx = ismember(tbim, opts.exclude);
                tbim(idx) = [];
                if any(idx) && opts.verbose
                    fprintf('%d variant were removed due to exclusion criteria\n', sum(idx))
                end
            end

            tbim = setdiff(tbim, queryIDs);
            if opts.verbose
                fprintf('%d more unannotated variants have been added.\n', numel(tbim))
            end

            queryIDs = union(queryIDs, tbim);
            writematrix(string(queryIDs), ...
                fullfile(opts.outdir, "tmpID." + randStr + ".txt"), ...
                'QuoteStrings', false, 'FileType', 'text')
            
        end
        
        if opts.regenie == 1 % don't read genotype calls for single-gene
            delete(fullfile(opts.outdir, "tmpID." + randStr + ".txt"))
            tbim = bimreader(fullfile(opts.home, opts.bed(ii) + ".bim"), 'struct', false);
            tbim(~ismember(tbim.snp, queryIDs), :) = [];
            bim{ii} = table2struct(tbim, 'ToScalar', true);
        else
            pcmd = '"' + opts.plinkdir + '"' + " --bfile " + '"' + ...
            fullfile(opts.home, opts.bed(ii)) + '"' + " --extract " + ...
                fullfile(opts.outdir, "tmpID." + randStr + ".txt") + ...
                " --make-bed --out " + fullfile(opts.outdir, "tmpBED" + ii + randStr);
            if opts.verbose
                fprintf('Running plink for annotated files %d of %d\n', ii, numel(opts.annfiles))
            end
            if opts.sexchr, pcmd = replace(pcmd, "plink.exe", "plink2.exe"); end
            [~, ~] = system(pcmd);
            warning off
            delete(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".nosex"))
            delete(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".log"))
            delete(fullfile(opts.outdir, "tmpID." + randStr + ".txt"))
            warning on
            
            if ~opts.regenie
                % Extract variant calls using bedreader
                if ~isempty(opts.variants) && all(opts.variants ~= "")
                    queryinfo_tmp = queryinfo(ismember(queryinfo.id, opts.variants), :);
                    if isempty(queryinfo_tmp)
                        queryinfo = queryinfo(contains(queryinfo.rsid, opts.variants), :);
                        opts.variants = queryinfo.id;
                    else
                        queryinfo = queryinfo_tmp;
                    end
                    [bed{ii}, bim{ii}, fam{ii}] = bedreader(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".bed"), opts.variants, parallel=opts.parallel);
                else
                    [bed{ii}, bim{ii}, fam{ii}] = bedreader(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".bed"), queryIDs, parallel=opts.parallel);
                end
                warning off
                delete(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".bed"))
                delete(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".bim"))
                delete(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".fam"))
                warning on
            else % regenie == 2 for multi-gene mode
                bim{ii} = bimreader(fullfile(opts.outdir, "tmpBED" + ii + randStr + ".bim"), 'struct', true);
                opts.bfile(ii) = fullfile(opts.outdir, "tmpBED" + ii + randStr);
            end
        end % end of opts.regenie if
        
    else
        if opts.verbose
            fprintf('Gene name was not found in file %d\n', ii)
        end
    end
    genedata{ii} = queryinfo;
end

empty_idx = cellfun(@isempty, genedata);
genedata(empty_idx) = [];
genedata = vertcat(genedata{:});
bim(empty_idx) = [];
tmpbim.chr = cellfun(@(x)x.chr, bim, 'uni', false); tmpbim.chr = vertcat(tmpbim.chr{:});
tmpbim.snp = cellfun(@(x)x.snp, bim, 'uni', false); tmpbim.snp = vertcat(tmpbim.snp{:});
tmpbim.pos = cellfun(@(x)x.pos, bim, 'uni', false); tmpbim.pos = vertcat(tmpbim.pos{:});
tmpbim.a1 = cellfun(@(x)x.a1, bim, 'uni', false); tmpbim.a1 = vertcat(tmpbim.a1{:});
tmpbim.a2 = cellfun(@(x)x.a2, bim, 'uni', false); tmpbim.a2 = vertcat(tmpbim.a2{:});
callData.bim = tmpbim;

if ~opts.regenie
    bed(empty_idx) = [];
    fam(empty_idx) = [];
    callData.bed = horzcat(bed{:});
    callData.fam = fam{1};
elseif opts.regenie == 2 % to be merged for multi-gene mode
    callData.bfile = opts.bfile;
end

end % END

%% ------------------------------------------------------------------------
function [geno, wesinfo] = geneCallMerger(geno, wesinfo, opts)

if opts.regenie ~= 1 % for 'regenie' single-gene only read bim file
    if ~isfile(fullfile(opts.plinkdir, 'plink.exe')) % should be also implemented for linux
        if ~isfile(fullfile(opts.plinkdir, 'plink2.exe'))
            error('couldn''t find plink within %s!', opts.plinkdir)
        else
            opts.plinkdir = fullfile(opts.plinkdir, 'plink2.exe');
        end
    else
        opts.plinkdir = fullfile(opts.plinkdir, 'plink.exe');
    end
end

% remove empty returned data (maybe wrong  gene symbols or absence of
% data in WES data)
emptIdx = cellfun(@isempty, wesinfo);
geno(emptIdx) = []; wesinfo(emptIdx) = [];

if numel(geno) > 1
    wesinfo = vertcat(wesinfo{:});

    if opts.regenie
        [tmpgeno.bim.chr, tmpgeno.bim.snp,...
            tmpgeno.bim.pos, tmpgeno.bim.a1, tmpgeno.bim.a2] = deal([]);
        for i = 1:numel(geno)
            tmpgeno.bim.chr = [tmpgeno.bim.chr; geno{i}.bim.chr];
            tmpgeno.bim.snp = [tmpgeno.bim.snp; geno{i}.bim.snp];
            tmpgeno.bim.pos = [tmpgeno.bim.pos; geno{i}.bim.pos];
            tmpgeno.bim.a1 = [tmpgeno.bim.a1; geno{i}.bim.a1];
            tmpgeno.bim.a2 = [tmpgeno.bim.a2; geno{i}.bim.a2];
        end

        % name of genotype files to be merged
        if opts.regenie == 2

            patt = getRandomName("t", 5);
            mergeFiles =  cellfun(@(x) x.bfile, geno);
            if ~iscolumn(mergeFiles); mergeFiles = mergeFiles'; end
            writematrix(mergeFiles, fullfile(opts.outdir, patt + ".plink.merge"), ...
                'FileType', 'text', 'QuoteStrings', false)
            
            % call PLINK
            [~, ~] = system('"' + opts.plinkdir + '"' + " --merge-list " +...
                '"' + fullfile(opts.outdir, patt + ".plink.merge") + '"' + ...
                " --keep-allele-order --allow-no-sex --allow-extra-chr --make-bed --out " + ...
                fullfile(opts.outdir, patt + ".merged"));
        
            % remove unnecessary files
            mergeFiles = mergeFiles + [".bed", ".bim", ".fam"];
            rmfiles = [fullfile(opts.outdir, patt + ".plink.merge"); mergeFiles(:); ...
                fullfile(opts.outdir, patt + ".merged." + ["nosex"; "log"])];
            delete(rmfiles{:})
        
            tmpgeno.bfile = fullfile(opts.outdir, patt + ".merged.bed");
        end

    else
        fam = cellfun(@(x)x.fam, geno, 'uni', false);
        tmpgeno.fam = fam{1};
        for i = 2:numel(geno)
            tmpgeno.fam = intersect(tmpgeno.fam, fam{i});
        end
        [tmpgeno.bed, tmpgeno.bim.chr, tmpgeno.bim.snp,...
            tmpgeno.bim.pos, tmpgeno.bim.a1, tmpgeno.bim.a2] = deal([]);
        for i = 1:numel(geno)
            [~, idx] = ismember(tmpgeno.fam, geno{i}.fam); idx(idx<1) = [];
            tmpgeno.bed = [tmpgeno.bed, geno{i}.bed(idx, :)];
            tmpgeno.bim.chr = [tmpgeno.bim.chr; geno{i}.bim.chr];
            tmpgeno.bim.snp = [tmpgeno.bim.snp; geno{i}.bim.snp];
            tmpgeno.bim.pos = [tmpgeno.bim.pos; geno{i}.bim.pos];
            tmpgeno.bim.a1 = [tmpgeno.bim.a1; geno{i}.bim.a1];
            tmpgeno.bim.a2 = [tmpgeno.bim.a2; geno{i}.bim.a2];
        end
    end
    geno = tmpgeno; clear tmpgeno
elseif isscalar(wesinfo)
    geno = geno{1}; wesinfo = wesinfo{1};
else % non of query genes/variants found in WES data
    error('gwasrunner:getWES', 'none of queries were found in WES data!')
end

end % END

%% ------------------------------------------------------------------------
function [geno, wesinfo] = readAnnoSetFiles(opts)

if ~isempty(opts.variants)
    if iscell(opts.variants), opts.variants = vertcat(opts.variants{:}); end
    opts.variants = unique(string(opts.variants));
end

opts.merger = false; % merge masks/anno/set files? only if they're on the same chrom and ";" or "," separated
if all(opts.gene ~= "")
    col = 1;
    patt = opts.gene;
    
    idx = contains(patt, [";", ","]);
    if any(idx) % merge genes (works only if they're on the same chrom)
        % which which gene sets to merge
        opts.merger_patt = patt(idx);
        patt = arrayfun(@(x)split(x, [";", ","]), patt, uni=false);
        patt = unique(vertcat(patt{:}));
        opts.merger = true;
    end
else % variants
    col = 4;
    patt = opts.variants;
end

if opts.loader == "bfilereader"
    tab.set = bfilereader(opts.set_list, "pattern", patt, ...
        "patternCol", col, "parallel", true, "verbose", "off", ...
        "header", false, "sep", char(9));
else % matlab
    [~, ~, ext] = fileparts(opts.set_list);
    tab.set = readall(tabularTextDatastore(opts.set_list, "ReadVariableNames", false, ...
        "TextType", "string", "FileExtensions", ext, "Delimiter", char(9)));
    tab.set(~contains(tab.set.(col), patt), :) = [];
end

if isempty(tab.set) % no info for this query
    if opts.verbose
        fprintf('getWES:query gene/variants cannot be found in the input set files!\n')
    end
     geno = []; wesinfo = table;
    return
end

% match exact gene names
% note that some genes maybe missing, we don't care about them!
frem = false(height(tab.set), 1);
if col == 1 && height(tab.set) > numel(patt) 
    for i = 1:height(tab.set)
        tmpgene = split(tab.set.(1)(i), ["(", ")", "_"]);
        tmpgene(tmpgene == "") = [];
        if ~any(ismember(tmpgene, patt))
            frem(i) = true;
        end
    end
end
tab.set(frem, :) = [];

ids = cell(height(tab.set), 1);
if ~isempty(opts.variants) || (isfield(opts, 'exclude') && all(opts.exclude ~= "")) % filter variants
    for i = 1:height(tab.set)
        ids{i} = split(tab.set.(4)(i), ",");
        if ~isempty(opts.variants)
            idx = ismember(ids{i}, opts.variants);
        end

        if all(opts.exclude ~= "")
            idx = idx & ~ismember(ids{i}, opts.exclude);
        end

        ids{i} = ids{i}(idx);
        tab.set.(4)(i) = join(ids{i}, ",");
    end
end

% tailor anno_file
if opts.loader == "bfilereader"
    apatt = tab.set.(1);
    apatt = replace(apatt, ["(", ")"], ["\(", "\)"]);
    tab.ann = bfilereader(opts.anno_file, "pattern", apatt, ...
        "patternCol", 2, "parallel", true, "verbose", "off", ...
        "header", false, "sep", char(9));
else % matlab
    [~, ~, ext] = fileparts(opts.set_list);
    tab.ann = readall(tabularTextDatastore(opts.anno_file, "ReadVariableNames", false, ...
        "TextType", "string", "FileExtensions", ext, "Delimiter", char(9)));
    
end
tab.ann(~ismember(tab.ann.(2), tab.set.(1)), :) = [];

idx = false(height(tab.ann), height(tab.set));
for i = 1:numel(ids)
    if ~isempty(ids{i})
        idx(:, i) = ~ismember(tab.ann.(1), ids{i}) & ismember(tab.ann.(2), tab.set.(1)(i));
    end
end
idx = any(idx, 2);
tab.ann(idx, :) = [];

patt = getRandomName("gene", 5);
geno.anno_file = fullfile(opts.outdir, patt + ".ann.txt");
geno.set_list = fullfile(opts.outdir, patt + ".set.txt");

% check if there are genes to be merged into one single region
if opts.merger
    % append ann and set tables: not that for merger to work, input must be
    % gene symbols, since the format in set file is name(id). See
    % geneWrapper docs.
    [set_merger, ann_merger] = deal(cell(numel(opts.merger_patt), 1));
    for k = 1:numel(opts.merger_patt)
        this_region = split(opts.merger_patt(k), [";", ","]);
        tmp_tab = tab.set.(1).startsWith(this_region + ("("|textBoundary("end")));
        
        if any(tmp_tab)
            tmp_tab = tab.set(tmp_tab, :);
            set_merger{k} = tmp_tab(1, :);
            set_merger{k}.(3) = mean(tmp_tab.(3)); % position
            set_merger{k}.(1) = replace(join(tmp_tab.(1), "_"), ["(", ")"], "_");
            tmp_vars = tmp_tab.(4).join(",");
            set_merger{k}.(4) = join(unique(tmp_vars.split(",")), ",");
            
            clear tmp_tab
            tmp_tab = tab.ann.(2).startsWith(this_region + ("("|textBoundary("end")));
            if any(tmp_tab)
                ann_merger{k} = tab.ann(tmp_tab, :);
                ann_merger{k}.(2)(:) = set_merger{k}.(1)(1);
            end
        end
    end
    set_merger(cellfun(@isempty, set_merger)) = [];
    ann_merger(cellfun(@isempty, ann_merger)) = [];
    tab.set = [tab.set; vertcat(set_merger{:})];
    tab.ann = [tab.ann; vertcat(ann_merger{:})];

end
writetable(tab.ann, geno.anno_file, "Delimiter", char(9), "WriteVariableNames", false);
writetable(tab.set, geno.set_list, "Delimiter", char(9), "WriteVariableNames", false);

tab.ann.Properties.VariableNames(1:3) = ["id", "gene", "consequence"];
% extract gene id (Ensembl)
idx = contains(tab.ann.gene, "(") & ~contains(tab.ann.gene, ","); % second logic for merger
tab.ann.ENSG = strings(height(tab.ann), 1);
tab.ann.ENSG(idx) = extractBetween(tab.ann.gene(idx), "(", ")");
tab.ann.gene(idx) = extractBefore(tab.ann.gene(idx), "(");

geno.bim.snp = tab.ann.id;
try
    tab.ann.pos = double(extractBetween(tab.ann.id, tab.set.(2)(1) + ":", ":"));
catch
    tab.ann.pos = nan(height(tab.ann), 1);
end

% tailor mask_def file
thesecats = unique(tab.ann.consequence);
maskfile = readlines(opts.mask_def);
maskfile(maskfile == "") = [];
maskfile = split(maskfile);
for i = 1:size(maskfile, 1)
    tmpcats = split(maskfile(i, 2), ",");
    tmpcats = intersect(tmpcats, thesecats);
    maskfile(i, 2) = join(tmpcats, ",");
end
[~, uidx] = unique(maskfile(:, 2), "stable");
maskfile = maskfile(uidx, :);
geno.mask_def = fullfile(opts.outdir, patt + ".mask.txt");
writematrix(join(maskfile, char(9)), geno.mask_def , "QuoteStrings", false);

% find genotype file
bfiles = getfilenames(opts.home, "bed").bed;
if isempty(bfiles)
    bfiles = getfilenames(opts.hom, "bgen").bgen;
    if isempty(bfiles)
        bfiles = getfilenames(opts.hom, "pgen").pgen;
    end
end

if any(contains(lower(bfiles), ["_cx_", "_cy_", "_cxy_"])) % sex chromosomes
    chrPatt = strshare(bfiles, 'pattern', true, 'patternstr', "(?<=_c)(.*?)(?=_)");
else
    chrPatt = strshare(bfiles, 'pattern', true);
end
chrPatt = replace(chrPatt, "%d", "%s"); 

chr = string(tab.set.(2)(1));
genofile = fullfile(opts.home, compose(chrPatt, chr));
if chr == "23" && ~isfile(genofile)
    genofile = fullfile(opts.home, compose(chrPatt, "X"));
    if ~isfile(genofile)
        genofile = fullfile(opts.home, compose(chrPatt, "x"));
    end
end

if chr == "24" && ~isfile(genofile)
    genofile = fullfile(opts.home, compose(chrPatt, "Y"));
    if ~isfile(genofile)
        genofile = fullfile(opts.home, compose(chrPatt, "y"));
    end
end

geno.bfile = genofile;
geno.regenie_merger = opts.merger;
wesinfo = {tab.ann};
geno = {geno};

end