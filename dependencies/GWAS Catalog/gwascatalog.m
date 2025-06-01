function [sums, stu] = gwascatalog(key, varargin)

% implementation of methods described here (GWAS Catalog database):
% https://www.ebi.ac.uk/gwas/rest/docs/api
% 
% Note: this function is still under development.
% 
% INPUTS
% key (required):        A string array (only for
%                        singleNucleotidePolymorphisms) or scalar
%                        corresponding to the search keyterm for
%                        entryPoint.
% entryPoint (optional): Scalar string or char array for entry points to be
%                        used: "associations", "efoTraits" (default),
%                        "studies", "singleNucleotidePolymorphisms" or
%                        "snpLocation".

% input parser ------------------------------------------------------------
[sums, stu] = deal([]);
p = inputParser;
p.CaseSensitive = false;
p.StructExpand = true;
addRequired(p, 'key', ...
    @(x) validateattributes(x,{'char','string'}, {'nonempty'}));
addOptional(p, 'entryPoint', 'singleNucleotidePolymorphisms',...
    @(x) validateattributes(x,{'char','string'},...
    {'nonempty', 'scalartext'}));

p.parse(key, varargin{:});
p = p.Results;

if ~isstring(p.key)
    p.key = string(p.key);
end

if p.entryPoint ~= "singleNucleotidePolymorphisms" ...
        && numel(p.key) > 1
    error('key size cannot be > 1!')
end
% -------------------------------------------------------------------------
% check the presence of "gwas-efo-trait-mappings.mat" and
% "gwas-catalog-studies_ontology-annotated.mat".
homedir = fileparts(which('gwascatalog.m'));
if ~exist(fullfile(homedir, 'gwas-efo-trait-mappings.mat'), 'file') ||...
    ~exist(fullfile(homedir, ...
    'gwas-catalog-studies_ontology-annotated.mat'), 'file')
    gwascatalogFTP(homedir);
end

fprintf('\n********************* GWAS Catalog *************************\n')
getURL = "https://www.ebi.ac.uk/gwas/rest/api/";
headers = {'Content-Type' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', 5000);

% % entry points into the API service
% gwasEntry = webread(getURL, options);
% gwasEntry = gwasEntry.x_links;
% gwasEntryFields = fieldnames(gwasEntry);
% if ~any(ismember(lower(gwasEntryFields), lower(p.entryPoint)))
%     error('entryPoint was not found! only %s are allowed!',...
%         strjoin(gwasEntryFields, ' | '))
% end

% if string(p.entryPoint) ~= "singleNucleotidePolymorphisms" && ...
%     string(p.entryPoint) ~= "snpLocation"
if string(p.entryPoint) == "efoTraits"
    % search within traits (EFO) gwas-efo-trait-mappings.mat
    efo = load(fullfile(homedir, 'gwas-efo-trait-mappings.mat'));
    efo = efo.efo; efo.id = (1:size(efo,1))';

    searchFlag = true;
    while searchFlag
        if isempty(p.key)
            p.key = input('search trait keyword (case insensitive): ', 's');
        end
        foundEFOIDX = contains(lower(efo.(2)), lower(p.key));
        if ~any(foundEFOIDX)
            contAns = input('nothing found :( wanna try again?[Y/n] ', 's');
            if strcmpi(contAns, 'n')
                searchFlag = false;
            end
            p.key = [];
        else

            foundEFO = (efo(foundEFOIDX, 1:3));
            if numel(unique(foundEFO.(3))) > 1 % multiple traits found
                foundEFO = groupsummary(foundEFO, {'EFO URI','EFO term'});
                foundEFO.id = (1:size(foundEFO,1))';
                foundEFO = movevars(foundEFO, 'id', 'before', 'EFO URI');
                foundEFO.GroupCount = [];
                disp(foundEFO);
                selKey = [];  
                while isempty(selKey) || numel(selKey) > 1 || ...
                    selKey > numel(foundEFO) || selKey < 1
                    selKey = input('? ');
                end
            else
                selKey = 1;
             end

            foundEFO = foundEFO.("EFO URI")(selKey);
            searchFlag = false;
        end
    end
    
else
    foundEFO = [];
end

if string(p.entryPoint) == "efoTraits" % get _all_ association data (summary stat) for foundEFO (found trait)
    
    sums = webread(getURL + "/" + p.entryPoint + "/" + foundEFO + ...
        "/associations", options);
    stu = webread(getURL + "/" + p.entryPoint + "/" + foundEFO + ...
        "/studies", options);
    sums = sums.x_embedded.associations;
    stu = stu.x_embedded.studies;
    stu = rmfield(stu, {'gxe', 'gxg', 'snpCount', 'qualifier', 'imputed',...
        'pooled', 'studyDesignComment', 'fullPvalueSet', 'platforms',...
        'genotypingTechnologies', 'userRequested'});
    if ~isempty(sums)
        [sums, stu] = pruneGWASstruct(sums, stu);
    end
    
elseif string(p.entryPoint) == "studies" % get association data per each study (on foundEFO)
    
    stu = load(fullfile(homedir, ...
        'gwas-catalog-studies_ontology-annotated.mat'));
    stu = stu.stu;

    if isempty(foundEFO)
        stu(stu.("STUDY ACCESSION") ~= key, :) = [];
        stu(:, [1:4, 11, 16]) = [];
        assert(~isempty(stu), "Study was not found!")
    else
        stu = stu(contains(stu.(14), foundEFO), :);
        stu(:, [1:4, 11, 16]) = [];
        fprintf('found %d studies for trait: %s\n', size(stu, 1), foundEFO)
    end
    
    sums = ({});
    progressBar = floor(linspace(1, size(stu, 1)+1, 11));
    fprintf('fetching association data [           ]')
    for i = 1:size(stu, 1)

        if isempty(foundEFO)
            sumsTemp = webread(getURL + p.entryPoint + "/" + ...
                stu.("STUDY ACCESSION")(i) + ...
                "/associations", options);
        else
            sumsTemp = webread(getURL + p.entryPoint + "/" + ...
                stu.("STUDY ACCESSION")(i) + ...
                "/associations?projection=associationByStudy", options);
        end
        progressTxt = [repmat('=', 1, sum(i >= progressBar)),'>',...
            repmat(' ', 1, 10-sum(i >= progressBar))];
        fprintf(repmat('\b', 1, 12))
        fprintf('%s]', progressTxt)
        sumsTemp = sumsTemp.x_embedded.associations;
        if ~isempty(sumsTemp)
            sumsTemp = pruneGWASstruct(sumsTemp, stu.("STUDY ACCESSION")(i));
            sumsTemp(sumsTemp.traitEFO ~= foundEFO, :) = [];
            sums{i, 1} = sumsTemp;
        end
    end
    fprintf('\n')
    
elseif string(p.entryPoint) == "singleNucleotidePolymorphisms" % fetch all available associations per SNP
    
    sums = ({});
    progressBar = floor(linspace(1, size(p.key, 1)+1, 11));
    fprintf('fetching association data [           ]')
    for i = 1:size(p.key, 1)
        try
            sumsTemp = webread(getURL + "/" + p.entryPoint + "/" + ...
                p.key(i) + ...
            "/associations?projection=associationBySnp", options);
        catch % not found
            sumsTemp.x_embedded.associations = [];
        end
        progressTxt = [repmat('=', 1, sum(i >= progressBar)),'>',...
            repmat(' ', 1, 10-sum(i >= progressBar))];
        fprintf(repmat('\b', 1, 12))
        fprintf('%s]', progressTxt)
        sumsTemp = sumsTemp.x_embedded.associations;
        if ~isempty(sumsTemp)
            sumsTemp = pruneSNPstruct(sumsTemp);
            sums{i, 1} = sumsTemp;
        end
        sumsTemp = [];
    end
    fprintf('\n')
    
elseif string(p.entryPoint) == "snpLocation" % fetch associations per location
    
    fprintf('under development!\n')
%     sums = ({});
%     progressBar = floor(linspace(1, size(p.key, 1), 11));
%     fprintf('fetching association data [           ]')
%     for i = 1:size(p.key, 1)
%         try
%             sumsTemp = webread(getURL + "/" + p.entryPoint + "/" + ...
%                 p.key(i), options);
%         catch % not found
%             sumsTemp.x_embedded.associations = [];
%         end
%         progressTxt = [repmat('=', 1, sum(i > progressBar)),'>',...
%             repmat(' ', 1, 10-sum(i > progressBar))];
%         fprintf(repmat('\b', 1, 12))
%         fprintf('%s]', progressTxt)
%         sumsTemp = sumsTemp.x_embedded.associations;
%         if ~isempty(sumsTemp)
%             sumsTemp = pruneSNPstruct(sumsTemp);
%             sums{i, 1} = sumsTemp;
%         end
%     end
%     fprintf('\n')
    
end

end % END

%% Subfunctions ===========================================================
function gwascatalogFTP(pth)
% prepares MAT files (annotations, studies, traits, etc.) for gwascatalog
% function. gwascatalogFTP connets to GWAS Catalog FTP at
% ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest and fetches
% details on studies/traits available on GWAS Catalog (e.g. EFO codes for
% traits). This function needs to be run before using gwascatalog, if user
% intends use keywords (e.g. fatty liver) to search for traits (key vargin
% in gwascatalog).

fprintf('fetching data from ftp.ebi.ac.uk...\n')
ftpobj = ftp('ftp.ebi.ac.uk');
goto = "/pub/databases/gwas/releases/latest";
cd(ftpobj, goto);
mget(ftpobj, 'gwas-efo-trait-mappings.tsv', pth);
try % old
    study_file = 'gwas-catalog-studies_ontology-annotated.tsv';
    mget(ftpobj, study_file, pth);
catch
    study_file = 'gwas-catalog-download-studies-v1.0.3.1.txt';
    mget(ftpobj, study_file, pth);
end
close(ftpobj);

% create mat file(s) for gwascatalog function
efo = readtable(fullfile(pth, 'gwas-efo-trait-mappings.tsv'), 'FileType', 'text',...
    'PreserveVariableNames', true, 'TextType', 'string');
tmp1 = regexp(efo.(3),'(?<=/)\w*$','match');
idx = ~cellfun(@isempty, tmp1);
efo.(3) = strings(height(efo), 1);
efo.(3)(idx) = vertcat(tmp1{:});

tmp1 = regexp(efo.(5),'(?<=/)\w*$','match');
idx = ~cellfun(@isempty, tmp1);
efo.(5) = strings(height(efo), 1);
efo.(5)(idx) = vertcat(tmp1{:});

efo.Properties.Description = string(datetime("today"));
save(fullfile(pth, 'gwas-efo-trait-mappings.mat'), 'efo');
fprintf('- gwas-efo-trait-mappings.mat was created.\n') 
delete(fullfile(pth, 'gwas-efo-trait-mappings.tsv'));
clear efo

stu = readtable(fullfile(pth, study_file), ...
    'FileType', 'text',...
    'PreserveVariableNames', true, 'TextType', 'string');
stu.Properties.Description = string(datetime("today"));
save(fullfile(pth, 'gwas-catalog-studies_ontology-annotated.mat'), 'stu')
fprintf('- gwas-catalog-studies_ontology-annotated.mat was created.\n') 
delete(fullfile(pth, study_file))

end % END

%% ------------------------------------------------------------------------
function [S, stu] = pruneGWASstruct(S, stu)

if isstruct(stu)
    efoFlag = true;
else
    efoFlag = false;
end

% extract loci info -------------------------------------------------------
loci = struct2table(arrayfun(@(x)x.loci.strongestRiskAlleles, S));
loci = string(loci.riskAlleleName);
lociSplit = reshape(split(loci, '-'), numel(loci), 2);

% fetch variant info ------------------------------------------------------
headers = {'Content-Type' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', 5000);

if efoFlag 
    snpsList = arrayfun(@(x)x.x_links.snps.href, stu, 'UniformOutput', false);
    stu_trait = string(arrayfun(@(x)join(unique(string(x.diseaseTrait.trait)),';'),...
        stu, 'UniformOutput', false));
    stu_pmid = double(string(arrayfun(@(x)x.publicationInfo.pubmedId, ...
        stu, 'UniformOutput', false)));
    stu = rmfield(stu, {'diseaseTrait', 'publicationInfo', 'x_links'});
    
    anc_struct = parseAncestryInfo(stu); % get ancestry data

    stu = rmfield(stu, 'ancestries');
    stu = struct2table(stu, 'AsArray', true);
    stu = [stu, struct2table(anc_struct, 'AsArray', true)];
    clear anc_struct
    stu.pmid = stu_pmid;
    stu.trait = stu_trait;
    clear anc_struct stu_pmid   
    
else % "studies"
    snpsList = "https://www.ebi.ac.uk/gwas/rest/api/studies/" + ...
        stu + "/snps";
   
end


[snps_chr, snps_pos, snps_rsId, snps_gene] = deal([]);
progressBar = floor(linspace(1, numel(snpsList), 11));
fprintf('fetching variant info [           ]')
    
for i = 1:numel(snpsList)
    snptemp = webread(snpsList{i}, options);
    snptemp = snptemp.x_embedded.singleNucleotidePolymorphisms;
    if isempty(snptemp)
        continue
    end
    
    snps_rsId_temp = string(arrayfun(@(x)x.rsId, snptemp, 'UniformOutput', false));
    
    % if field "locations" is empty, it's gonna be 0 or 1. otherwise I
    % should change it! this is just a sanity check!
    icheck = (sum(arrayfun(@(x)isempty(x.locations), snptemp))/...
        numel(arrayfun(@(x)isempty(x.locations), snptemp)));
%     if ~ismember(icheck, [0, 1])
%         error('some locations are missing but not all!')
%     end
    
    if ~any(arrayfun(@(x)isempty(x.locations), snptemp))
        snptemp2 = {snptemp.locations};
        snps_chr_tmp = repmat("-", numel(snps_rsId_temp), 1);
        snps_pos_tmp = repmat("NaN", numel(snps_rsId_temp), 1);
        for j = 1:numel(snptemp2)
            try
                snps_chr_tmp(j) = string(snptemp2{j}.chromosomeName);
                snps_pos_tmp(j) = string(snptemp2{j}.chromosomePosition);
            catch
            end
        end

        snps_chr = [snps_chr; snps_chr_tmp];
        snps_pos = [snps_pos; snps_pos_tmp];
    else
        snps_chr = [snps_chr; repmat("-", numel(snps_rsId_temp), 1)];
        snps_pos = [snps_pos; repmat("NaN", numel(snps_rsId_temp), 1)];
    end
    
    snps_rsId = [snps_rsId; snps_rsId_temp];
    
    icheck = (sum(arrayfun(@(x)isempty(x.genomicContexts), snptemp))/...
        numel(arrayfun(@(x)isempty(x.genomicContexts), snptemp)));
    if ~ismember(icheck, [0, 1])
%         warning('some genomicContexts are missing but not all!')
    end   
    
    if all(arrayfun(@(x)isempty(x.genomicContexts), snptemp))
        snps_gene_temp = repmat("-", numel(snps_rsId_temp), 1);
    else
        snptemp = {snptemp.genomicContexts};
        snps_gene_temp = repmat("-", numel(snps_rsId_temp), 1);
        for j = 1:numel(snptemp)
            try
                [~, nearestIdx] = min([snptemp{j}.distance]); % nearest gene
                snps_gene_temp(j, 1) = snptemp{j}(nearestIdx).gene.geneName;
            catch
            end
        end
    end
  
    snps_gene = [snps_gene; snps_gene_temp];
    clear snptemp
    
    progressTxt = [repmat('=', 1, sum(i > progressBar)),'>',...
            repmat(' ', 1, 10-sum(i > progressBar))];
    fprintf(repmat('\b', 1, 12))
    fprintf('%s]', progressTxt)
end
fprintf('\n')

snps = table(snps_chr, snps_pos, snps_rsId, snps_gene, ...
    'VariableNames', {'chr', 'pos', 'varid', 'gene'});
clear snps_chr snps_pos snps_rsId snps_gene snps_gene_temp

[f_loci, f_keep] = ismember(lociSplit(:, 1), snps.varid);
if any(f_keep < 1) || any(~f_loci)
    if efoFlag
        % there is an issue with GWAS catalog: there may be some variants
        % in efoTrait data which do not exist in studies associated with
        % that efoTrait (e.g. variant rs123 is in GCST1 which is not listed
        % under efoTrait studies in API).
        f_loci = find(~f_loci);
        snps_append = strings;
        for i = 1:numel(f_loci)
            snps_append(i, 3) = lociSplit(f_loci(i), 1);
            
            try
                snpstemp = webread("https://www.ebi.ac.uk/gwas/rest/api/" + ...
                    "singleNucleotidePolymorphisms/" + ...
                    lociSplit(f_loci(i), 1), options);
            catch
                snps_append(i, 1) = "-"; snps_append(i, 2) = "NaN";
                snps_append(i, 4) = "-";
                continue
            end
            
            if ~any(arrayfun(@(x)isempty(x.locations), snpstemp))
                snps_append(i, 1) = string(arrayfun(...
                    @(x)x.locations.chromosomeName, snpstemp,...
                    'UniformOutput', false));
                snps_append(i, 2) = string(arrayfun(...
                    @(x)x.locations.chromosomePosition, snpstemp,...
                    'UniformOutput', false));
            else
                snps_append(i, 1) = "-"; snps_append(i, 2) = "NaN";
            end
            
            if any(arrayfun(@(x)isempty(x.genomicContexts), snpstemp))
                snps_append(i, 4) = "-";
            else
                snpstemp = snpstemp.genomicContexts;
                [~, nearestIdx] = min([snpstemp.distance]); % nearest gene
                snps_append(i, 4) = snpstemp(nearestIdx).gene.geneName;
            end  
        end
        snps_append = array2table(snps_append);
        snps_append.Properties.VariableNames = {'chr', 'pos', 'varid', 'gene'};
        snps = [snps; snps_append];
        [~, f_keep] = ismember(lociSplit(:, 1), snps.varid);
        
    else % associations per study "must" match fetched variants
        % something went wrong, mismatch between reported and found variants
        error('could not match reported and found variants!')
    end
end

snps = snps(f_keep, :);
snps.ea = lociSplit(:, 2);
clear lociSplit
snps.pos = double(snps.pos);
snps.varid = regexprep(snps.varid, '^chr', ''); % remove chr tag from variant id

% prune S -----------------------------------------------------------------
loci = table;
if ~efoFlag % "studies"
    loci.trait = string(arrayfun(@(x)x.efoTraits.trait, S, ...
        'UniformOutput', false));
    loci.traitEFO = string(arrayfun(@(x)x.efoTraits.shortForm, S, ...
        'UniformOutput', false));
    S = rmfield(S, {'efoTraits', 'loci','study', 'multiSnpHaplotype', ...
        'pvalueExponent', 'pvalueMantissa', 'description', 'x_links'});

else %efoTraits
    S = rmfield(S, {'loci', 'multiSnpHaplotype', 'pvalueExponent',...
        'pvalueMantissa', 'description', 'x_links', 'lastMappingDate',...
        'lastUpdateDate'});
end

S = struct2table(S, 'AsArray', true); % convert to table

colNames = S.Properties.VariableNames;
S.Properties.VariableNames(ismember(colNames, 'riskFrequency')) = {'eaf'};
S.Properties.VariableNames(ismember(colNames, 'orPerCopyNum')) = {'or2beta'};
S.Properties.VariableNames(ismember(colNames, 'betaNum')) = {'beta'};
S.Properties.VariableNames(ismember(colNames, 'standardError')) = {'se'};

% take care of numeric variables 
editCols = {'eaf', 'or2beta', 'beta', 'se', 'pvalue'};
check_double = varfun(@isnumeric, S, 'InputVariables', editCols);
editCols = editCols(~table2array(check_double));
check_empty = varfun(@(x) cellfun(@isempty, x), S, 'InputVariables', ...
    editCols); 
check_empty = check_empty.Variables;

for i = 1:numel(editCols)
    S.(editCols{i})(check_empty(:, i)) = {'NaN'};
    S.(editCols{i}) = double(string(S.(editCols{i})));
end

if any(S.beta)
    betaSign = lower(unique(S.betaDirection(~isnan(S.beta))));
    if any(ismember(betaSign, {'decrease', 'increase'}))
        S.betaDirection(isnan(S.beta)) = {'NA'};
        f_negative = ismember(S.betaDirection, 'decrease');
        S.beta(f_negative) = -1.*S.beta(f_negative);
        S.betaDirection = [];
    end
end
if any(S.or2beta)
    S.or2beta(~isnan(S.or2beta)) = log(S.or2beta(~isnan(S.or2beta)));
end

S = [snps, S, loci]; 
end % END

%% ------------------------------------------------------------------------
function S = pruneSNPstruct(S)

% fetch study info --------------------------------------------------------
stu = struct;
for i = 1:numel(S)
    stu_temp = webread(S(i).x_links.study.href);
    diseaseTrait = string({stu_temp.diseaseTrait.trait});
    pubmedId = double(string(stu_temp.publicationInfo.pubmedId));
    stu_temp = rmfield(stu_temp, {'gxe', 'gxg', 'snpCount', ...
        'qualifier', 'imputed',...
        'pooled', 'studyDesignComment', 'fullPvalueSet', 'platforms',...
        'genotypingTechnologies', 'userRequested', 'publicationInfo',...
        'diseaseTrait', 'x_links'});
    anc_struct = parseAncestryInfo(stu_temp);
    stu_temp = rmfield(stu_temp, 'ancestries');
    fnames1 = fieldnames(stu_temp);
    fnames2 = fieldnames(anc_struct);
    for j = 1:numel(fnames1)
        stu(i).(fnames1{j}) = stu_temp.(fnames1{j});
    end
    for j = 1:numel(fnames2)
        stu(i).(fnames2{j}) = anc_struct.(fnames2{j});
    end
    stu(i).diseaseTrait = diseaseTrait;
    stu(i).pubmedId = pubmedId;
    clear anc_struct diseaseTrait pubmedId stu_temp fnames1 fnames2
end
stu = struct2table(stu, 'AsArray', true);

% extract loci/snp info ---------------------------------------------------
snps = {S.snps};
loci = {S.loci};
traits = {S.efoTraits};
S = rmfield(S, {'loci', 'snps', 'multiSnpHaplotype', 'pvalueExponent',...
    'pvalueMantissa', 'description', 'x_links', ...
    'efoTraits', 'pvalueDescription'});

snps_info = strings;
for j = 1:numel(snps)
    if isempty(snps{j})
        snps_info(j, 1:7) = repmat("NR", 1, 7);
        continue
    end
    snps_temp = snps{j};
    loci_temp = {loci{j}.strongestRiskAlleles};
    traits_temp = traits{j};
    snps_temp = rmfield(snps_temp, ...
        {'functionalClass', 'lastUpdateDate', 'merged'});
    
    snps_info(j, 1) = join(unique(arrayfun(@(x)string(x.rsId), snps_temp), 'stable'), ';');
    snps_locations = {snps_temp.locations};
    snps_genomicContexts = {snps_temp.genomicContexts};
    
    snps_locations_struct = strings;
    for k = 1:numel(snps_locations)
        try
            riskAlleleName = loci_temp{k}.riskAlleleName;
        catch
            riskAlleleName = loci_temp{1}(k).riskAlleleName;
        end
        riskAlleleName = split(riskAlleleName, '-');
        snps_locations_struct(k, 1) = riskAlleleName(2);
        snps_locations_struct(k, 2:3) = [snps_locations{k}(1).chromosomeName, ...
            string(snps_locations{k}(1).chromosomePosition)];
        [~, nearestGeneIdx] = min([snps_genomicContexts{k}.distance]);
        snps_locations_struct(k, 4) = ...
            snps_genomicContexts{k}(nearestGeneIdx).gene.geneName;
    end
    
    if size(snps_locations_struct, 1) > 1
        snps_locations_struct = join(snps_locations_struct', ';')';
    end
    
    if ~isempty(traits_temp)
        snps_locations_struct(1, 5) = join(string({traits_temp.trait}), '; ');
        snps_locations_struct(1, 6) = join(string({traits_temp.shortForm}), '; ');
    else
        snps_locations_struct(1, 5:6) = "";
    end
    
    snps_info(j, 2:7) = snps_locations_struct;
end
clear snps_locations_struct snps_genomicContexts snps_locations snps_temp snps loci_temp loci traits traits_temp

snps_info(:, 1) = regexprep(snps_info(:, 1), '^chr', ''); % remove chr tag from variant id
snps_info = array2table(snps_info, 'VariableNames', {'varid', 'ea', 'chr', ...
    'pos', 'gene', 'trait', 'traitEFO'});

% prune S -----------------------------------------------------------------
S = struct2table(S, 'AsArray', true); % convert to table

colNames = S.Properties.VariableNames;
S.Properties.VariableNames(ismember(colNames, 'riskFrequency')) = {'eaf'};
S.Properties.VariableNames(ismember(colNames, 'orPerCopyNum')) = {'or2beta'};
S.Properties.VariableNames(ismember(colNames, 'betaNum')) = {'beta'};
S.Properties.VariableNames(ismember(colNames, 'standardError')) = {'se'};

% take care of numeric variables 
editCols = {'eaf', 'or2beta', 'beta', 'se', 'pvalue'};
check_double = varfun(@isnumeric, S, 'InputVariables', editCols);
editCols = editCols(~table2array(check_double));
check_empty = varfun(@(x) cellfun(@isempty, x), S, 'InputVariables', ...
    editCols); 
check_empty = check_empty.Variables;

for i = 1:numel(editCols)
    S.(editCols{i})(check_empty(:, i)) = {'NaN'};
    S.(editCols{i}) = double(string(S.(editCols{i})));
end

if any(S.beta)
    betaSign = lower(unique(S.betaDirection(~isnan(S.beta))));
    if any(ismember(betaSign, {'decrease', 'increase'}))
        S.betaDirection(isnan(S.beta)) = {'NA'};
        f_negative = ismember(S.betaDirection, 'decrease');
        S.beta(f_negative) = -1.*S.beta(f_negative);
    end
end
if ismember('betaDirection', S.Properties.VariableNames)
    S.betaDirection = [];
end
if any(S.or2beta)
    S.or2beta(~isnan(S.or2beta)) = log(S.or2beta(~isnan(S.or2beta)));
end

S = [snps_info, S, stu]; 
S.pos = double(S.pos);
S = movevars(S, {'beta' ,'or2beta', 'se', 'pvalue', 'range', 'eaf'}, ...
    'Before', 'snpType');
end % END

%% ------------------------------------------------------------------------
function anc_struct = parseAncestryInfo(stu)
% get ancestry info for each study
anc_struct = struct;
for i = 1:numel(stu)
    anc = stu(i).ancestries;
    countryOfRecruitment = {anc.countryOfRecruitment};
    ancestralGroups = {anc.ancestralGroups};
    anc = rmfield(anc, {'countryOfOrigin', 'countryOfRecruitment'});
    anc = struct2table(anc, 'AsArray', true);

    anc_country = strings;
    for j = 1:numel(countryOfRecruitment)
        if isempty(countryOfRecruitment{j})
            anc_country(j, 1:3) = ["NR", "NR", "NR"];
            continue
        end
        countryOfRecruitment_temp = ...
            struct2table(countryOfRecruitment{j}, 'AsArray', true);
        emp = varfun(@(x)cellfun(@isempty, x), countryOfRecruitment_temp);
        emp = emp.Variables;
        countryOfRecruitment_temp = countryOfRecruitment_temp.Variables;
        countryOfRecruitment_temp(emp) = {'NR'};
        countryOfRecruitment_temp = string(countryOfRecruitment_temp);
        if size(countryOfRecruitment_temp, 1) > 1
            countryOfRecruitment_temp = ...
                join(countryOfRecruitment_temp', ';')';
        end
        anc_country(j, 1:3) = countryOfRecruitment_temp;
    end

    for j = 1:numel(ancestralGroups)
        if isempty(ancestralGroups{j})
            anc_country(j, 4) = "NR";
            continue
        end
        anc_country(j, 4) = ...
            join(string({ancestralGroups{j}.ancestralGroup}), ';');
    end

    anc.majorArea = anc_country(:, 1);
    anc.region = anc_country(:, 2);
    anc.country = anc_country(:, 3);
    anc.ancestralGroups = anc_country(:, 4);
    clear anc_country countryOfRecruitment countryOfRecruitment_temp ancestralGroups
    anc = groupsummary(anc, 'type', @(x)join(string(x),' | '));
    anc.GroupCount = [];
    anc.Properties.VariableNames = ...
        replace(anc.Properties.VariableNames, 'fun1_', '');

    initial_idx = ismember(anc.type, 'initial');
    if any(initial_idx)
        anc_struct(i).initial_N = anc.numberOfIndividuals(initial_idx);
        anc_struct(i).initial_ancestralGroups = anc.ancestralGroups(initial_idx);
        anc_struct(i).initial_majorArea = anc.majorArea(initial_idx);
        anc_struct(i).initial_region = anc.region(initial_idx);
        anc_struct(i).initial_country = anc.country(initial_idx);
    else
        anc_struct(i).initial_N = "NR";
        anc_struct(i).initial_ancestralGroups = "NR";
        anc_struct(i).initial_majorArea = "NR";
        anc_struct(i).initial_region = "NR";
        anc_struct(i).initial_country = "NR";
    end
    
    rep_idx = ismember(anc.type, 'replication');
    if any(rep_idx)
        anc_struct(i).rep_N = anc.numberOfIndividuals(rep_idx);
        anc_struct(i).rep_ancestralGroups = anc.ancestralGroups(rep_idx);
        anc_struct(i).rep_majorArea = anc.majorArea(rep_idx);
        anc_struct(i).rep_region = anc.region(rep_idx);
        anc_struct(i).rep_country = anc.country(rep_idx);
    else
        anc_struct(i).rep_N = "NR";
        anc_struct(i).rep_ancestralGroups = "NR";
        anc_struct(i).rep_majorArea = "NR";
        anc_struct(i).rep_region = "NR";
        anc_struct(i).rep_country = "NR";
    end

end

end % END