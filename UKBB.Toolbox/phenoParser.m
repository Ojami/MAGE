function UKB_STRUCT_ALL = phenoParser(opts)
% parses UK Biobank processed data basket (output of UKBBasketParser).
% Accepts phenotype definitions in a file or can search within UKB
% dictionary in an interactively. 
% 
% Oveis Jamialahmadi, University of Gothenburg, April 2021.
% 
% @14/08/2022: a bug was fixed.
% @27/08/2022: 'merge' flag was added, while merges all queried phenotypes
%               into one unified trait. This is useful for categorical
%               traits where multiple data-fields are used for the
%               definition of the target phenotype (e.g. in-patient and
%               self-reported resources). When using 'merge' flag, first
%               'tag' element is used for the merged phenotype and there is
%               no need to input a vector for 'tag' argument.
% @07/10/2022: 'desc' option was added for additional descriptions for each
%               trait. This is written to 'info.desc' subfield. 
% @31/10/2022: minor modifications for HESIN-DEATH data.
% @31MAR2023: a bug was fixed.
% @16MAY2024: 'surv' falg was added for definition of time-to-even data.
%             Note that phenoParser defines a comprehensive trait, meaning
%             it writes also competing risk/death information so they can
%             be properly analyzed in downstream steps (e.g.
%             subdistribution or cause-specific hazard models). Moreover,
%             it's possible to consider landmark time:
%             https://www.ahajournals.org/doi/full/10.1161/CIRCOUTCOMES.110.957951
%             https://www.sciencedirect.com/science/article/pii/S0002929720303232?via%3Dihub
%             https://link.springer.com/article/10.1186/s12916-018-1063-1
%             https://www.sciencedirect.com/science/article/pii/S016882782200023X
%             generally a landmark period of 2 yrs is considered (above),
%             however there are other cut-off points (e.g. 90 days) which
%             can be used:
%             https://www.nature.com/articles/s41467-023-44512-4
%             For genetic studies, however the starting time can be set to
%             0 and there is no reverse-causality problem:
%             https://stats.stackexchange.com/questions/580701/choice-of-time-origin-for-survival-analysis-when-no-specific-event-or-beginning
%             'landmark' flag can be used for this purpose (default is 0.5
%             year)
%            
%             Overall a good approach can be found here:
%             https://www.bmj.com/content/363/bmj.k4168.full
% @25MAY2024: 'landmark' option is deprecated now because it's easy for
%              calculate in downstream analyses and results in unnecessary
%              memory use for saving/loading phenotypes.
% @02JUNE2024: To define a time-to-event phenotype for genome-wide studies,
%              one should consider age-of-onset (AAO), since the exposure
%              in this case is a genetic variant existing at time of birht.
%              In this case, one should not adjust the model for age:
%              
%              https://www.nature.com/articles/s41467-021-22538-w#Sec9:
%              "We did not include any covariates of age or year of birth
%              because these are directly associated with our phenotypes."
% 
%              Other useful refereces:
%              - https://www.nature.com/articles/s41398-022-01782-8#Sec2
%              - https://github.com/WenjianBI/SPACox
%              - https://www.nature.com/articles/s41467-023-41210-z
%              
%              on left-truncation (if using AAO), adjusted for age also as
%              covariate: Surv(enrollment_age, AAO, event) = Surv(tt0-tt, tt, censor==1)
%              - Main reference: https://academic.oup.com/biomedgerontology/article/53A/5/M337/588259
%              - https://www.sciencedirect.com/science/article/pii/S2666247721000580
%              - https://www.nature.com/articles/s41416-019-0465-y
%              - On age-of-onset: https://www.nature.com/articles/s41467-022-32885-x
% 
%             on case ascertainment and left-truncation, not adjusted for age as covariate:
%             - https://www.nature.com/articles/s41467-020-17374-3#Sec7
%               "All participants diagnosed at dates prior to enrollment in
%               the UK Biobank were considered prevalent at baseline, while
%               participants diagnosed after enrollment were considered
%               incident"
% 
%             When the outcome is death, left-truncation is underpowered,
%             and if not handled, is biased. Because a subject who died
%             before enrollment could not be in the study, ergo not at
%             risk.
%             - https://www.sciencedirect.com/science/article/pii/S2405844024010065
%               "In survival analysis, left truncation bias [46] was
%               corrected by considering only the time between patient
%               sampling and time of death. Left truncation bias occurs
%               when risk of death is measured over a time in which it
%               could not have occurred; by definition a patient who had
%               already died could not have been recruited into a study.
%               The rationale is to exclude survival time that occurred
%               before sampling, because the patient could not, by
%               definition, have died during this period. Analysis of left
%               truncated data can lead to false positive associations
%               [46], however, left truncation bias correction can lead to
%               under-powered analysis because of the loss of information."
% 
% @19SEP2024: 'ukbrap' flag was added to fetch data from UKB-RAP on
%             DNANexus. See dx_extract_dataset function for more details.
% @22APR2025: 'threads' option was added to use parallele loop (parfor) for
%             query matching to inpatient data. Ideally this should be used
%             when ukbrap is used. Recommended threads < 10.


arguments
    opts.method {mustBeTextScalar, mustBeMember(opts.method, ["merge", "file", "manual"])} = "manual"
    opts.file {mustBeTextScalar, mustBeFile} % should be set if method is "file"
    opts.phenodir {mustBeFolder} = fullfile(fileparts(which('phenoParser.m')), 'UKB_PHENO')
    opts.basketdir {mustBeFolder} = fileparts(which('phenoParser.m')) % basket folders have prefix 'UKBFileParser_'
    opts.dictdir {mustBeFolder} = fullfile(fileparts(which('phenoParser.m')), "UKB_DICTIONARY") % dictionary/codings folder
    % query to search for, can be data-field or (part of) the description
    % of data-field(s). If more than 1 
    opts.query {mustBeText, mustBeVector}
    opts.tag {mustBeText, mustBeVector} % phenotype tag (name) for each element of query, effective when only one phenotype is generated for each query (e.g. ICD10 codes)
    opts.desc {mustBeText, mustBeVector} % additional descriptions (e.g. definitions, citations, etc). should be of the same size as 'tag'
    
    % exclusion criteria: codes (for each element of query) to be excluded
    % when creating a categorical phenotype (case/control) in this case,
    % 'remove' must be a text vector. if 'surv' is used, this should be a
    % struct with each fields named as "dfxxx_n" where xxx is the
    % data-field number from UKBB and n is the instance, or alternatively
    % "cdxxx_n" where xxx is the coding string (similar to the 'coding'
    % argument) and n is the instance. In case of multiple phenotype
    % definitions when 'surv' is true, 'remove' should be a struct array
    % with the same size (size(, 2)) as input query. If n is not provided,
    % all instances will be used. Example of such struct:
    %   ex = struct;
    %   ex.df20002_0 = "1136;1155"; % non-cancer self reported at baseline
    %   ex.df20001_0 = "1070"; % cancer self reported at baseline
    %   ex.cd19 = "C810"; % coding-19 (HESIN, death, cancer registry) at all instances
    %   
    % for multiple pheno, expand the struct, e.g. ex(2).dfxxx_n ...

    opts.remove 
    opts.coding {mustBeText, mustBeVector} % codings to be used for each element of query. e.g. if coding is 19, only data-fields having this coding are kept.
    opts.df {mustBeText, mustBeVector} % similar to 'codings' but keeps only these data-fields for each element of query (e.g. query is ICD-10 code, an df can be 41270)  
    opts.array double {mustBeVector, mustBeNonnegative, mustBeNonNan} % arrays to be read (default: all)
    opts.instance double {mustBeVector, mustBeNonnegative, mustBeNonNan} % instances to be read (default: all)
    opts.all (1,1) logical = false % load/parse all data-fields found in the data dictionary (true), or ask user to select (default: false)
    opts.save (1,1) logical = true % write the phenotype to a mat file (only if 'method' is "manual")
    opts.merge (1,1) logical = false 

    %for prospective/survival analysis
    opts.surv (1,1) logical = true % to create the phenotype ready for time-to-event analysis
    % opts.landmark (1,1) double = 0.5 % in years to exclude cases who developed the outcome before this value

    %@19SEP2024
    opts.ukbrap (1,1) logical = true
    opts.ukbrap_dir {mustBeTextScalar} = "UKBFileParser_RAP" % name of UKB-RAP pheno directory (not path)
    
    %@22APR2025
    opts.threads (1,1) double = nan
end

% check inputs
if opts.method == "manual"
    if ~isfield(opts, 'query')
        error('phenoParser:wrongInput', 'query cannot be empty with "manual" method!')
    end
end

if isfield(opts, "remove")
    if opts.surv % check if exclusion criteria is struct or struct array (multi pheno)
         if ~isstruct(opts.remove)
            error('phenoParser:wrongInput', 'remove option must be a struct/struct array for time-to-event phenotypes. See the doc.')
        end
    else
        if ~isvector(opts.remove)
            error('phenoParser:wrongInput', 'remove option must be a string vector. See the doc.')
        end
    end
end

if ~isfield(opts, 'array'); opts.array = []; end
if ~isfield(opts, 'instance'); opts.instance = []; end

if opts.method == "file"
    fprintf('~~~~~~~~~~~~~~~~~~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~\n')
    fprintf('1-Each phenotype should be stored in a single line\n')
    fprintf('2-Phenotype information should be entered as:\n')
    fprintf('\ttag name|key terms|secondary terms\n')
    fprintf('If no secondary term is available, leave it empty\n')
    fprintf('Multiple key terms should be separated with ;\n')
    fprintf('Multiple secondary terms should be separated wtih ,\n')
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    
    pheno_info = readlines(opts.file);
    pheno_info(pheno_info == "") = [];
    try
        pheno_info = split(pheno_info, '|');
    catch
        fprintf('Error:wrong phenotype info file: missing vertical bar(|)!\n')
        return
    end
    inputTerms.tags = pheno_info(:, 1);
    inputTerms.keyTerms = strtrim(cellfun(@(x)strsplit(x, ';')', pheno_info(:, 2), 'UniformOutput', false));
    inputTerms.secondaryTerms = strtrim(cellfun(@(x)strsplit(x, ';')', pheno_info(:, 3), 'UniformOutput', false));
    empty_secondary_terms = cellfun(@(x)isempty(x{1}), inputTerms.secondaryTerms);
    inputTerms.secondaryTerms(empty_secondary_terms) = {{'NOSECONDARYTERMISSELECTED'}};
    inputTerms = printInputterms(inputTerms);
    % Define options for UKBtraitParser
    UKB_STRUCT_ALL = getUKBtrait_prop(inputTerms);
    
else % "merge" or "manual:
    
    [data_dict, codings] = getUKBdictionary(dir=opts.dictdir, ukbrap=opts.ukbrap);

    if opts.method == "merge" % query terms are merged and then a single phenotype will be generated 
        if ~isfield(opts, 'query')
            fprintf('Note: to define search terms that appear in > 1 data-field, use ; to separate terms\n')
            fprintf('      for instance, diagnos;icd10 identifies both\n')
            fprintf('      Diagnoses - main ICD10 and Diagnoses - secondary ICD10\n')
            fprintf('Note: You can also insert data-field IDs of interest as well:\n')
            fprintf('      for instance, 41202;41204 corresponds to the same as above.\n')
            opts.query = input('Insert search terms: ', 's');
            opts.query = strtrim(strsplit(opts.query, ';'));
        end
    
        inputTerms = searchUKBdictionary(opts.query, data_dict, codings);
        inputTerms = printInputterms(inputTerms);
        UKB_STRUCT_ALL = getUKBtrait_prop(inputTerms);

    else % manual
        UKB_STRUCT_ALL = parsePhenoQuery(opts.query, data_dict, codings, opts);
        return
    end

end

% check if is empty
if isempty(UKB_STRUCT_ALL)
    return
end

% Add/save new traits
for i = 1:numel(UKB_STRUCT_ALL)
    UKB_STRUCT_ALL(i).tag = strtrim(char(UKB_STRUCT_ALL(i).tag));
    add2savedtraits(UKB_STRUCT_ALL(i), opts.phenodir)
end

end % END

%% subfunctions -----------------------------------------------------------
function inputTerms = printInputterms(inputTerms)
Ncoded = cellfun(@numel, inputTerms.secondaryTerms);
for ii = 1:numel(Ncoded)
    inputTerms.coded{ii,1} = repmat({'1'}, Ncoded(ii), 1);
end
fprintf('----------------------------------\n')
fprintf('Input trait properties:\n')
for ii = 1:numel(Ncoded)
    fprintf('Info for pheno %d\n', ii)
    fprintf('tag: %s\n', inputTerms.tags{ii})
    if isfield(inputTerms, 'fieldID')
        for jj = 1:numel(inputTerms.fieldID{ii})
            fprintf('field ID %d: %s\n', jj, inputTerms.fieldID{ii}{jj})
        end
    end
    for jj = 1:numel(inputTerms.keyTerms{ii})
        fprintf('Key term %d: %s\n', jj, inputTerms.keyTerms{ii}{jj})
    end
    for jj = 1:numel(inputTerms.secondaryTerms{ii})
        fprintf('Secondary term %d: %s\n', jj, inputTerms.secondaryTerms{ii}{jj})
    end
    if ii ~= numel(Ncoded)
        fprintf('\n')
    end
end
fprintf('----------------------------------\n')
end

%% ------------------------------------------------------------------------
function UKB_STRUCT_ALL = getUKBtrait_prop(inputTerms)
UKBFileParserFlag = true;
inputUKBfile = 'ukb29150'; %ukb29150 ukb26112

wd = fileparts(which('phenoParser.m'));
% check available ukb baskets
ukb_baskets = dir(fullfile(wd, 'UKBFileParser_*'));
ukb_baskets_datenum = [ukb_baskets.datenum]';
ukb_baskets = {ukb_baskets.name}';

fields = inputTerms.fieldID;
[wrong_baskets_idx, found_baskets_idx] = deal(false(numel(ukb_baskets), 1)); 
for i = 1:numel(ukb_baskets)
    filePath = fullfile(wd, ukb_baskets{i}, 'variableMapper.mat');
    if ~exist(filePath, 'file')
        wrong_baskets_idx(i) = true;
        continue
    end
    idx = load(fullfile(wd, ukb_baskets{i}, 'variableMapper.mat'));
    
    ubiq = false; cnt = 1;
    for j = 1:numel(fields)
        infields = fields{j};
        infields = strcat('x', infields, '_');
        for k = 1:numel(infields)
            if any(contains(idx.variableMapper(:, 2), infields{k}))
                ubiq(cnt) = true;
            else
                ubiq(cnt) = false;
            end
            cnt = cnt + 1;
        end
    end
     % all input data-feilds should be present in this basket. This may not
     % be too specific, as only a portion of input data-fields may be
     % present in a basket. 
    if all(ubiq)
        found_baskets_idx(i) = true;
    end
end
ukb_baskets(wrong_baskets_idx) = [];
ukb_baskets_datenum(wrong_baskets_idx) = [];
ukb_baskets = replace(ukb_baskets, 'UKBFileParser_', '');

% sort by date 
[~, sorted] = sort(ukb_baskets_datenum);
ukb_baskets = ukb_baskets(sorted);
found_baskets_idx = found_baskets_idx(sorted);

if all(~found_baskets_idx)
    fprintf('WARNING:no existing ukb basket covers all input data-fields!\n')
end
found_baskets = repmat("No", numel(ukb_baskets), 1);
found_baskets(found_baskets_idx) = "Yes";
basket_choice = table;
basket_choice.ID = (1:numel(ukb_baskets))';
basket_choice.Basket = string(ukb_baskets);
basket_choice.("Has input DF?") = found_baskets;
disp(basket_choice)

ukbSel = input(sprintf('enter UKB basket ID(default:%s): ', inputUKBfile));
if ~isempty(ukbSel)
    while ukbSel > size(basket_choice, 1)
        fprintf('wrong index!\n')
        ukbSel = input(...
            sprintf('enter UKB basket ID(default:%s): ', inputUKBfile));
    end
    inputUKBfile = ukb_baskets{ukbSel};
    [~,inputUKBfile] = fileparts(inputUKBfile);
end

if strcmp(inputUKBfile, 'HESIN_DEATH') % use directly downloaded hospitalized/death data
    UKB_STRUCT_ALL = ({});
    for i = 1:numel(inputTerms.secondaryTerms{1})
        if numel(inputTerms.secondaryTerms{1}) > 1
            fprintf('secondary term: %d of %d\n', i, numel(inputTerms.secondaryTerms{1}))
        end
        inputTermstmp = inputTerms;
        inputTermstmp.secondaryTerms{1} = inputTerms.secondaryTerms{1}(i);
        inputTermstmp.secondaryTermsMeaning{1} = inputTerms.secondaryTermsMeaning{1}(i);
        inputTermstmp.coded{1} = inputTerms.coded{1}(i);
        UKB_STRUCT_ALL{i} = getHESIN_DEATH(inputTermstmp);
        if numel(inputTerms.secondaryTerms{1}) > 1
            fprintf('-------------------------------\n')
            UKB_STRUCT_ALL{i}.tag = cellstr(join([unique(UKB_STRUCT_ALL{i}.rawUKB).',UKB_STRUCT_ALL{i}.tag{1}], '_'));
        end
    end
    UKB_STRUCT_ALL = vertcat(UKB_STRUCT_ALL{:});
    return
end
inputUKBfile = [inputUKBfile, '.csv'];
optionStruct.matsave = false;
optionStruct.mergeInstances = true;
optionStruct.instanceNum = [];

UKB_eidSTRUCT = getUKBtrait(inputTerms, optionStruct, inputUKBfile, UKBFileParserFlag);

% check if is empty
if numel(UKB_eidSTRUCT) < 2 && all(cellfun(@isempty, UKB_eidSTRUCT{1}))
    fprintf('ERROR: nothing to save!\n')
    UKB_STRUCT_ALL = [];
    return
end
UKB_STRUCT_ALL = combine_UKB_STRUCTs(UKB_eidSTRUCT, inputTerms.tags);

UKB_STRUCT_ALL = modifyUKB_STRUCT_ALLtags(UKB_STRUCT_ALL); % @8/29/2019:Modify tags

end % END

%% ------------------------------------------------------------------------
function add2savedtraits(inTraits, pheno_dir)

for ii = 1:numel(inTraits)
    UKB_STRUCT_ALL = inTraits(ii);
    save_name = string(UKB_STRUCT_ALL.tag);
    save_name = matlab.lang.makeValidName(save_name);
    if ~exist(save_name + ".mat", 'file')
        save(fullfile(pheno_dir, save_name + ".mat"), 'UKB_STRUCT_ALL')
    end
end

end % END

%% ------------------------------------------------------------------------
function inputTerms = searchUKBdictionary(keyterms, data_dict_showcase, codings_showcase)
% searches keyterms within UKB dictionary/encoding files

% Search for key terms
if ~isstring(keyterms)
    keyterms = string(keyterms);
end

if all(~isnan(double(keyterms))) % data-field IDs are provided
    f_data_dict = ismember(data_dict_showcase.FieldID, keyterms);
    if ~any(f_data_dict)
        error('Data-field(s) was not found in UKB')
    end
else
    for ii = 1:numel(keyterms) % loop over key terms
        f_data_dict_all = contains(lower(data_dict_showcase{:, :}), lower(keyterms(ii)));
        f_coding_all = contains(lower(codings_showcase{:, :}), lower(keyterms(ii)));
        if ii < 2
            f_data_dict = f_data_dict_all;
            f_coding = f_coding_all;
        else
            f_data_dict = f_data_dict & f_data_dict_all;
            f_coding = f_coding & f_coding_all;
        end
    end
    
    f_data_dict = any(f_data_dict, 2);
    f_coding = any(f_coding, 2);
    if ~any(f_data_dict) && ~any(f_coding)
        fprintf('Terms were not found in UKB dic/enc files\n')
        return
    end
end

h_dict = string(data_dict_showcase.Properties.VariableNames); % Header of showcase dictionary
h_coding = string(codings_showcase.Properties.VariableNames); % Header of coding dictionary

included_info_dict = ["Path", "Field", "Coding", "Notes"];
info_dict_idx = find(ismember(h_dict, included_info_dict));

s_dict = data_dict_showcase(f_data_dict, :);
s_coding = codings_showcase(f_coding, :);

s_dict = fillmissing(s_dict, 'constant', "NA");
s_coding = fillmissing(s_coding, 'constant', "NA");

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
if ~isempty(s_dict)
    fprintf('Found info in showcase:\n')
    for jj = 1:size(s_dict, 1)
        fprintf('%d-', jj)
        for kk = 1:numel(info_dict_idx)
            fprintf('%s: %s\n', h_dict(info_dict_idx(kk)), s_dict{jj, info_dict_idx(kk)})
        end   
        fprintf('-------------------------------------------\n')
    end
end

if ~isempty(s_coding)
    fprintf('\nFound info in coding data:\n')
    for jj = 1:size(s_coding, 1)
        fprintf('%d-', jj)
        for kk = 1:numel(h_coding)
            fprintf('%s: %s\n', h_coding(kk), s_coding{jj, kk})
        end
        fprintf('-------------------------------------------\n')
    end
end

% Check showcase or coding approach -----------------------------------
if isempty(s_coding)
    check_dict = 1;
elseif isempty(s_dict)
    check_dict = 2;
else
    check_dict = [];
    while isempty(check_dict)
        check_dict = input('Continue according to 1-showcase 2-coding [1]? ');
        if isempty(check_dict) || ~ismember(check_dict, 1:2)
            check_dict = 1;
        end
    end
end

if check_dict == 2
    s_dict = s_coding;
    header_info = h_coding;
    info_dict_idx = 1:numel(header_info);
else
    header_info = h_dict;
end
% Show again options based on selected approach -----------------------
if ~isempty(s_coding) && ~isempty(s_dict)
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    for jj = 1:size(s_dict, 1)
        fprintf('%d-', jj)
        for kk = 1:numel(info_dict_idx)
            fprintf('%s: %s\n', header_info(info_dict_idx(kk)), s_dict{jj, info_dict_idx(kk)})
        end
        fprintf('------------------------------------------------------\n')
    end
end

check_term = [];
while isempty(check_term)
    if size(s_dict, 1) == 1
        check_term = 1;
        break
    end
    check_term = input('Insert term index(indices): ');
    if ~isempty(setdiff(check_term, 1:size(s_dict, 1)))
        check_term = [];
    end
end
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

% merge showcase-coding dictionaries --------------------------------------
s_dict = s_dict(check_term, :);
s_dict = s_dict.Variables;
coding_hdr_id = lower(header_info) == "coding";
if isscalar(unique(s_dict(:, coding_hdr_id)))
    s_dict_codinginfo = s_dict(1, coding_hdr_id);
else
    s_dict_codinginfo = s_dict(:, coding_hdr_id);
    fprintf(['\t\t*************\nWarning: >1 unique coding found!\nAlthough I ',...
        'generate a trait for you:\nthis step has not been widely tested!\n',...
        '\t\t*************\n'])
end

if check_dict == 2
    header_info = h_dict;
    coding_hdr_id_2 = lower(h_dict) == "coding";
    coding_hdr_id_2 = lower(data_dict_showcase{:, coding_hdr_id_2}) == s_dict_codinginfo;
    connected_info_all = data_dict_showcase{coding_hdr_id_2, :};
    included_info_dict = ["Path", "Field", "Coding", "Notes"];
    info_dict_idx = find(ismember(h_dict, included_info_dict));
else
    header_info = h_coding;
    coding_hdr_id_2 = lower(h_coding) == "coding";
    coding_hdr_id_2 = ismember(lower(codings_showcase{:, coding_hdr_id_2}), s_dict_codinginfo); 
    connected_info_all = codings_showcase{coding_hdr_id_2, :};
    info_dict_idx = 1:numel(h_coding);
end

% Display connected info to selected term
connected_info_single = ({});
addmore_flag = true;
ct = 1;
while addmore_flag
    if ~isempty(connected_info_all)
        fprintf('\n~~~~~~~~~~~~~~~~~~~ Connected info ~~~~~~~~~~~~~~~~~~~~~\n')
        keyword_sel = input('Do you want to search using keyword ([y]/n)? ', 's');
        if ~strcmpi(keyword_sel, 'n')
            keyword_sel = input('Insert key term(s). Use ; to separate multiple choices: ', 's');
            keyword_sel = string(strsplit(keyword_sel, ';'));
            [row_id, ~] = find(ismember(lower(connected_info_all), lower(string(keyword_sel))));
            connected_info = connected_info_all(unique(row_id), :);
        else
            connected_info = connected_info_all;
        end

        for jj = 1:size(connected_info, 1)
            fprintf('%d-', jj)
            for kk = 1:numel(info_dict_idx)
                fprintf('%s: %s\n', header_info(info_dict_idx(kk)), connected_info(jj, info_dict_idx(kk)))
            end
            fprintf('------------------------------------------------------\n')
        end
        check_term_flag = true;
        while check_term_flag
            if check_dict == 1
                check_term = input('Insert secondary term id(s) [press Enter to skipp]: ');
                if isempty(setdiff(check_term, 1:size(connected_info, 1))) || isempty(check_term)
                    check_term_flag = false;
                end
            else % showcase_dictionary data cannot be skipped!
                check_term = input('Insert term id(s): ');
                if isempty(setdiff(check_term, 1:size(connected_info, 1)))
                    check_term_flag = false;
                    if isempty(check_term)
                        check_term_flag = true;
                    end
                end
            end
        end

        if ~isempty(check_term)
            connected_info_single{ct} = connected_info(check_term, :); % Selected connected info
            if check_dict == 1 % Check if user wants to add more secondary terms
                addmore_sel = input('Add more secondary terms to create another trait using the same df (e.g. CLD and cirrhosis) (y/[n])?', 's');
                if strcmpi(addmore_sel, 'y')
                    addmore_flag = true;
                else
                    addmore_flag = false;
                end
            else
                addmore_flag = false;
            end
        else
            connected_info_single{ct} = [];
        end

    else
        connected_info_single{ct} = [];
        addmore_flag = false;
    end
    ct = ct + 1;
end

% Define input terms --------------------------------------------------

if check_dict == 1
    layer_1_data = s_dict;
    layer_1_h = h_dict;
    layer_2_data = connected_info_single;
    layer_2_h = h_coding;
else
    layer_2_data = s_dict;
    layer_2_h = h_coding;
    layer_1_data = connected_info_single{1};
    layer_1_h = h_dict;
end

% Extract additional dictionary for UKBtraitParser
coding_hdr_id = lower(layer_1_h) == "coding";
layer_1_coding = layer_1_data(:, coding_hdr_id);

coding_hdr_id = lower(layer_2_h) == "coding";
additionalDict_id = ismember(codings_showcase{:, coding_hdr_id}, unique(layer_1_coding));
meaning_hdr_ids = contains(lower(layer_2_h), ["meaning", "value"]);
inputTerms.additionalDict = codings_showcase{additionalDict_id, meaning_hdr_ids};

fieldID_idx = layer_1_h == "FieldID";
field_idx = layer_1_h == "Field";

tag_change = input('Do you want to rename the tag y/[n]? ', 's');
if strcmpi(tag_change, 'y')
    tag_change = input('Insert tag of intereset: ', 's');
    inputTerms.tags = cellstr(tag_change);
else
    inputTerms.tags = cellstr(join(layer_1_data(:, field_idx), "_"));
end
inputTerms.fieldID = {cellstr(layer_1_data(:, fieldID_idx))};

inputTerms.keyTerms = {cellstr(layer_1_data(:, field_idx))};

value_idx = layer_2_h == "Value";
meaning_idx = layer_2_h == "Meaning";

if ~isstring(layer_2_data)
    if ~isempty(layer_2_data{1})
        for jj = 1:numel(layer_2_data)
            layer_2_eachdata = layer_2_data{jj};
            inputTerms.secondaryTerms{jj, 1} = char(join(layer_2_eachdata(:, value_idx), ","));
            inputTerms.secondaryTermsMeaning{jj, 1} = char(join(layer_2_eachdata(:, meaning_idx), ","));
        end
    else
        inputTerms.secondaryTerms = {'NOSECONDARYTERMISSELECTED'};
        inputTerms.secondaryTermsMeaning = {'NOSECONDARYTERMISSELECTED'};
    end
else
    inputTerms.secondaryTerms{1} = char(join(layer_2_data(:, value_idx), ","));
    inputTerms.secondaryTermsMeaning{1} = char(join(layer_2_data(:, meaning_idx), ","));
end
inputTerms.secondaryTerms = {inputTerms.secondaryTerms};
inputTerms.secondaryTermsMeaning = {inputTerms.secondaryTermsMeaning};

end % END

%% ------------------------------------------------------------------------
function UKB_STRUCT_ALL = parsePhenoQuery(query, dict, codings, opts)
% searches for the query in data dictionary/codings. It's different from
% searchUKBdictionary in the sense that it doesn't merge different query
% terms into one (see above), but instead stores each match to a distinct

if opts.surv, opts.merge = true; end
if opts.surv && ~opts.save, opts.mergesave = false; end

% phenotype struct.
query = string(query);

if opts.merge && isfield(opts, "tag")
    opts.tag = string(opts.tag);
    % if opts.surv
    %     opts.tag = opts.tag(1);
    % else
        opts.tag = repmat(opts.tag(1), 1, numel(query));
    % end
    opts.mergesave = opts.save;
    opts.save = false; % don't save individual phenos
end

if opts.merge && isfield(opts, "desc")
    opts.desc = string(opts.desc);
    opts.desc = repmat(opts.desc(1), 1, numel(query));
end


% initialize pheno tags/exclusion criteria
fi = ["df", "tag", "remove", "coding", "desc"];
for i = 1:numel(fi)
    if isfield(opts, fi(i))
        if opts.surv && fi(i) == "remove" % check for struct
            ff = string(fieldnames(opts.(fi(i))));
            for k = 1:numel(ff)
                tmp = opts.(fi(i)).(ff(k)); % if struct array
                misstagidx = ismissing(tmp) | (tmp == "");
                tmp(misstagidx) = missing;
                opts.(fi(i)).(ff(k)) = tmp.split([",", ";"]);
            end

            % if size(opts.(fi(i)), 2) ~= numel(query)
            %     error('size of %s must be equal to the size of query!', fi(i))
            % end
        else

            opts.(fi(i)) = string(opts.(fi(i)));
            misstagidx = ismissing(opts.(fi(i))) | (opts.(fi(i)) == "");
            opts.(fi(i))(misstagidx) = missing;
            if numel(opts.(fi(i))) ~= numel(query)
                error('size of %s must be equal to the size of query!', fi(i))
            end
        end
    elseif opts.surv && fi(i) == "remove" 
        opts.(fi(i)) = struct;
    else
        opts.(fi(i)) = strings(numel(query), 1);
        opts.(fi(i))(:) = missing;
    end
end

removecodes = opts.remove;
opts.remove2 = cell(numel(opts.remove), 1); %@12JUNE2024: used only when 'surv' is false
UKB_STRUCT_ALL = cell(numel(query), 1);
for i = 1:numel(query)
    opts.dfonly = false;
    
    opts.mortality = false; %@28MAY2024: overall mortality
    % is query is missed, fetch data-fields in 'df'
    if query(i).lower == "mortality"
        opts.mortality = true; % coding: 19
    elseif isempty(query(i)) || ismissing(query(i)) || query(i) == ""
        if ismissing(opts.df(i))
            fprintf('%d (of %d)-missing query/df\n', i, numel(query))
            continue
        else
            opts.dfonly = true; % ignore coding values
            query(i) = opts.df(i);
        end
    end

    fprintf('%d (of %d)-query: %s\n', i, numel(query), query(i))
    txt = split(query(i), [",", ";"]).lower;
    
    % @22FEB2024: 'remove' can be a struct with different df/cd fields. In
     % this case we treat it independently of the input query
    if isstruct(opts.remove)
        
        % each struct per each query, with the exception of 'surve' == true
        if opts.surv
            removecodes = opts.remove; 
        else
            removecodes = opts.remove(i);
        end
        
        if ~isempty(fieldnames(removecodes))
            % appened coding/df/instances for exclusion criteria to
            % opts.remove struct
            ff = setdiff(string(fieldnames(removecodes)), ["meta", "metac", "instance"]);
            re_df_cd = arrayfun(@(x)split(x, "_"), ff, uni=false);
            % [opts.remove, opts.remove_instance] = deal(strings(numel(ff), 1));
            % idxcr = cell(numel(ff), 1); % a cell of containing matching indices between opts.remove fields and codings/value
            for k = 1:numel(re_df_cd)
                re_cd_df_val = extractAfter(re_df_cd{k}(1), "df"|"cd");
                if ff(k).startsWith("df")
                    removecodes.meta{k, 1} = dict(dict.FieldID == re_cd_df_val, :);
                else
                    removecodes.meta{k, 1} = dict(dict.Coding == re_cd_df_val, :);
                end
    
                % include coding table for each field
                removecodes.metac{k, 1} = codings(ismember(codings.Coding, removecodes.meta{k}.Coding) & ismember(codings.Value, removecodes.(ff(k))), :);
        
                if numel(re_df_cd{k}) > 1 % instance is provided
                    removecodes.instance(k) = re_df_cd{k}(2);
                else % use all instances for exclusion criteria of ff(k)
                    removecodes.instance(k) = ""; 
                end
            end
            removecodes.instance(removecodes.instance == "") = missing;
    
            if opts.surv
                opts.remove = removecodes; 
            else
                opts.remove2{i} = removecodes;
            end
        end

    else % text vector
        opts.remove2{i} = split(removecodes(i), [";", ","]).lower;
    end


    % look for field-ids/descriptions or coding meanings
    if all(~isnan(double(txt))) % field id or coding value
        idxd = ismember(dict.FieldID, txt);
        idxc = ismember(lower(codings.Value), txt);
        if ~any(idxc), idxc = contains(lower(codings.Value), txt); end

        if ~opts.surv && ~ismissing(opts.remove2{i}) && ~isstruct(opts.remove2{i}) % codes to be removed: survival treats them differently
            idxcr = ismember(lower(codings.Value), opts.remove2{i});
            if ~any(idxcr)
                idxcr = contains(lower(codings.Value), opts.remove2{i});
            end
        else
            idxcr = false(height(codings), 1);
        end
        

        % idxcr should have the same coding as idxc (e.g. 1065 for
        % hypertension in df 20002). @22FEB2024: this is true only when
        % 'remove' is a text vector
        idxc = idxc | idxcr;
    else % text

        if opts.mortality
            idxd = dict.Coding == "19";
            [idxc, idxcr] = deal(false(height(codings), 1));

        else
            idxd = ismember(lower(dict.Field), txt);
            if ~any(idxd), idxd = contains(lower(dict.Field), txt); end
            idxc = ismember(lower(codings.Meaning), txt) | ismember(lower(codings.Value), txt);
            if ~any(idxc)
                idxc = contains(lower(codings.Meaning), txt) | contains(lower(codings.Value), txt);
            end
    
            if ~opts.surv && ~all(ismissing(opts.remove2{i})) && ~isstruct(opts.remove2{i}) % codes to be removed: survival treats them differently
                idxcr = ismember(lower(codings.Meaning), opts.remove2{i}) |...
                    ismember(lower(codings.Value), opts.remove2{i});
                if ~any(idxcr)
                    idxcr = contains(lower(codings.Meaning), opts.remove2{i}) |...
                        contains(lower(codings.Value), opts.remove2{i});
                end
            else
                idxcr = false(height(codings), 1);
            end
            % idxcr should ideally have the same coding as idxc (e.g.
            % ICD-10 codes) @22FEB2024: this is true only when
            % 'remove' is a text vector
            idxc = idxc | idxcr;
        end
    end
    
    [tdic, tcod] = deal([]);
    if any(idxd); tdic = dict(idxd, :); end
    if any(idxc); tcod = codings(idxc, :); end
    if any(idxcr); tcodr = codings(idxcr, :); else; tcodr.Value = missing; end
    
    if ~isempty(tdic) && ~isempty(tcod) % found query in both (i.e. description and meaning)
        if ~ismissing(opts.coding(i))
            tdic(~ismember(tdic.Coding, opts.coding(i)), :) = [];
            tcod(~ismember(tcod.Coding, opts.coding(i)), :) = [];
        end

        if ismissing(opts.df(i)) && ~ismissing(opts.coding(i))
            % @22APR2025: if 'df' field is empty and 'coding' is used, it
            % naturally means coding is preferred over 'df', so we clear
            % tdic
            fprintf("\tboth dictionary/coding are non-empty, chose coding since 'df' is empty\n")
            tdic = [];
            
        elseif ~ismissing(opts.df(i)) % only keep selected dfs
            tdic(~ismember(tdic.FieldID, split(opts.df(i), ';')), :) = [];
        end
        
        if isempty(tdic)
            if ~ismissing(opts.coding(i))
                tcod(~ismember(tcod.Coding, opts.coding(i)), :) = [];
            end
             opts.categorical = 1; % temp
             tdic = dict(ismember(dict.Coding, tcod.Coding), :);
        elseif opts.dfonly % only dfs (no query, no codings)
            opts.categorical = [];
        else
            format("compact")
            disp("found in data dictionary:")
            disp(tdic(:, ["Field", "FieldID", "Coding"]))
            fprintf('\n')
            disp("found in data codings:")
            disp(tcod)
            sel = input("[1]dictionary or [2]codings [1]? ");
            if isempty(sel) || sel ~= 2
                opts.categorical = [];
            else
                opts.categorical = 1; % temp
                tdic = dict(ismember(dict.Coding, tcod.Coding), :);
            end
            format("default")
        end

    elseif ~isempty(tdic)
        if ~ismissing(opts.coding(i))
            tdic(~ismember(tdic.Coding, opts.coding(i)), :) = [];
        end
        opts.categorical = [];
        
    elseif ~isempty(tcod)
        if ~ismissing(opts.coding(i))
            tcod(~ismember(tcod.Coding, opts.coding(i)), :) = [];
        end

        % find Coding within dict
        opts.categorical = 1; % temp
        tdic = dict(ismember(dict.Coding, tcod.Coding), :);

    else
        fprintf('term(s): %s not found\n', join(txt, ";"))
        continue
    end

    if ~ismissing(opts.df(i)) % only keep these dfs
        tdic(~ismember(tdic.FieldID, split(opts.df(i), ';')), :) = [];
    end

    if isempty(tdic) % codings don't correspond to any df
        fprintf('input codings do not correspond to any data-fields!\n')
        continue
    end

    % check tdic
    if height(tdic) > 1 && ~opts.all
        tdic = sortrows(tdic, 'Coding');
        showtdic = tdic(:, ["Field", "FieldID", "Coding"]);
        showtdic.index = (1:height(showtdic))';
        showtdic = movevars(showtdic, 'index', 'Before', 1);
        disp(showtdic)
        sel = input("index(indices) [enter to select all]? ");
        if isempty(sel) || any(sel > height(tdic)) || any(sel < 0)
            sel = 1:height(tdic);
        end
        tdic = tdic(sel, :);
    end

    if ~isempty(opts.categorical)
        tcod(~ismember(tcod.Coding, tdic.Coding), :) = [];
        opts.categorical = unique(tcod.Value); % categorical data have codings only
        if opts.surv
            opts.categorical_surv{i} = opts.categorical;
        end
        if all(~ismissing(tcodr.Value)) && all(~ismissing(opts.remove2{i}))
            opts.remove2{i} = intersect(opts.categorical, tcodr.Value);
        end
    end
    
    % @22FEB2024: add date corresponding to each data-field: note that many
    % of these dates for diagnosis and death data are now merged in HESIN
    % and DEATH data tables from UKBB.
    if opts.surv
        % this info can be found at: https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=141140
        % 20008 is interpolated year of non-cancer diagnosis
        % 20006 same as above for cancer self reported
        d2t = struct;
        d2t.df = ["41202", "41203", "41270", "41271", "20002", "20001", "40006", "40001"];
        d2t.date = ["41262", "41263", "41280", "41281", "20008", "20006", "40005", "40000"];

        tdic = getTimeDatedfs(tdic, d2t, dict, opts);

        % do the same for 'remove' fields
        if ~isempty(fieldnames(opts.remove))
            for k = 1:numel(opts.remove.meta)
                opts.remove.meta{k} = getTimeDatedfs(opts.remove.meta{k}, d2t, dict, opts);
                opts.remove.meta{k}.fnc = strings(height(opts.remove.meta{k}), 1);
                opts.remove.meta{k}.fnc(ismember(lower(opts.remove.meta{k}.ValueType), ["date", "time"])) = "datetime";
                opts.remove.meta{k}.fnc(ismember(lower(opts.remove.meta{k}.ValueType), ["continuous", "integer"])) = "double";
            end
        end

        % also fetch the date leaving follow-up
        if ~isfield(opts, "dlf")
            opts.dlf = table2struct(phenoParser("query", "", "df", "191", ...
                "instance", 0, ...
                "save", false, "surv", false, "ukbrap", opts.ukbrap));
        end

        % also fetch date of attending
        if ~isfield(opts, "regdate")
            opts.regdate = table2struct(phenoParser("query", "", "df", "53", ...
                "save", false, "ukbrap", opts.ukbrap, ...
                "instance", 0, "surv", false, "all", true));
        end

        % feth date of birth
        if ~isfield(opts, "dob")
            % read date of birht based on year and month of birth
                yob = table2struct(phenoParser("query", "", "df", "34", "save", false, "surv", false, "ukbrap", opts.ukbrap));
                mob = table2struct(phenoParser("query", "", "df", "52", "save", false, "surv", false, "ukbrap", opts.ukbrap));
                [f1, f2] = ismember(yob.eid, mob.eid); f2(f2<1) = [];
                opts.dob = table(yob.eid(f1), datetime(yob.rawUKB(f1), double(mob.rawUKB(f2)), ones(sum(f1), 1), ...
                    Format="dd/MM/uuuu"), VariableNames=["eid", "dob"]);
                clear yob mob f1 f2
        end

    else
        tdic.dateField = cell(height(tdic), 1);

    end
    
    % append data transformation funcs
    tdic.fnc = strings(height(tdic), 1);
    tdic.fnc(ismember(lower(tdic.ValueType), ["date", "time"])) = "datetime";
    tdic.fnc(ismember(lower(tdic.ValueType), ["continuous", "integer"])) = "double";

    % fetch data from UKB-RAP
    if opts.ukbrap

        % add 'remove' data-fields.
        if isfield(opts, "remove") && isfield(opts.remove, "meta")
            tmp_tdic = [tdic; vertcat(opts.remove.meta{:})];
        else
            tmp_tdic = tdic;
        end

        dx_extract_dataset(tmp_tdic);
    end

    % load tdic df(s) from latest baskets
    if opts.surv
        UKB_STRUCT_ALL{i} = tdic; % (:, ["FieldID", "ValueType", "fnc", "dateField"]);
        continue
    end

    opts.ct = i;

    basket = loadDFfromUKBBasket(tdic, opts);

    if isempty(basket)
        error('phenoParser:noEID', 'query has no observations!')
    end

    if ~isempty(opts.categorical)
        % only when query is found within codings table
        basket.meaning = tcod(ismember(tcod.Value, basket.codes), :);
        if ~ismissing(opts.tag(i))
            basket.tag = opts.tag(i);
        end
        try basket.desc = opts.desc(i); catch, end
        if isempty(basket.eid)
            disp("query codes could not be found within the basket(s)")
            continue
        end
        UKB_STRUCT_ALL{i, 1} = writePheno2Mat(basket, opts);
        fprintf('Done.\n\n')
        continue
    end
    
    % for non-categorical variables (only when query is not found within
    % codings table)
    % re-arrange returned df(s) into pheno struct
    outpheno_basket = cell(height(basket), 1);
    for j = 1:height(basket)
        if isempty(basket.chunk{j})
            continue
        end
        phenoindex.fi = string(setdiff(fieldnames(basket.chunk{j}), 'eid'));

        if opts.ukbrap
            phenoindex.df = tdic.FieldID(ismember(tdic.name, phenoindex.fi));
            [~, idx] = ismember(phenoindex.fi, tdic.name);
        else
            phenoindex.df = string(extractBetween(phenoindex.fi, 'x', '_'));
            [~, idx] = ismember(phenoindex.df, tdic.FieldID);
        end

        phenoindex.raw = basket.chunk{j};
        basket.chunk{j} = [];
        phenoindex.basket = string(basket.name);
        try phenoindex.desc = opts.desc(i); catch, end
        if height(basket) == 1 && isscalar(unique(phenoindex.df)) && ~ismissing(opts.tag(i))
            phenoindex.tag = repmat(opts.tag(i), numel(idx), 1);
        else
            phenoindex.tag = tdic.Field(idx);
        end
        phenoindex.coding = tdic.Coding(idx);
        phenoindex.numeric = true(numel(phenoindex.fi), 1);
        
        if opts.ukbrap
            [~, idx] = ismember(phenoindex.fi, tdic.name);
        else
            [~, idx] = ismember(phenoindex.df, tdic.FieldID);
        end
        phenoindex.numeric(tdic.fnc(idx) == "") = false;

        if any(colnames(tdic) == "Units")
            phenoindex.unit = tdic.Units(idx);
        end

        % get term meanings for categorical traits
        tmpcodes = unique(phenoindex.coding);
        tmpcodes(tmpcodes == "" | ismissing(tmpcodes)) = [];
        phenoindex.meaning = cell(numel(phenoindex.fi), 1);
        if ~isempty(tmpcodes)
            
            for k = 1:numel(tmpcodes)
                idx = codings.Coding == tmpcodes(k);
                tmpterms = codings(idx, :);
                
                idx = phenoindex.coding == tmpcodes(k);
                phenoindex.meaning(idx) = {tmpterms};
            end
        end
        
        % save phenotypes
        outpheno_basket{j} = writePheno2Mat(phenoindex, opts);
    end
    
    UKB_STRUCT_ALL{i, 1} = vertcat(outpheno_basket{:});
    fprintf('Done.\n\n')
end


if opts.surv
    % add codes (inclusion criteria)
    if opts.mortality
        tdic = vertcat(UKB_STRUCT_ALL{:});
        tdic.codes(:) = "";
    else
        for k = 1:numel(UKB_STRUCT_ALL)
            UKB_STRUCT_ALL{k}.codes(:) = join(opts.categorical_surv{k}, ",");
        end
        tdic = vertcat(UKB_STRUCT_ALL{:});
        opts.categorical = unique(vertcat(opts.categorical_surv{:}));
    end
    
    % fetch date-fields if necessary
    if opts.ukbrap
        dt_tmp1 = tdic.dateField;
        dt_tmp1(cellfun(@isempty, dt_tmp1)) = [];
        dt_tmp1 = vertcat(dt_tmp1{:});

        % add 'remove' data-fields.
        if isfield(opts, "remove") && isfield(opts.remove, "meta")
            dt_rem_tmp = vertcat(opts.remove.meta{:});
            dt_rem_tmp = dt_rem_tmp.dateField;
            dt_rem_tmp(cellfun(@isempty, dt_rem_tmp)) = [];
            dt_rem_tmp = vertcat(dt_rem_tmp{:});
            
            dt_tmp2 = [dt_tmp1; dt_rem_tmp];
        else
            dt_tmp2 = dt_tmp1;
        end
        
        dt_tmp2 = convertvars(dt_tmp2, ["array", "instance"], @string);
        dt_tmp2 = fillmissing(dt_tmp2, "constant", "", ...
            DataVariables=["primary_key_type", ...
            "is_sparse_coding", "linkout", "Units", "array", "instance"]);
        dt_tmp2 = unique(dt_tmp2, "rows");
        assert(isempty(duplicates(dt_tmp2.name)), "date-field names cannot have duplicated names!")

        dx_extract_dataset(dt_tmp2);
    end
    
    [UKB_STRUCT_ALL, df] = loadDFfromUKBBasket(tdic, opts);
    
    assert(all(UKB_STRUCT_ALL.tt(UKB_STRUCT_ALL.tt0 < 0) < 0), "Some age-of-onset are undefined while they're defined for time-to-event!")

    %@07JUNE2024: remove those with an age-of-onset < 0 (nonsensical)
    UKB_STRUCT_ALL(UKB_STRUCT_ALL.tt0 < 0, :) = [];

    %@10JUNE2024: for mortality, individuals died before enrollment
    %(nonsensical) should be removed
    if opts.mortality
        UKB_STRUCT_ALL(UKB_STRUCT_ALL.tt < 0, :) = [];
    end

    % idx = UKB_STRUCT_ALL.tt < 0; % those before the start up: see fitSurvCox
    % if any(idx)
    %     before_startup = UKB_STRUCT_ALL(idx, :);
    %     UKB_STRUCT_ALL(idx, :) = [];
    %     ex_eid_before = before_startup.eid;
    %     fprintf("\t%d individuals with a diagnoses before start-up date were excluded (but available in 'ex' field).\n\n", numel(ex_eid_before))
    %     % before_startup = groupsummary(before_startup, ["ex", "censor", "died", "df"]); % for descriptive purposes
    % else
    %     ex_eid_before = [];
    %     before_startup = table;
    % end

    % baseline individuals, remove if necessary: see fitSurvCox
    % @13JUNE2024: not a precise approach, because a tt still can be < 0.
    % This is handled properly in 'fitSurvCox' function though.
    df_baseline = df(df.instance == "0" & ~df.inc & df.codes ~= "", :);
    UKB_STRUCT_ALL.base(:) = false;
    if ~isempty(df_baseline)
        baseline_codes = arrayfun(@(x)split(x, ","), df_baseline.codes, uni=false);
        baseline_codes = unique(vertcat(baseline_codes{:}));
        idx_base = find(~ismissing(UKB_STRUCT_ALL.rawUKB));
        rawcodes = arrayfun(@(x)split(x, ","), UKB_STRUCT_ALL.rawUKB(idx_base), uni=false);
        idx_base_true = cellfun(@(x)any(ismember(x, baseline_codes)), rawcodes);
        idx_base(~idx_base_true) = [];
        
        UKB_STRUCT_ALL.base(idx_base) = true;
        % base_tab = UKB_STRUCT_ALL(idx_base, :); % for descriptive purposes
        % UKB_STRUCT_ALL(idx_base, :) = [];
        fprintf("\t%d individuals with the pre-existing codes at baseline (available in 'base' field).\n\n", numel(idx_base))
        % ex_eid_baseline = base_tab.eid;
        % base_tab = groupsummary(base_tab, ["ex", "censor", "died", "df"]);
    % else
    %     base_tab = table;
    %     ex_eid_baseline = [];
    end
    
    %25MAY2024: deprecated
    % % check for landmark cut-off point
    % UKB_STRUCT_ALL.landmark(:) = false; % for downstream sensitivity analyses
    % if opts.landmark > 0
    %     idx = UKB_STRUCT_ALL.censor == 1 & UKB_STRUCT_ALL.tt <= opts.landmark;
    %     UKB_STRUCT_ALL.landmark(idx) = true;
    %     if any(idx)
    %         fprintf("\t%d individuals have a diagnosis before landmark time-point (%.2f yr)." + ...
    %             "\n\tNote that they're not excluded and are useful for sensitivity analyses.\n\n", sum(idx), opts.landmark);
    %     end
    % end
    
    % prepare a struct for writing
    UKB_STRUCT_ALL(:, ["ex", "died"]) = []; % not needed anymore. Downstream analyses are based on censor field
    UKB_STRUCT_ALL = renamevars(UKB_STRUCT_ALL, "df", "source");
    cen_date = UKB_STRUCT_ALL.Properties.UserData;
    UKB_STRUCT_ALL.Properties.UserData = [];
    UKB_STRUCT_ALL = table2struct(UKB_STRUCT_ALL, "ToScalar", true);
    UKB_STRUCT_ALL.info = struct;
    UKB_STRUCT_ALL.info.basket = string(join(extractAfter(df.Properties.UserData.name, "UKBFileParser_"), ","));
    UKB_STRUCT_ALL.info.date = datetime("now");
    UKB_STRUCT_ALL.info.desc = opts.desc;
    UKB_STRUCT_ALL.info.help = table(["0", "1", "2", "3", "tt", "tt0"]', ...
        ["after censoring date", "event", "death before event", ...
        "competing risk before event", "time-to-event", "age-of-onset"]', ...
        VariableNames=["code", "meaning"]);
    UKB_STRUCT_ALL.info.cen_date = cen_date;
    UKB_STRUCT_ALL.info.df = join(unique(df.df), ",");

    if ~isfield(opts, "tag") && (opts.tag == "" || ismissing(opts.tag))
        opts.tag = replace(getRandomName("pheno", ".", "_"));
    end

    if opts.mortality && ismissing(opts.tag), opts.tag = "Mortality"; end
    UKB_STRUCT_ALL.tag = opts.tag;
    UKB_STRUCT_ALL.numericFlag = false;
    
    % % add excluding eids
    % ex_eid = union(ex_eid_baseline, ex_eid_before);
    % if ~isempty(ex_eid)
    %     UKB_STRUCT_ALL.exeid = ex_eid;
    %     UKB_STRUCT_ALL.ex = struct;
    %     UKB_STRUCT_ALL.ex.before_startup = before_startup;
    %     UKB_STRUCT_ALL.ex.before_startup_eid = ex_eid_before;
    %     UKB_STRUCT_ALL.ex.baseline = base_tab;
    %     UKB_STRUCT_ALL.ex.baseline_eid = ex_eid_baseline;
    % 
    % end

    % write to a mat file
    if opts.save || opts.mergesave
        % for a lighter version, we write less usefull fields: 'rawUKB',
        % 'excode', 'source' to metadata saved together with pheno
        meta = struct;
        meta.rawUKB = UKB_STRUCT_ALL.rawUKB;
        meta.excode = UKB_STRUCT_ALL.excode;
        meta.source = UKB_STRUCT_ALL.source;
        UKB_STRUCT_ALL = rmfield(UKB_STRUCT_ALL, ["rawUKB", "excode", "source"]);
        out_name = matlab.lang.makeValidName(UKB_STRUCT_ALL.tag) + ".surv.mat";
        save(fullfile(opts.phenodir, out_name), "UKB_STRUCT_ALL", "meta")
        UKB_STRUCT_ALL.meta = meta;
        clear meta
    end

    return
end

% merge all phenos into a single one. This should be done for categorical
% traits
if opts.merge
    if isfield(opts, "mergesave")
        opts.save = opts.mergesave;
    end
    UKB_STRUCT_ALL = mergePhenos(UKB_STRUCT_ALL, opts);
else
    try UKB_STRUCT_ALL = vertcat(UKB_STRUCT_ALL{:}); catch; end
end

end % END

%% ------------------------------------------------------------------------
function [basket, df] = loadDFfromUKBBasket(df, opts)
% first finds latest basked id that contains input data fields (df), and
% then loads df(s) themselves.

df = renamevars(df, ["FieldID", "ValueType", "fnc"], ["df", "vt", "fnc"]); % data-field, variable type

if ~opts.ukbrap
    [~, idx] = unique(df.df, "stable");
    df = df(idx, :);
end

if opts.surv
    df = renamevars(df, "dateField", "dt"); % date data field associated with 'df'
    
    dt = df.dt; df.dt = [];

    % merge df with dt into a single df
    tmp_dt = cell(numel(dt), 1);
    for k = 1:numel(dt)
        if istable(dt{k})
            tmp_dt{k} = dt{k};
            tmp_dt{k}.Field = [];
            if opts.ukbrap
                tmp_dt{k}.relation(:) = df.name(k);
            else
                tmp_dt{k}.relation(:) = df.df(k);
            end
        end
    end
    tmp_dt(cellfun(@isempty, tmp_dt)) = [];
    tmp_dt = vertcat(tmp_dt{:});
    if ~isempty(tmp_dt)
        tmp_dt = renamevars(tmp_dt, ["FieldID", "ValueType"], ["df", "vt"]);
    end
    tmp_dt.codes(:) = "";
    df.relation(:) = "";
    df = [df(:, colnames(tmp_dt)); tmp_dt];
    % df.codes(:) = "";
    df.inc(:) = true; % inclusion
    
    if any(colnames(df) == "instance")
        df.instance = string(df.instance);
        df.ii = df.instance; % for future use: original instance
    end

    if isempty(opts.instance)
        df.instance(:) = "";
    else
        df.instance(:) = string(opts.instance);
    end

    % check if exclusion criteria are also provided
    if ~isempty(fieldnames(opts.remove))
        rst = opts.remove;
        [rdf, rdt] = deal(cell(numel(rst.meta), 1));
        for k = 1:numel(rst.meta)
            tmp = rst.meta{k};

            if opts.ukbrap
                tmp = renamevars(tmp, ["FieldID", "ValueType", "fnc", "dateField"], ...
                    ["df", "vt", "fnc", "dt"]);
            else
                tmp = tmp(:, ["FieldID", "ValueType", "fnc", "dateField"]);
                tmp.Properties.VariableNames = ["df", "vt", "fnc", "dt"];
            end

            if any(colnames(tmp) == "instance")
                tmp.instance = string(tmp.instance);
                tmp.ii = tmp.instance;
            end

            tmp.codes(:) = rst.metac{k}.Value.join(",");
            tmp.relation(:) = "";
            tmp.instance(:) = rst.instance(k);

            if opts.ukbrap
                rdf{k} = tmp;
            else
                rdf{k} = tmp(:, ["df", "vt", "fnc", "relation", "codes", "instance"]);
            end

            % date-time table for this field
            dt = tmp.dt;
            tmp_dt = cell(numel(dt), 1);
            for j = 1:numel(dt)
                if istable(dt{j})
                    tmp_dt{j} = dt{j};
                    if any(colnames(tmp_dt{j}) == "instance")
                        tmp_dt{j}.instance = string(tmp_dt{j}.instance);
                        tmp_dt{j}.ii = tmp_dt{j}.instance;
                    end
                    tmp_dt{j}.Field = [];

                    if opts.ukbrap
                        tmp_dt{j}.relation(:) = rdf{k}.name(j);
                    else
                        tmp_dt{j}.relation(:) = rdf{k}.df(j);
                    end
                    tmp_dt{j}.instance(:) = rst.instance(k);
                end
            end

            tmp_dt(cellfun(@isempty, tmp_dt)) = [];
            tmp_dt = vertcat(tmp_dt{:});
            if ~isempty(tmp_dt)
                tmp_dt = renamevars(tmp_dt, ["FieldID", "ValueType"], ["df", "vt"]);
            end
            rdt{k} = tmp_dt;
        end

        rdt(cellfun(@isempty, rdt)) = [];
        rdt = vertcat(rdt{:});
        rdf = vertcat(rdf{:});

        if ~isempty(rdt)
            rdt.codes(:) = "";
            keep_cols = intersect(colnames(rdt), colnames(rdf));
            rdf = [rdf(:, keep_cols); rdt(:, keep_cols)];
        end
        
        rdf.inc(:) = false; % exclusion criteria
        keep_cols = intersect(colnames(rdf), colnames(df));
        df = [df(:, keep_cols); rdf(:, keep_cols)];
    end

    df.instance(df.instance == "") = missing;
    
end
    
basket = struct2table(dir(fullfile(opts.basketdir, "UKBFileParser_*")));
basket(~basket.isdir, :) = [];

if isempty(basket)
    fprintf('no basket was found in %s\n', opts.basketdir)
    return
end

if opts.ukbrap
    % remove other static baskets
    idx = contains(basket.name, "_ukb" + digitsPattern) & ~ismember(basket.name, opts.ukbrap_dir);
    basket(idx, :) = [];
else
    basket(ismember(basket.name, opts.ukbrap_dir), :) = [];
end

basket.folder = fullfile(basket.folder, basket.name);
basket.date = NaT(height(basket), 1); % date each basket was created
basket.dfcoverage = zeros(height(basket), 1);


for i = 1:height(basket)
    vmapper_file = fullfile(basket.folder{i}, 'variableMapper.mat');
    if ~isfile(vmapper_file)
        basket.empty(i, 1) = true;
        continue
    else
        basket.empty(i, 1) = false;
    end

    index = load(vmapper_file);

    %@21SEP2024
    if numel(index.datev) > 1 % UKB-RAP
        basket.date(i, 1) = max(index.datev);
        basket.idx{i, 1} = intersect(df.name, index.variableMapper(:, 2));
    else
        basket.date(i, 1) = index.datev;
        index = string(index.variableMapper(:, 2));
        index = unique(extractBetween(index, "x", "_"));
        basket.idx{i, 1} = intersect(df.df, index);
        basket.dfcoverage(i) = numel(basket.idx{i, 1})/numel(unique(df.df)); % ratio of df(s) found in this basket
    end
end

basket(basket.empty, :) = [];
if isempty(basket)
    disp("all baskets are empty!")
    return
end

if ~opts.ukbrap
    % keep only baskets containing query data-fields
    basket(~basket.dfcoverage, :) = [];
end

if isempty(basket)
    fprintf('query data-fields cannot be found in neither of available baskets!\n')
    return
end

basket = sortrows(basket, 'date', 'descend');

basket.rem = false(height(basket), 1);
for i = 1:height(basket)-1
    if basket.rem(i); continue; end
    for k1 = i+1:height(basket)
        if basket.rem(k1); continue; end
        % what data-fields will remain after excluding those in the
        % most updated basket?
        rem_dfs = setdiff(basket.idx{k1}, basket.idx{i});

        % if there is nothing new in an older basket. However, if both
        % created on the same date, both should be checked.
        if isempty(rem_dfs) && (basket.date(i) > basket.date(k1))
            basket.rem(k1) = true;
        elseif basket.date(i) > basket.date(k1)
            basket.idx{k1} = rem_dfs;
        end
    end
end

basket(basket.rem, :) = [];
basket(cellfun(@isempty, basket.idx), :) = [];

if opts.mortality
    basket = basket(endsWith(basket.name, "_HESIN_DEATH"), :);
    if isempty(basket)
        fprintf('There is no basket with HESIN_DEATH tag!\n')
        return
    end
end

if opts.surv
    df.Properties.UserData = basket;
end

% load df(s) per basket
basket.chunk = cell(height(basket), 1);
for i = 1:height(basket) % loop over each basket
    fprintf("%d (of %d) reading basket: %s", i, height(basket), string(basket.name{i}))
    
     % subset of df overlapping the basket df idx
    if opts.ukbrap && ismember(opts.ukbrap_dir, basket.name{i})
        df_in = df(ismember(df.name, basket.idx{i}), :);
    else
        df_in = df(ismember(df.df, basket.idx{i}), :);
    end

    if endsWith(basket.name{i}, 'HESIN_DEATH')
        % HESIN_DEATH data are processed differently
        indata.fieldID{1} = basket.idx{i};

        if opts.surv
            % check if there are exclusion criteria: note that for
            % HESIN/death data these exclusion codes should be later
            % treated carefully, i.e. if they are competing cause of events
            % (in this case Fine-Gray models are more appropriate). Here,
            % we consider that competing risk are independent of the main
            % event (not optimal). This can be later further tunned when
            % conducting the survival analysis. Also note that, df_in
            % contains ICD10 codes, and df_in_ex also contains ICD10 codes
            % (we ignore ICD9 therefore), so it's safe to merge/join them
            % together (e.g. cancer + viral hepatitis codes), which in the
            % end should be censored (or perform a competing risk analysis).
            % 
            % Note that we ignore date data-fields (relation ~= "") since
            % HESIN data already contains this info.
            df_in_ex = df_in(~df_in.inc, :);
            if ~isempty(df_in_ex)
                ex_codes = df_in_ex.codes(df_in_ex.codes ~= "" & ~ismissing(df_in_ex.codes));
                ex_codes = arrayfun(@(x)split(x, ","), ex_codes, uni=false);
                ex_codes = join(unique(vertcat(ex_codes{:})), ",");
            else
                ex_codes = "";
            end
            
            inc_codes = df_in.codes(df_in.inc);
            inc_codes(inc_codes == "" | ismissing(inc_codes)) = [];
            if isempty(inc_codes) % only exclusion criteria provided for HESIN/DEATH?
                opts.censor = true; % change censoring to 
                inc_codes = ""; % ex_codes;
            else
                inc_codes = arrayfun(@(x)split(x, ","), inc_codes, uni=false);
                inc_codes = join(unique(vertcat(inc_codes{:})), ",");
                opts.censor = false;
            end
            indata.secondaryTerms{1} = inc_codes;

            % if isempty(opts.categorical), opts.categorical = ""; end
            % indata.secondaryTerms{1} = join(opts.categorical, ",");
            % if all(opts.categorical == ""), opts.categorical = []; end
            basket.chunk{i} = getHESIN_DEATH(indata, ...
                'basketdir', basket.folder{i}, 'meaning', false, ...
                'verbose', false, 'surv', true, 'ex', ex_codes, ...
                "dlf", opts.dlf, "mortality", opts.mortality);

        else % case/control
            if all(~ismissing(opts.remove2{opts.ct}))
                % current implementation of getHESIN_DEATH keeps only one code
                % per each eid, so two calls should be made, one with exclusion
                % and one with inclusion criteria. 
                indata.secondaryTerms{1} = join(setdiff(opts.categorical, opts.remove2{opts.ct}), ",");
            else
                if isempty(opts.categorical), opts.categorical = ""; end
                indata.secondaryTerms{1} = join(opts.categorical, ",");
                if all(opts.categorical == ""), opts.categorical = []; end
            end
            basket.chunk{i} = getHESIN_DEATH(indata, ...
                'basketdir', basket.folder{i}, 'meaning', false, 'verbose', false);

            if all(~ismissing(opts.remove2{opts.ct}))
                indata.secondaryTerms{1} = join(opts.remove2{opts.ct}, ",");
                remchunk = getHESIN_DEATH(indata, 'basketdir', ...
                    basket.folder{i}, 'meaning', false, 'verbose', false);
                idx1 = ismember(basket.chunk{i}.eid, remchunk.eid);
                basket.chunk{i}.eid(idx1) = [];
                basket.chunk{i}.rawUKB(idx1) = [];
                basket.chunk{i}.source(idx1) = [];
                basket.chunk{i}.ex.eid = remchunk.eid;
                basket.chunk{i}.ex.rawUKB = remchunk.rawUKB;
                basket.chunk{i}.ex.source = remchunk.source;
            end
        end

        fprintf('\n')
        continue
    end

    % load eid 
    if ~opts.ukbrap % UKB-RAP chunks have their own eids
        eid = load(fullfile(basket.folder{i}, 'UKB_eid.mat')).UKB_eid;
    else
        if opts.surv
            % remove death/HESIN data: should be handeled by HESIN, otherwise
            % censoring would be wrong
            df_in(ismember(df_in.df, ["40000", "40001", "40002"]), :) = [];
            basket.idx{i} = intersect(basket.idx{i}, df_in.name);
        end
    end

    index = string(load(fullfile(basket.folder{i}, 'variableMapper.mat')).variableMapper);

    if opts.ukbrap % dx_extract_dataset and UKBBasketParser use 'name' field
        index(:, 3) = index(:, 2);
    else
        index(:, 3) = extractBetween(index(:, 2), "x", "_");
    end

    index(~ismember(index(:, 3), basket.idx{i}), :) = [];
    chunkidx = unique(index(:, 1));
    progbar = progressGen(numel(chunkidx)); 

    basket.chunk{i} = cell(numel(chunkidx), 1);
    for k1 = 1:numel(chunkidx) % loop over df chunks in basket i
        progressGen(progbar, k1)

        chunkvars = index(ismember(index(:, 1), chunkidx(k1)), 2);

        % keep only input array(s)/instance(s)
        [inst_idx, arr_idx] = deal(true(numel(chunkvars), 1));
        if opts.ukbrap
            %@21SEP2024: array/instance is different for UKB-RAP dfs
            [~, ia_idx] = ismember(chunkvars, df_in.name); ia_idx(ia_idx < 1) =[];
            df_in_chunk = df_in(ia_idx, :);
            
            if ~isempty(opts.instance) && ~opts.surv
                inst_idx = ismember(df_in_chunk.instance, opts.instance);
                inst_idx(ismissing(df_in_chunk.instance)) = true; % keep df if there is no instance defined
            end

            if ~isempty(opts.array) && ~opts.surv
                arr_idx = ismember(df_in_chunk.array, opts.array);
                arr_idx(ismissing(df_in_chunk.array)) = true; % keep df if there is no array defined
            end

        else
            check_inst_arr = chunkvars.split("_");
            if isvector(check_inst_arr)
                check_inst_arr = check_inst_arr(2:3);
                if iscolumn(check_inst_arr); check_inst_arr = check_inst_arr'; end
            else
                check_inst_arr = check_inst_arr(:, 2:3);
            end
    
            if ~isempty(opts.instance) && ~opts.surv
                inst_idx = ismember(double(check_inst_arr(:, 1)), opts.instance);
            end
    
            if ~isempty(opts.array) && ~opts.surv
                arr_idx = ismember(double(check_inst_arr(:, 2)), opts.array);
            end
        end

        chunkvars = chunkvars(inst_idx & arr_idx);

        if isempty(chunkvars); continue; end

        if opts.ukbrap
            chunkvars = union(chunkvars, "eid"); % eid per each chunk
        end

        basket.chunk{i}{k1} = load(fullfile(basket.folder{i}, ...
            "UKB_" + chunkidx(k1) + ".mat"), chunkvars{:});
        
        if opts.ukbrap
            chunk_eid = basket.chunk{i}{k1}.eid;
            basket.chunk{i}{k1} = rmfield(basket.chunk{i}{k1}, "eid");
        end
        
        cfis = string(fieldnames(basket.chunk{i}{k1}));
        for j = 1:numel(cfis)
            if isstring(basket.chunk{i}{k1}.(cfis(j)))
                basket.chunk{i}{k1}.(cfis(j)) = standardizeMissing(basket.chunk{i}{k1}.(cfis(j)), "");
            end
        end

        if opts.ukbrap
            basket.chunk{i}{k1}.eid = chunk_eid;
            clear chunk_eid
        else
            basket.chunk{i}{k1} = struct2table(basket.chunk{i}{k1}, 'AsArray', true);
        end

    end

    basket.chunk{i}(cellfun(@isempty, basket.chunk{i})) = [];
    if isempty(basket.chunk{i})
        fprintf('instance or arrays cannot be found in this basket!\n')
        continue
    else

        if opts.ukbrap
            [basket.chunk{i}, opts] = rapChunkMerger(basket.chunk{i}, opts);
        else
            basket.chunk{i} = table2struct(horzcat(basket.chunk{i}{:}));
        end

    end

    % convert df(s) to its(their) respective data type(s)?
    if isempty(opts.categorical)

        if opts.ukbrap
            eid = basket.chunk{i}.eid;
            basket.chunk{i} = rmfield(basket.chunk{i}, "eid");
        end

        finames = fieldnames(basket.chunk{i});

        if opts.ukbrap
            [~, vtidx] = ismember(finames, df.name); vtidx(vtidx < 1) = [];
        else
            fidfs = string(extractBetween(finames, "x", "_"));
            [~, vtidx] = ismember(fidfs, df.df); vtidx(vtidx < 1) = [];
        end

        fifnc = df.fnc(vtidx);
        for k1 = 1:numel(finames)
            if ismissing(fifnc(k1)) || fifnc(k1) == ""
                continue
            end

            missidx = ismissing(basket.chunk{i}.(finames{k1}));

            if fifnc(k1) == "datetime"
                tmpval = feval(fifnc(k1), basket.chunk{i}.(finames{k1})(~missidx), Format="dd/MM/uuuu");
                basket.chunk{i}.(finames{k1}) = NaT(numel(basket.chunk{i}.(finames{k1})), 1, Format="dd/MM/uuuu");

            elseif fifnc(k1) == "double"
                tmpval = feval(fifnc(k1), basket.chunk{i}.(finames{k1})(~missidx));
                basket.chunk{i}.(finames{k1}) = nan(numel(basket.chunk{i}.(finames{k1})), 1);
            else
                error('unkonw conversion function! %s', fifnc(k1))
            end
            basket.chunk{i}.(finames{k1})(~missidx) = tmpval;
        end
        basket.chunk{i}.eid = eid;

    else  % categorical data: only keep entries which have query terms
        basket.chunk{i} = struct2table(basket.chunk{i});
                   
        if opts.ukbrap
            eid = basket.chunk{i}.eid;
            basket.chunk{i}.eid = [];
        end

        if opts.surv
            % main codes
            df_in_inc = df_in(df_in.inc, :);
            % inc_tab = subsetUKBBtab(df_in_inc, basket.chunk{i}, eid, opts.categorical);
            inc_tab = subsetUKBBtab(df_in_inc, basket.chunk{i}, eid, opts);
            if ~isempty(inc_tab)
                inc_tab.ex(:) = 0;
                inc_tab.excode(:) = "";
                inc_tab.censor(:) = 1;
            end

            % exclusion codes
            df_in_ex = df_in(~df_in.inc, :);
            if ~isempty(df_in_ex)
                ex_tab = subsetUKBBtab(df_in_ex, basket.chunk{i}, eid, opts);
                if ~isempty(ex_tab)
                    ex_tab.ex(:) = 2;
                    ex_tab.excode(:) = ex_tab.rawUKB;
                    ex_tab.censor(:) = 3;
    
                    % merge with inc_tab and add censoring based on competing
                    % risk (exclusion criteria): same as in getHESIN_DEATH
                    % function
                    
                    if ~isempty(inc_tab)
                        [f1, f2] = ismember(inc_tab.eid, ex_tab.eid); f2(f2<1) = [];
                        if any(f1) % overlapping individuals
                            inc_tab.excode(f1) = ex_tab.rawUKB(f2);
                            inc_tab.ex(f1) = 1;
        
                            % for censoring: check if exclusion (or competing) codes
                            % happened before the outcome diagnosis
                            idx_ex_before_outcome = find(f1);
                            f12 = ex_tab.date(f2) < inc_tab.date(f1);
                            idx_ex_before_outcome = idx_ex_before_outcome(f12);
                            inc_tab.censor(idx_ex_before_outcome) = 3;
                            inc_tab.date(idx_ex_before_outcome) = ex_tab.date(f2(f12));
        
                            ex_tab(f2, :) = []; % overlapping
                        end
        
                        ex_tab.excode = ex_tab.rawUKB;
                        basket.chunk{i} = [inc_tab; ex_tab];
                    else
                        basket.chunk{i} = ex_tab;
                    end
                else
                    basket.chunk{i} = inc_tab;
                end
            else
                basket.chunk{i} = inc_tab;

            end
            
            % long format datasets (e.g. hesin_oper) are no longer needed
            if isfield(opts, "schunk"), opts = rmfield(opts, "schunk"); end
            
            if ~isempty(basket.chunk{i})
                basket.chunk{i}.died(:) = false; % on the assumption that HESIN and DEATH data are being used
                basket.chunk{i}.excode(basket.chunk{i}.excode == "") = missing;
    
                % add registration date
                [f1, f2] = ismember(basket.chunk{i}.eid, opts.regdate.eid);
                % assert(all(f1), "Some individuals don't have registration date!!")

                if ~all(f1)
                    fprintf("\n\tWARNING:%d individuals don't have registration date!\n", sum(~f1))
                    fprintf("\tthis may happen when regdate and chunk data are from different baskets.\n")
                    fprintf("\tthose individuals were dropped.\n")
                    basket.chunk{i}(~f1, :) = []; f2(f2 < 1) = [];
                end
                
                basket.chunk{i}.regdate = opts.regdate.rawUKB(f2);

                % add date of birth
                [f1, f2] = ismember(basket.chunk{i}.eid, opts.dob.eid);
                assert(all(f1), "Some individuals don't have date of birth!!")
                basket.chunk{i}.dob = opts.dob.dob(f2);
    
                % last day of follow-up: unlike HESIN and DEATH data tables
                % (see getHESIN_DEATH function), the source of data for these
                % individuals is unclear. Therefore, we use England HESIN
                % censoring date, as the last day of follow-up.
                cen_path = fullfile(fileparts(which('phenoParser.m')), 'UKBFileParser_HESIN_DEATH');
                cen_date = load(fullfile(cen_path, 'censoring_info.mat')).cd.h.HES;
                idx0 = cen_date < basket.chunk{i}.date;
                basket.chunk{i}.censor(idx0) = 0; % censored because diagnoses after the last day of censoring date
                basket.chunk{i}.date(idx0) = cen_date; % replace censoring date for time-to-even calculation
    
                % censor based on date lost follow-up (overrules last day of
                % follow up censoring date)
                if isfield(opts, "dlf")
                    [f1, f2] = ismember(basket.chunk{i}.eid, opts.dlf.eid); f2(f2<1) = [];
                    if any(f1)
                        idx_lf_before_outcome = find(f1); % lost follow-up before outcome
                        idx = opts.dlf.rawUKB(f2) < basket.chunk{i}.date(f1);
                        idx2 = idx_lf_before_outcome(idx);
                        basket.chunk{i}.censor(idx2) = 0;
                        basket.chunk{i}.date(idx2) = opts.dlf.rawUKB(f2(idx)); % used the date lost-follow up for these individuals
                    end
                end
    
                % calculate time to event
                basket.chunk{i}.tt = years(basket.chunk{i}.date - basket.chunk{i}.regdate);

                % calculate age-of-onset
                basket.chunk{i}.tt0 = years(basket.chunk{i}.date - basket.chunk{i}.dob);
            end

        else
            
            % match input strings to each basket dataset
            [basket.chunk{i}, eid, matchidx] = filterChunkData(basket.chunk{i}, eid, opts);

            if opts.ukbrap
                basket.chunk{i} = rapPruneDelimiter(basket.chunk{i}, opts.categorical, opts.threads);
            end

            basket.chunk{i}.eid = eid(matchidx);
            basket.chunk{i} = table2struct(basket.chunk{i}, 'ToScalar', true);
        end
    end

    if ~isempty(basket.chunk{i}) && isstring(basket.chunk{i}.eid)
       basket.chunk{i}.eid = double(basket.chunk{i}.eid); 
    end

    fprintf('\n')
end

% merge survival data-fields
if opts.surv
    basket = basket.chunk;
    basket(cellfun(@isempty, basket)) = [];
    basket = vertcat(basket{:});
    
    % first unify the death status of individuals across different
    % data-fields, or data baskets
    if any(basket.died)
        tmp = groupsummary(basket, "eid", "sum", "died");
        tmp(tmp.sum_died < 1 | tmp.GroupCount < 2, :) = [];
        basket.died(ismember(basket.eid, tmp.eid)) = true;
    end


    dup_eid = groupsummary(basket, "eid");

    % merge excodes (competing risk) for those duplicated individuals with
    % a different diagnosis date across baskets. Note that, those with the
    % same diagnosis date are being handled later (see below). This is done
    % because one individual may have a diagnosis without a competing risk
    % in one basket, but no diagnosis but competing risk in another basket.
    % Keeping the soonest date results in removing that competing risk,
    % example, where K7 is risk of interest, and CR is competing risk:
    % EID | code |   excode   | tt
    %  1  |  K7  |  <missing> | 2.1
    %  1  |  CR  |     CR     | 4.5
    % So, to keep CR and merged it with K7 in a single row after
    % groupfilter, merge them beforehand (only for those with a different
    % tt)
    % EID | code | excode | tt
    %  1  |  K7  |   CR   | 2.1
    % Note that above doesn't add new information to the regression
    % analysese (cause-specific or subdistribution) since tt is
    % time-to-even for main code K7 which happened before CR. It's only
    % helpful for descriptive overview of individuals that will develop
    % also CR of interests.
    dup_eid_tmp = dup_eid.eid(dup_eid.GroupCount > 1);
    if ~isempty(dup_eid_tmp)
        idx = ismember(basket.eid, dup_eid_tmp);
        tmp = basket(idx, :);
        [G, ueid] = findgroups(tmp.eid);
        excodes_joined = splitapply(@innerSplitSurv, tmp.excode, G);
    
        basket = groupfilter(basket, "eid", @(x) x == min(x), "date"); 
    
        % replace excodes with excodes_joined for intersection of ueid and
        % unique basket.eid
        dup_eid_tmp = groupsummary(basket, "eid");
        dup_eid_tmp = dup_eid_tmp.eid(dup_eid_tmp.GroupCount < 2);
        idx = ismember(ueid, dup_eid_tmp);
        ueid(~idx) = []; excodes_joined(~idx) = [];
        excodes_joined(excodes_joined == "") = missing;
        [idx1, idx2] = ismember(basket.eid, ueid); idx2 = idx2(idx1);
        basket.excode(idx1) = excodes_joined(idx2);
    end

    dup_eid = groupsummary(basket, "eid");
    % find duplicated individuals that have the same date of diagnoses
    % across different data-fields.
    if ~isempty(dup_eid)
        dup_eid = dup_eid.eid(dup_eid.GroupCount > 1);
        idx = ismember(basket.eid, dup_eid);
        tmp = basket(idx, :);
        basket(idx, :) = [];
        
        tmp = fillmissing(tmp, "constant", "", "DataVariables", ["rawUKB", "excode"]);

        %@16MAY2024: deemed unnecessary and problematic if data has more
        %than 2 baskets (i.e. part of eid can end up in t_same and some in t_diff).
        % % find different patterns based on exclusion/censoring
        % tmp2 = groupsummary(tmp, ["eid", "ex", "censor"]);
        
        % idx = tmp2.GroupCount > 1;
        % idx = ismember(tmp.eid, tmp2.eid(idx));
        % t_same = tmp(idx, :); % same ex/censor: pick one and merge codes/df
        % t_diff = tmp(~idx, :); % different ex/censor: should be decided separately

        % if ~isempty(t_same)
        %     t_same_codes = groupsummary(t_same, "eid",...
        %         @(x)join(unique(x, "stable"), ","),...
        %         ["rawUKB", "excode", "df"]);
        %     t_same_codes.GroupCount = [];
        %     t_same_codes.Properties.VariableNames = regexprep(t_same_codes.Properties.VariableNames, "^fun1_", "");
        %     [~, idx] = ismember(t_same_codes.eid, t_same.eid);
        %     dcols = setdiff(colnames(t_same), colnames(t_same_codes));
        %     t_same = [t_same_codes, t_same(idx, dcols)];
        %     t_same = t_same(:, colnames(tmp));
        % end
        
        t_diff = tmp;

        if ~isempty(t_diff)

            if opts.ukbrap
                % only for RAP this may happen
                zero_idx = t_diff.censor == 0;
                if any(zero_idx)
                    t_diff_0 = t_diff(zero_idx, :);

                    % keep only data from HESIN
                    check_hesin_0 = groupfilter(t_diff_0, "eid", @(x) x == "HESIN", "df");
                    assert(all(ismember(t_diff_0.eid, check_hesin_0.eid)), "inconsistency between HESIN and RAP basket! Is HESIN outdated?")
                    basket = [basket; check_hesin_0];
                    t_diff(zero_idx, :) = [];
                end

            end

            assert(~any(t_diff.censor == 0), "Some baskets seem to have different cutpoints for last follow-up!")
            ueid = unique(t_diff.eid);
            t_diff.ex(t_diff.ex < 1) = 3; % temporary (see below)
            t_diff_u = cell(numel(ueid), 1);
            for k = 1:numel(ueid)
                tt = t_diff(t_diff.eid == ueid(k), :);

                % censoring values can be either 1, 2 or 3. Different scenarios can happen:
                %   1 and 3: 1 is case happened before ex code, 3 is excode
                %       happened before: 1 is preferred ONLY is excodes are the
                %       same (e.g. main code comes from HESIN)
                %   If 2 and 3: 3 is chosen, because died of 3 and absent from death registry.
                %   otherwise if at least one 2, then it's 2: other causes
                %   of death


                if all(ismember([2, 3], tt.censor))
                    tt.censor(:) = 3;
                elseif all(ismember([1, 3], tt.censor))
                    tt1 = tt.excode(tt.censor == 1);
                    tt3 = tt.excode(tt.censor == 3);
                    tt1(ismissing(tt1)) = []; tt3(ismissing(tt3)) = [];
                    tt1 = tt1.split(","); tt3 = tt3.split(",");

                    if any(~ismember(tt3, tt1)) % is any excode new that was missing in censoring 1? if so, set censoring to 3
                        tt.censor(:) = 3;
                    else
                        tt.censor(:) = 1;
                    end
                elseif any(tt.censor == 2) && all(tt.excode == "" | ismissing(tt.excode))
                    tt.censor(:) = 2;
                else
                    assert(isscalar(unique(tt.censor)), "Censoring cannot be reconciled for one individual")
                end
                
                % ex code can be 1 (main/ex), 2 (only ex) or 3 (main code).
                % if there is at least one 1 --> 1 (both)
                % if there is no 1, but 2 and 3 --> 1 (both)
                % otherwise --> unique 
                if any(tt.ex == 1) 
                    tt.rawUKB(:) = join(unique(split(tt.rawUKB(tt.ex == 1), ",")), ",");
                    tt.ex(:) = 1;
                elseif all(ismember([2, 3], tt.ex))
                    tt.rawUKB(:) = join(unique(split(tt.rawUKB(tt.ex == 3), ",")), ",");
                    tt.ex(:) = 1;
                else
                    assert(isscalar(unique(tt.ex)), "exclusion criteria cannot be reconciled for one individual")
                end
                
                % join df, excode, rawUKB
                cols = ["rawUKB", "excode", "df"];
                for m = 1:numel(cols)
                    ttcodes = arrayfun(@(x)split(x, ","), tt.(cols(m)), uni=false);
                    ttcodes = unique(vertcat(ttcodes{:}));
                    ttcodes(ttcodes == "" | ismissing(ttcodes)) = [];
                    if isempty(ttcodes)
                        tt.(cols(m))(:) = "";
                    else
                        tt.(cols(m))(:) = ttcodes.join(",");
                    end

                end
                % gtt = groupsummary(tt, "eid", @(x)join(unique(x.split(",")), ","), ["rawUKB", "excode", "df"]);
                % gtt.fun1_rawUKB = regexprep(gtt.fun1_rawUKB, ["\,{1,}", "^\,{1,}", "\,{1,}$"] , [",", "", ""]);
                % gtt.fun1_excode = regexprep(gtt.fun1_excode, ["\,{1,}", "^\,{1,}", "\,{1,}$"] , [",", "", ""]);
                tt = tt(1, :);
                % tt(:, ["rawUKB", "excode", "df"]) = gtt(:, "fun1_" + ["rawUKB", "excode", "df"]);

                t_diff_u{k} = tt;

            end

            t_diff = vertcat(t_diff_u{:});
            t_diff.ex(t_diff.ex == 3) = 0;
        end

        % basket = [basket; t_diff; t_same];
        basket = [basket; t_diff];

    end
    
    basket(:, ["date", "regdate", "dob"]) = [];
    assert(height(basket) == numel(unique(basket.eid)), "Basket should have unique IDs!")

    % add time to events for healthy controls -----------------------------
    % last day of follow-up: unlike HESIN and DEATH data tables (see
    % getHESIN_DEATH function), the source of data for these individuals is
    % unclear. Therefore, we use last available censoring date as the last
    % day of follow-up.
    cen_path = fullfile(fileparts(which('phenoParser.m')), 'UKBFileParser_HESIN_DEATH');
    cen_date = load(fullfile(cen_path, 'censoring_info.mat')).cd;
    cen_date = max(struct2array(struct2array(cen_date)));
    idx = ismember(opts.regdate.eid, basket.eid);
    hc = struct;
    hc.eid = opts.regdate.eid(~idx);
    hc.rawUKB = strings(numel(hc.eid), 1);
    hc = struct2table(hc);
    hc.rawUKB(:) = missing;
    hc.date(:) = cen_date;
    hc.regdate = opts.regdate.rawUKB(~idx);

    % add birth date
    [~, idx] = ismember(hc.eid, opts.dob.eid);
    hc.dob = opts.dob.dob(idx);

    hc.ex(:) = nan;
    hc.excode = strings(numel(hc.eid), 1);
    hc.excode(:) = missing;
    hc.censor(:) = 0;
    hc.died(:) = false;

    % censor based on date lost follow-up (overrules last day of
    % follow up censoring date)
    if isfield(opts, "dlf")
        [f1, f2] = ismember(hc.eid, opts.dlf.eid); f2(f2<1) = [];
        if any(f1)
            idx_lf_before_outcome = find(f1); % lost follow-up before outcome
            idx = opts.dlf.rawUKB(f2) < hc.date(f1);
            idx2 = idx_lf_before_outcome(idx);
            hc.date(idx2) = opts.dlf.rawUKB(f2(idx)); % used the date lost-follow up for these individuals
        end
    end

    % calculate time to event
    hc.tt = years(hc.date - hc.regdate);
    
    % calculate age-of-onset
    hc.tt0 = years(hc.date - hc.dob);
    hc.df(:) = "Control";
    hc(:, ["date", "regdate", "dob"]) = [];
    basket = [basket; hc];
    basket.Properties.UserData = cen_date; % last day of follow-up
    return

end

% if there are baskets created on the same date, keep the one that has
% largest sample size (eid)
dupdates = datetime(duplicates(basket.date));

if ~isempty(dupdates)
    idx = ismember(basket.date, dupdates);
    rembasket = basket(idx, :);
    basket(idx, :) = [];
    
    dupbasket = cell(numel(dupdates), 1);
    for i = 1:numel(dupdates)
        dupdateidx = ismember(rembasket.date, dupdates(i));
        dupbasket{i} = rembasket(dupdateidx, :);
        rembasket(dupdateidx, :) = [];

        for k1 = 1:height(dupbasket{i})-1
            for k2 = k1+1:height(dupbasket{i})
                [idx1, idx2] = ismember(dupbasket{i}.chunk{k1}.eid, dupbasket{i}.chunk{k2}.eid);
                fi1 = setdiff(fieldnames(dupbasket{i}.chunk{k1}), 'eid');
                fi2 = setdiff(fieldnames(dupbasket{i}.chunk{k2}), 'eid');
                fi12 = intersect(fi1, fi2);
                for j = 1:numel(fi12)
                    df1 = dupbasket{i}.chunk{k1}.(fi12{j})(idx1);
                    df2 = dupbasket{i}.chunk{k2}.(fi12{j})(idx2(idx1));
                    if isnumeric(df1)
                        df1(isnan(df1)) = [];  df2(isnan(df2)) = [];
                    else
                        df1(ismissing(df1)) = []; df2(ismissing(df2)) = [];
                    end

                    if numel(df1) > numel(df2) % keep df in chunk k1
                        dupbasket{i}.chunk{k2} = rmfield(dupbasket{i}.chunk{k2}, fi12{j});
                        dupbasket{i}.idx{k2} = setdiff(dupbasket{i}.idx{k2}, extractBetween(fi12{j}, "x", "_"));
                    elseif numel(df1) < numel(df2)
                        dupbasket{i}.chunk{k1} = rmfield(dupbasket{i}.chunk{k1}, fi12{j});
                        dupbasket{i}.idx{k1} = setdiff(dupbasket{i}.idx{k1}, extractBetween(fi12{j}, "x", "_"));

                        % check eid size if df1 and df2 are of the same
                        % size
                    elseif numel(dupbasket{i}.chunk{k1}.eid) >= numel(dupbasket{i}.chunk{k2}.eid)
                        dupbasket{i}.chunk{k2} = rmfield(dupbasket{i}.chunk{k2}, fi12{j});
                        dupbasket{i}.idx{k2} = setdiff(dupbasket{i}.idx{k2}, extractBetween(fi12{j}, "x", "_"));
                    else
                        dupbasket{i}.chunk{k1} = rmfield(dupbasket{i}.chunk{k1}, fi12{j});
                        dupbasket{i}.idx{k1} = setdiff(dupbasket{i}.idx{k1}, extractBetween(fi12{j}, "x", "_"));
                    end
                end
            end
        end
        
        dupbasket{i}(cellfun(@isempty, dupbasket{i}.idx), :) = [];
        
        % remove baskets with no fields (may not be the same as idx)
        dup_emp_rm_idx = cellfun(@isempty, ...
            cellfun(@(x)setdiff(fieldnames(x), "eid"), ...
            dupbasket{i}.chunk, 'uni', false));
        dupbasket{i}(dup_emp_rm_idx, :) = [];

        if ~isempty(opts.categorical) % remove missing entries
            for j = 1:height(dupbasket{i})
                dupbasket{i}.chunk{j} = struct2table(dupbasket{i}.chunk{j});
                tmpeid = dupbasket{i}.chunk{j}.eid;
                dupbasket{i}.chunk{j}.eid = [];
                [dupbasket{i}.chunk{j}, rmidx] = rmmissing(dupbasket{i}.chunk{j}, 1, 'MinNumMissing', width(dupbasket{i}.chunk{j}));
                dupbasket{i}.chunk{j} = rmmissing(dupbasket{i}.chunk{j}, 2, 'MinNumMissing', height(dupbasket{i}.chunk{j}));
                dupbasket{i}.chunk{j} = table2struct(dupbasket{i}.chunk{j}, 'ToScalar', true);
                dupbasket{i}.chunk{j}.eid = tmpeid(~rmidx);
            end
        end
    end

    basket = [basket; vertcat(dupbasket{:})];
end

basket = sortrows(basket, 'date', 'descend');

% if HESIN_DEATH basket is present and 'categorical' field is empty, this
% means the user wants all the info in the queried categorical data-fields,
% and query is empty. To avoid errors when writing to a phenotype struct,
% the chunk struct from HESIN_DEATH basket should be modified and
% consistent with other baskets (a struct with eid and df fields)
hesinidx = find(endsWith(lower(basket.name), "hesin_death"));
if isempty(opts.categorical) && ~isempty(hesinidx)
    for i = 1:numel(hesinidx) % can multiple HESIN_DEATH baskets be present?!
        thisdf = string(basket.idx{hesinidx(i)});
        chunk = basket.chunk{hesinidx(i)};
        fis = setdiff(fieldnames(chunk), ["eid", "rawUKB"]);
        chunk = rmfield(chunk, fis);
        for j = 1:numel(thisdf), chunk.("x" + thisdf(j) + "_0_0") = chunk.rawUKB; end
        basket.chunk{hesinidx(i)} = rmfield(chunk, "rawUKB");
    end
end

if ~isempty(opts.categorical)
    % merge all baskets: categorical queries need to merged into a single
    % variable. For instance, an ICD-10 code among different data-fields
    % (cancer, death or hesin) an over all existing instances/arrays needs
    % a single phenotype for eids having this ICD-10 code.
    
    % HESIN_DEATH should be handled differently
    hd_basket_idx = endsWith(basket.name, '_HESIN_DEATH');
    hd_basket = basket(hd_basket_idx, :);
    basket(hd_basket_idx, :) = [];

    % @31MAR2023: remove empty baskets (none of query codes were found in this basket)
    miss_idx = cellfun(@(x)isempty(x.eid), basket.chunk);
    basket(miss_idx, :) = [];

    if isempty(basket) 
        % only hesin is available
        if ~isempty(hd_basket)
            basket = hd_basket.chunk{1};
            basket.chunk = basket.rawUKB;
            basket.df = hd_basket.idx{1};
            basket.codes = unique(basket.chunk);
            if any(contains(basket.codes, ","))
                basket.codes = arrayfun(@(x)split(unique(x), ","), basket.codes, uni=false);
                basket.codes = unique(vertcat(basket.codes{:}));
            end
            basket.basket = string(hd_basket.name);
            if isfield(basket, 'ex')
                basket.ex.chunk = basket.ex.rawUKB;    
            end
        end
        return
    end
    
    eid = cellfun(@(x) x.eid, basket.chunk, 'uni', false);
    eid = unique(vertcat(eid{:}));

    mergedb.chunk = cell(height(basket), 1);
    for i = 1:height(basket)
        basket.chunk{i} = struct2table(basket.chunk{i});
        [idx1, idx2] = ismember(eid, basket.chunk{i}.eid);
        idx2(idx2 < 1) = [];
        
        cols = ~ismember(basket.chunk{i}.Properties.VariableNames, "eid");
        mergedb.chunk{i} = strings(numel(eid), sum(cols));
        mergedb.chunk{i}(idx1, :) = basket.chunk{i}{idx2, cols};
        mergedb.chunk{i} = array2table(mergedb.chunk{i}, 'VariableNames', ...
            basket.chunk{i}.Properties.VariableNames(cols));
    end
    
    mergedb.chunk = horzcat(mergedb.chunk{:});
    
    % merge dfs
    if opts.ukbrap
        tmp_df = ismember(df.name, colnames(mergedb.chunk));
        mergedb.df = unique(df.df(tmp_df)).';
    else
        mergedb.df = extractBetween(string(mergedb.chunk.Properties.VariableNames), "x", "_");
    end
    mergedb.source = repmat(mergedb.df, numel(eid), 1);

    if ~opts.ukbrap
        mergedb.source(ismissing(mergedb.chunk)) = "";
    end
    mergedb.source = join(mergedb.source, ",", 2);
    mergedb.source = regexprep(mergedb.source, '(\,+)', ',');
    mergedb.source = regexprep(mergedb.source, '(^,$|,$|^,)', '');
    % idx = contains(mergedb.source, ",");
    % tmpvec = arrayfun(@(x)join(unique(split(x, ",")), ","), mergedb.source(idx));
    % mergedb.source(idx) = tmpvec;
    mergedb.df = unique(mergedb.df);
    
    % merge all codings (comma separated)
    mergedb.chunk = fillmissing(mergedb.chunk, 'constant', "");
    mergedb.chunk = join(mergedb.chunk.Variables, ",", 2);
    mergedb.chunk = regexprep(mergedb.chunk, '(\,+)', ',');
    mergedb.chunk = regexprep(mergedb.chunk, '(^,$|,$|^,)', '');
    if isempty(hd_basket) % otherwise, will do at once when merging with HESIN_DEATH basket
        mergedb.codes = arrayfun(@(x)split(x, ","), unique(mergedb.chunk), 'uni', false);
        mergedb.codes = unique(vertcat(mergedb.codes{:}));
        % idx = contains(mergedb.chunk, ",");
        % tmpvec = arrayfun(@(x)join(unique(split(x, ",")), ","), mergedb.chunk(idx));
        % mergedb.chunk(idx) = tmpvec;
    end
    mergedb.eid = eid;
    mergedb.basket = string(basket.name);   
    
    if ~isempty(hd_basket) % merge with HESIN_DEATH data
        eid = union(mergedb.eid, hd_basket.chunk{1}.eid);
        oldsource = mergedb.source;
        [chunk, mergedb.source] = deal(strings(numel(eid), 2));
    
        [idx1, idx2] = ismember(eid, mergedb.eid); idx2(idx2 < 1) = [];
        chunk(idx1, 1) = mergedb.chunk(idx2);
        mergedb.source(idx1, 1) = oldsource;
    
        [idx1, idx2] = ismember(eid, hd_basket.chunk{1}.eid); idx2(idx2 < 1) = [];
        chunk(idx1, 2) = hd_basket.chunk{1}.rawUKB(idx2);
        mergedb.source(idx1, 2) = hd_basket.chunk{1}.source(idx2);
        mergedb.chunk = chunk;
        
        joinfi = ["chunk", "source"];
        for i = 1:numel(joinfi)
            mergedb.(joinfi(i)) = mergedb.(joinfi(i)).join(","); 
            mergedb.(joinfi(i)) = regexprep(mergedb.(joinfi(i)), '(^,$|,$|^,)', '');
            if joinfi(i) == "chunk"
                idx = contains(mergedb.chunk, ",");
                tmpvec = arrayfun(@(x)join(unique(split(x, ",")), ","), mergedb.chunk(idx));
                mergedb.chunk(idx) = tmpvec;
            end
        end
        mergedb.codes = arrayfun(@(x)split(x, ","), unique(mergedb.chunk), 'uni', false);
        mergedb.codes = unique(vertcat(mergedb.codes{:}));
    
        mergedb.eid = eid;
        mergedb.basket = [mergedb.basket; string(hd_basket.name)];
        mergedb.df = union(mergedb.df, hd_basket.idx{1});
    
    end
    
    % apply exclusion criteria
    if all(~ismissing(opts.remove2{opts.ct}))
        if ~isempty(hd_basket) && isfield(hd_basket.chunk{1}, 'ex')
            rmidx2 = ismember(mergedb.eid, hd_basket.chunk{1}.ex.eid);
            mergedb.eid(rmidx2) = [];
            mergedb.chunk(rmidx2) = [];
            mergedb.source(rmidx2) = [];
        end

        rmidx = contains(mergedb.chunk, regexpPattern("(,|^)" + opts.remove2{opts.ct} + "(,|$)"));
        mergedb.ex.eid = mergedb.eid(rmidx);
        mergedb.ex.chunk = mergedb.chunk(rmidx);
        mergedb.ex.source = mergedb.source(rmidx);
    
        mergedb.eid(rmidx) = [];
        mergedb.chunk(rmidx) = [];
        mergedb.source(rmidx) = [];
        
        if ~isempty(hd_basket) && isfield(hd_basket.chunk{1}, 'ex')
            mergedb.ex.eid = [mergedb.ex.eid; hd_basket.chunk{1}.ex.eid];
            mergedb.ex.chunk = [mergedb.ex.chunk; hd_basket.chunk{1}.ex.rawUKB];
            mergedb.ex.source = [mergedb.ex.source; hd_basket.chunk{1}.ex.source];
            if numel(mergedb.ex.eid) ~= numel(unique(mergedb.ex.eid))
                error('something went wrong with exclusion criteria')
            end
        end

        mergedb.codes = arrayfun(@(x)split(x, ","), ...
            union(unique(mergedb.chunk), unique(mergedb.ex.chunk)), 'uni', false);
        mergedb.codes = unique(vertcat(mergedb.codes{:}));
    end
    
    basket = mergedb; clear mergedb
end

end % END

%% ------------------------------------------------------------------------
function x = innerSplitSurv(x)
x = arrayfun(@(y)split(y, ","), x, uni=false);
x = vertcat(x{:});
x(ismissing(x) | x == "") = [];
if isempty(x)
    x = "";
else
    x = join(unique(x), ",");
end

end

%% ------------------------------------------------------------------------
function out = subsetUKBBtab(df, tab, eid, opts, codes)
% subsets the full table "tab" and creates a table of individuals with
% codes in either "codes" or in "df" struct. Only used when 'surv' is true.
% 
if nargin < 5
    codes = []; % used only for main inclusion criteria
end

% keep only dfs in df tab
cols = colnames(tab);
if opts.ukbrap
    df.df = df.name;

    ins_idx = ~ismissing(df.instance);
    if any(ins_idx) % keep only queried instances
        df_tmp1 = df(ins_idx, :);
        df_tmp1(df_tmp1.instance ~= df_tmp1.ii, :) = [];
        df_tmp2 = df(~ins_idx, :);
        tab = tab(:, intersect(cols, union(df_tmp1.df, df_tmp2.df)));

    else
        tab = tab(:, intersect(cols, df.df));
    end
else
    for k = 1:height(df)
        df.patt(k) = "x" + df.df(k) + "_";
        if ~ismissing(df.instance(k))
            df.patt(k) = df.patt(k) + df.instance(k) + "_";
        end
    
        % also add codes
        if ~isempty(codes) && df.relation(k) == ""
            df.codes(k) = join(codes, ",");
        end 
    end
    cols(~startsWith(cols, df.patt)) = [];
    tab = tab(:, cols);
end

tab.eid = eid;
cols = colnames(tab);

cflag = true; ct = 1;
out = ({});
while cflag
    % detect date-field and it's related data-field
    if isempty(df), cflag = false; end
    if ~cflag, break; end
    didx = find(df.relation ~= "", 1);

    if ~isempty(didx) % is there any related date-field?

        % for some fields in UKB-RAP number of date dfs are more than data
        % dfs (e.g. 41202 for main ICD10 codes). data dfs are '|' joined
        % into a single df by DNANexus. For these fields we skip matching
        % and rely on HESIN data (which cover these fields).
        dup_date_idx = df.relation == df.relation(didx);
        if opts.ukbrap && sum(dup_date_idx) > 1
            df1 = df(dup_date_idx, :);
            df2 = df(ismember(df.df, df1.relation), :);
            df(ismember(df.df, union(df1.df, df2.df)), :) = [];
            fprintf("\t%s field has multiple date dfs, skipped.\n", df2.df)
            continue
        end

        df1 = df(didx, :); % date_field
        df2 = df(df.df == df.relation(didx), :); % data_field
        df(ismember(df.df, union(df1.df, df2.df)), :) = [];
        
        if opts.ukbrap
            keep_date_col_idx = ismember(cols, union(df1.df, "eid"));
            keep_data_col_idx = ismember(cols, union(df2.df, "eid"));

        else
            keep_date_col_idx = startsWith(cols, union(df1.patt, "eid"));
            keep_data_col_idx = startsWith(cols, union(df2.patt, "eid"));
        end

        date_tab = tab(:, keep_date_col_idx);
        data_tab = tab(:, keep_data_col_idx);
        
        long_tab_flag = false;
        if width(data_tab) < 2 % i.e. df is not in 'tab' but in opts.schunk (e.g. hesin_oper datasets)
            if ~isfield(opts, "schunk")

                % @16NOV2024: rapChunkMerger removes missing vectors, so
                % sometimes the reason is that there is no observation.
                % Error was changed to a warning, and skipping
                % error("long format tables are missing!")
                fprintf("\tWarning:subsetUKBBtab: long format tables are missing for %s\n", df1.df)
                fprintf("\tskipping %s and %s\n", df1.df, df2.df)

                continue
            end

            % table which has both date/data columns
            matched_schunk = cellfun(@(x)all(ismember(union(df1.df, df2.df), colnames(x))), opts.schunk);
            assert(sum(matched_schunk) == 1, "multiple matches to 'schunk' tables!")
            date_tab = opts.schunk{matched_schunk}(:, union(df1.df, "eid"));
            data_tab = opts.schunk{matched_schunk}(:, union(df2.df, "eid"));
            long_tab_flag = true;
        end

        assert(width(date_tab) == width(data_tab), "data/date tabe must have the same size!")

        if long_tab_flag
            [data_tab, ~, matched_idx] = matchDataTab(data_tab, df2.codes);
            date_tab(~matched_idx, :) = [];
            assert(height(data_tab) == height(date_tab), "something went wrong with long date-/data-table matching!")
            clear matched_idx

        else
            data_tab = matchDataTab(data_tab, df2.codes);
            if isempty(data_tab)
                % df(ismember(df.df, union(df1.df, df2.df)), :) = [];
                continue; 
            end
            date_tab(~ismember(date_tab.eid, data_tab.eid), :) = [];
        end

        % convert date-time to string for consistency 
        for k = 1:width(date_tab)
            if isdatetime(date_tab.(k))
                date_tab.(k) = string(date_tab.(k));
            end
        end

        dcols = setdiff(colnames(data_tab), "eid");

        if opts.ukbrap
            [~, keep_idx] = ismember(dcols, df1.relation);
            dcols = df1.df(keep_idx);
        else
            dcols = dcols.replace(textBoundary("start") + "x" + digitsPattern + "_", "x" + df1.df + "_");
        end

        date_tab = date_tab(:, union("eid", dcols));

        % unify column names and order
        dcols = colnames(date_tab);
            
        if opts.ukbrap
            dcols(dcols == "eid") = [];
            [~, match_idx] = ismember(dcols, df1.df);
            dcols_new = df1.relation(match_idx);
            dcols_new = union("eid", dcols_new, "Stable");
            dcols = union("eid", dcols, "Stable");
        else
            dcols_new = dcols.replace(textBoundary("start") + df1.patt, df2.patt);
        end

        date_tab = renamevars(date_tab, dcols, dcols_new);
        data_tab = data_tab(:, dcols_new);
        date_tab = date_tab{:, :};
        date_tab(ismissing(data_tab)) = missing;
        date_tab = array2table(date_tab, VariableNames=colnames(data_tab));
        date_tab.eid = double(date_tab.eid);

        if long_tab_flag
            assert(all(date_tab.eid == data_tab.eid), "data/date tab eid mismatch for long table!")
            
            % row indices to keep track of long formatted tables with
            % duplicated eids
            [date_tab.Row, data_tab.Row] = deal(string(1:height(data_tab)));
        end

        % pick the date that the code happened first
        date_tab = movevars(date_tab, "eid", Before=1);
        date_tab = stack(date_tab, 2:width(date_tab));
        date_tab = rmmissing(date_tab);
        date_tab.Properties.VariableNames = ["eid", "df", "date"];
        date_tab.df = string(date_tab.df);
        if df1.fnc == "datetime"
            date_tab.date = datetime(date_tab.date, Format="dd/MM/uuuu");
        else % must be year?
            date_tab.date = feval(df1.fnc, date_tab.date);

            % e.g. df20006 can have negative values -1, or -3, these ids
            % will be removed later (age was not clear). We check if there
            % are subjects with multiple dates and at least one -1 or -3,
            % negative values will be removed here. 
            dups = groupsummary(date_tab, "eid", @min, "date");
            dups = dups.eid(dups.GroupCount > 1 & dups.fun1_date < 0);
            if ~isempty(dups)
                % there should be at least one positive value
                dups2 = groupfilter(date_tab(ismember(date_tab.eid, dups), :), "eid", @(x) x > 0, "date");
                rm_idx = ismember(date_tab.eid, dups2.eid) & date_tab.date < 0;
                date_tab(rm_idx, :) = [];
            end

            date_tab.date = datetime(0,1,0, "Format", "dd/MM/uuuu") + years(date_tab.date); 
        end

        if long_tab_flag % duplicated eids
            % under a reasonable assumption that for long tables (e.g.
            % hesin_oper) there is only one date field, otherwise
            % phenoParser would throw an error earlier (see where date
            % dfs are fetched).
            assert(all(ismember(date_tab.Row, data_tab.Row)), "date_tab in long format dataset has more than date column")

             % subset data_tab if missing values removed from date_tab
            if height(data_tab) ~= height(date_tab)
                [~, idx] = ismember(date_tab.Row, data_tab.Row);
                data_tab = data_tab(idx, :);
            end
        end

        date_tab = groupfilter(date_tab, "eid", @(x) x == min(x), "date");

        % pick the first array_indices if there are individuals with the
        % same date (arbitrary filter)
        date_tab = groupfilter(date_tab, "eid", @natsort2, "df");

        % add codes to date_tab now with the same array_indices
        data_tab = movevars(data_tab, "eid", Before=1);
        data_tab = stack(data_tab, 2:width(data_tab));
        data_tab = rmmissing(data_tab);
        data_tab.Properties.VariableNames = ["eid", "df", "code"];
        data_tab = convertvars(data_tab, ["df", "code"], @string);
        
        % groupfilter on date_tab picks the first occurrence only 
        if long_tab_flag && height(data_tab) ~= height(date_tab) 
            [~, idx] = ismember(date_tab.Row, data_tab.Row);
            data_tab = data_tab(idx, :);
        end
        
        data_tab = join(date_tab, data_tab, "Keys", ["eid", "df"]);

    else % no date_field
        df2 = df(1, :);
        df(1, :) = [];
        
         if opts.ukbrap
            keep_data_col_idx = ismember(cols, union(df2.df, "eid"));
         else
            keep_data_col_idx = startsWith(cols, union(df2.patt, "eid"));
         end

        data_tab = tab(:, keep_data_col_idx);
        data_tab = matchDataTab(data_tab, df2.codes);
        if isempty(data_tab), continue; end

        data_tab = movevars(data_tab, "eid", Before=1);
        
        if opts.ukbrap
            data_tab = rapPruneDelimiter(data_tab, df2.codes.split(","), opts.threads);
        end

        data_tab = stack(data_tab, 2:width(data_tab));
        data_tab = rmmissing(data_tab);
        data_tab.Properties.VariableNames = ["eid", "df", "code"];
        data_tab = convertvars(data_tab, ["df", "code"], @string);
        data_tab = groupfilter(data_tab, "eid", @natsort2, "df");
        data_tab.date = NaT(height(data_tab), 1);
        data_tab = movevars(data_tab, "date", After="df");
    end
    
    out{ct, 1} = data_tab;
    ct = ct + 1;

end

out(cellfun(@isempty, out)) = [];
out = vertcat(out{:});

if ~isempty(out) && any(~ismissing(out.date))

    % for duplicated individuals among different data-fields, pick the
    % soonest date 
    out = groupfilter(out, "eid", @(x) x == min(x), "date");
    
    % if there are still duplicated individuals because of same date of
    % diagnosis across different data-fields, merge their codes
    if ~isempty(duplicates(out.eid))
        out = groupsummary(out, "eid", @(x)join(string(unique(x, "stable")), ","));
        out.GroupCount = [];
        out.Properties.VariableNames = regexprep(out.Properties.VariableNames, "^fun1_", "");
        out.date = datetime(out.date, Format="dd/MM/uuuu");
    end

    out = renamevars(out, "code", "rawUKB");
else
    out = table;
end

end % END

%% ---------------------------------------------------------------------
function input = natsort2(input)
    [~, ndx] = natsort(input);
    input = false(numel(input), 1);
    input(ndx(1)) = true; % keep the first
    
end

%% ---------------------------------------------------------------------
function [tab, matchidx, matchidx_long] = matchDataTab(tab, codings)
eid = tab.eid;
tab.eid = [];
codings = split(codings, ",");

matchidx = ~ismissing(tab);
for k = 1:size(matchidx, 2)

    % 24SEP2024: UKB-RAP aggregated dfs
    if any(tab.(k).contains("|")) 
        patt = (textBoundary("start") | "|") + ...
            codings + ...
            (textBoundary("end") | "|");
        fix = contains(tab.(k)(matchidx(:, k)), patt);
        
    else
        fix = ismember(tab.(k)(matchidx(:, k)), codings);
    end

    matchidx(matchidx(:, k), k) = fix;
end

[matchidx, matchidx_long] = deal(any(matchidx, 2));
tab(~matchidx, :) = [];
eid(~matchidx) = [];
tab = rmmissing(tab, 2, 'MinNumMissing', height(tab));
matchidx = ~ismissing(tab);
for k = 1:size(matchidx, 2)

    % 24SEP2024: UKB-RAP aggregated dfs
    if any(tab.(k).contains("|")) 
        patt = (textBoundary("start") | "|") + ...
            codings + ...
            (textBoundary("end") | "|");
        fix = contains(tab.(k)(matchidx(:, k)), patt);
    else
        fix = ismember(tab.(k)(matchidx(:, k)), codings);
    end

    matchidx(matchidx(:, k), k) = fix;
end
tab = fillmissing(tab, 'constant', "", 'MissingLocations', ~matchidx);
tab = standardizeMissing(tab, "");
matchidx2 = any(matchidx, 2);
tab(~matchidx2, :) = [];
tab = rmmissing(tab, 2, 'MinNumMissing', height(tab));
tab.eid = eid(matchidx2);
matchidx = ~ismissing(tab);

end % END

%% ------------------------------------------------------------------------
function outpheno = writePheno2Mat(phenoindex, opts)

UKB_STRUCT_ALL = struct;
UKB_STRUCT_ALL.info = struct;
phenoindex.basket = regexprep(phenoindex.basket, "^UKBFileParser_", "");

if ~isempty(opts.categorical)
    
    UKB_STRUCT_ALL.eid = phenoindex.eid;
    UKB_STRUCT_ALL.rawUKB = phenoindex.chunk;
    UKB_STRUCT_ALL.source = phenoindex.source;
    UKB_STRUCT_ALL.numericFlag = false;
    UKB_STRUCT_ALL.info.basket = join(phenoindex.basket, ",");
    UKB_STRUCT_ALL.info.date = datetime('now');
    UKB_STRUCT_ALL.info.df = join(phenoindex.df, ",");
    try UKB_STRUCT_ALL.info.desc = phenoindex.desc; catch, end

    if isfield(phenoindex, 'ex')
        UKB_STRUCT_ALL.exeid = phenoindex.ex.eid;
        UKB_STRUCT_ALL.ex.rawUKB = phenoindex.ex.chunk;
        UKB_STRUCT_ALL.ex.source = phenoindex.ex.source;
    end
    
    % write term meanings .................................................
    UKB_STRUCT_ALL.termMeaning = strings(numel(phenoindex.chunk), 1);
    if isfield(phenoindex, 'ex')
        UKB_STRUCT_ALL.ex.termMeaning = strings(numel(phenoindex.ex.chunk), 1);
    end

    for i = 1:height(phenoindex.meaning)
        if ismember(phenoindex.meaning.Value(i), opts.remove2{opts.ct}) && isfield(phenoindex, 'ex')
            codidx = contains(phenoindex.ex.chunk, ...
                regexpPattern("(,|^)" + phenoindex.meaning.Value(i) + "(,|$)"));
            UKB_STRUCT_ALL.ex.termMeaning(codidx) = ...
                UKB_STRUCT_ALL.ex.termMeaning(codidx) + ";" + ...
                phenoindex.meaning.Meaning(i);
            continue
        end

        codidx = contains(phenoindex.chunk, ...
            regexpPattern("(,|^)" + phenoindex.meaning.Value(i) + "(,|$)"));
            UKB_STRUCT_ALL.termMeaning(codidx) = ...
                UKB_STRUCT_ALL.termMeaning(codidx) + ";" + ...
                phenoindex.meaning.Meaning(i);
    end
    UKB_STRUCT_ALL.termMeaning = regexprep(UKB_STRUCT_ALL.termMeaning, '(\;+)', ';');
    UKB_STRUCT_ALL.termMeaning = regexprep(UKB_STRUCT_ALL.termMeaning, '(^;$|;$|^;)', '');
    
    if isfield(phenoindex, 'ex')
        UKB_STRUCT_ALL.ex.termMeaning = regexprep(UKB_STRUCT_ALL.ex.termMeaning, '(\;+)', ';');
        UKB_STRUCT_ALL.ex.termMeaning = regexprep(UKB_STRUCT_ALL.ex.termMeaning, '(^;$|;$|^;)', '');
    end

    % .....................................................................
    if isfield(phenoindex, 'tag') && phenoindex.tag ~= ""
        UKB_STRUCT_ALL.tag = phenoindex.tag;
        phenoindex.out = matlab.lang.makeValidName(phenoindex.tag);
    else
        UKB_STRUCT_ALL.tag = join(phenoindex.codes, ",");
        phenoindex.out = UKB_STRUCT_ALL.tag + ...
            "_df" + join(phenoindex.df, ".");
        phenoindex.out = matlab.lang.makeValidName(phenoindex.out);
    end

     % if already exists, append date time to the file name.
    if isfile(fullfile(opts.phenodir, phenoindex.out + ".mat")) && opts.save
        % fprintf('%s is already present, new name is generated.\n', phenoindex.out)

        oans = input(sprintf("\t%s is already present, overwrite? [y/N]: ", phenoindex.out), "s");
        oans = lower(string(oans));
        
        if oans == "" || oans == "n"
            phenoindex.out = phenoindex.out + ...
                string(string(datetime)).replace(["-"," ", ":"], ".");
            fprintf("\t\tnew name was generated %s\n", phenoindex.out)
        end
    end
    phenoindex.out = phenoindex.out + ".mat";
    
    if opts.save
        if isfile(fullfile(opts.phenodir, phenoindex.out))
            fprintf("\t\toverwriting %s \n", phenoindex.out)
        % else % still file is present
            % fprintf('cannot save file %s to %s\n', phenoindex.out, opts.phenodir)
        end
        save(fullfile(opts.phenodir, phenoindex.out), "UKB_STRUCT_ALL");
    end
    outpheno = UKB_STRUCT_ALL;

    return
end

% for queries not in codings (e.g. numeric phenotypes, or QC flags, or
% categorical traits with a data-field query rather a coding value)

if numel(phenoindex.tag) ~= numel(unique(phenoindex.tag))
    % if multiple arrays/instances are present per each phenotype
    inst_arr = split(phenoindex.fi, "_");
    if iscolumn(inst_arr); inst_arr = inst_arr'; end
    phenoindex.instance = "i" + inst_arr(:, 2);

    try
        phenoindex.array = "a" + inst_arr(:, 3);
        phenoindex.out = matlab.lang.makeValidName(phenoindex.tag) + "." + ...
            phenoindex.instance + "." + phenoindex.array;
    catch % no array
        phenoindex.out = matlab.lang.makeValidName(phenoindex.tag) + "." + ...
            phenoindex.instance;
    end
    
    
else
    phenoindex.out = matlab.lang.makeValidName(phenoindex.tag);
end

outpheno = cell(numel(phenoindex.fi), 1);
for i = 1:numel(phenoindex.fi)
    UKB_STRUCT_ALL.eid = phenoindex.raw.eid;
    UKB_STRUCT_ALL.rawUKB = phenoindex.raw.(phenoindex.fi(i));
    UKB_STRUCT_ALL.numericFlag = phenoindex.numeric(i);
    UKB_STRUCT_ALL.tag = phenoindex.tag(i);

    if isfield(phenoindex, "unit")
        UKB_STRUCT_ALL.unit = phenoindex.unit(i);
    end

    UKB_STRUCT_ALL.info.basket = phenoindex.basket;
    UKB_STRUCT_ALL.info.date = datetime('now');
    UKB_STRUCT_ALL.info.df = phenoindex.df(i);
    UKB_STRUCT_ALL.info.dfraw = phenoindex.fi(i);
    try UKB_STRUCT_ALL.info.desc = phenoindex.desc(i); catch, end

    if ~phenoindex.numeric(i) && isnumeric(UKB_STRUCT_ALL.rawUKB)
        UKB_STRUCT_ALL.rawUKB = string(UKB_STRUCT_ALL.rawUKB);
    end

    if phenoindex.numeric(i)
        if isdatetime(UKB_STRUCT_ALL.rawUKB)
            rmidx = isnat(UKB_STRUCT_ALL.rawUKB);
        else
            rmidx = isnan(UKB_STRUCT_ALL.rawUKB);
        end
    else
        rmidx = ismissing(UKB_STRUCT_ALL.rawUKB);
    end

    UKB_STRUCT_ALL.eid(rmidx) = [];
    UKB_STRUCT_ALL.rawUKB(rmidx) = [];

    if isempty(phenoindex.meaning{i})
        UKB_STRUCT_ALL.termMeaning = '';
    else
        try
            if isnumeric(UKB_STRUCT_ALL.rawUKB)
                [~, idx] = ismember(UKB_STRUCT_ALL.rawUKB, double(phenoindex.meaning{i}.Value));
            else
                [~, idx] = ismember(UKB_STRUCT_ALL.rawUKB, phenoindex.meaning{i}.Value);
            end
            UKB_STRUCT_ALL.termMeaning = phenoindex.meaning{i}.Meaning(idx);
        catch % e.g. df 36 which doesn't have a coding table
            UKB_STRUCT_ALL.termMeaning = '';
        end
    end

    % if already exists, append date time to the file name.
    if isfile(fullfile(opts.phenodir, phenoindex.out(i) + ".mat")) && opts.save

        oans = input(sprintf("\t%s is already present, overwrite? [y/N]: ", phenoindex.out(i)), "s");
        oans = lower(string(oans));
        
        if oans == "" || oans == "n"
            phenoindex.out(i) = phenoindex.out(i) + ...
                string(string(datetime)).replace(["-"," ", ":"], ".");
            fprintf("\t\tnew name was generated %s\n", phenoindex.out(i))
        end
    end
    
    phenoindex.out(i) = phenoindex.out(i) + ".mat";
    if opts.save
        if isfile(fullfile(opts.phenodir, phenoindex.out(i)))
            fprintf("\t\toverwriting %s \n", phenoindex.out(i))   
        % else % still file is present
        %     fprintf('cannot save file %s to %s\n', phenoindex.out(i), opts.phenodir)
        end

        save(fullfile(opts.phenodir, phenoindex.out(i)), "UKB_STRUCT_ALL");
    end
    outpheno{i, 1} = UKB_STRUCT_ALL;
end

try
    outpheno = struct2table(vertcat(outpheno{:}));
catch
    outpheno = struct2table(vertcat(outpheno{:}), "AsArray", true);
end

end

%% ------------------------------------------------------------------------
function tout = mergePhenos(tin, opts)
% merges all queried phenotypes into one unified trait. This is useful for
% categorical traits where multiple data-fields are used for the definition
% of the target phenotype (e.g. in-patient and self-reported resources).

if numel(tin) < 2, tout = tin; return; end

tout = struct;
tout.tag = tin{1}.tag;
tout.numericFlag = tin{1}.numericFlag; % first trait numeric flag is only used

tout.eid = cellfun(@(x)x.eid, tin, 'uni', false);
tout.eid = unique(vertcat(tout.eid{:}));
[tout.rawUKB, tout.termMeaning, ...
    tout.source] = deal(strings(numel(tout.eid), 1));

for j = 1:numel(tin)
    if j < 2
        tout.info.basket = tin{j}.info.basket;
        tout.info.date = tin{j}.info.date;
        tout.info.df = tin{j}.info.df;
        try tout.info.desc = tin{j}.info.desc; catch, end
    else
        tout.info.basket = union(tout.info.basket, tin{j}.info.basket);
        tout.info.df = union(tout.info.df, tin{j}.info.df);
        try tout.info.desc = union(tout.info.desc, tin{j}.info.desc); catch, end
    end
    
    [idx1, ih2] = ismember(tout.eid, tin{j}.eid); ih2(ih2 < 1) = [];
    fi = ["rawUKB", "source", "termMeaning"];
    for k = 1:numel(fi)
        tout.(fi(k))(idx1) = tin{j}.(fi(k))(ih2) + "," + tout.(fi(k))(idx1);
    end
    
end

tout.rawUKB = regexprep(tout.rawUKB, '(^,$|,$|^,)', '');
tout.termMeaning = regexprep(tout.termMeaning, '(^,$|,$|^,)', '');
tout.source = regexprep(tout.source, '(^,$|,$|^,)', '');
tout.info.df = join(tout.info.df, ",");
tout.info.basket = join(tout.info.basket, ",");

if opts.save
    UKB_STRUCT_ALL = tout;
    savename = matlab.lang.makeValidName(UKB_STRUCT_ALL.tag);
    
     % if already exists, append date time to the file name.
    if isfile(fullfile(opts.phenodir, savename + ".mat"))

        oans = input(sprintf("\t%s is already present, overwrite? [y/N]: ", savename), "s");
        oans = lower(string(oans));
        
        if oans == "" || oans == "n"
            savename = savename + ...
                string(string(datetime)).replace(["-"," ", ":"], ".");
            fprintf("\t\tnew name was generated %s\n", savename)
        end
        
    end
    savename = savename + ".mat";
    if isfile(fullfile(opts.phenodir, savename))
        fprintf("\t\toverwriting %s \n", savename)
        
    % else % still file is present
    %     fprintf('cannot save file %s to %s\n', savename, opts.phenodir)
    end

    save(fullfile(opts.phenodir, savename), "UKB_STRUCT_ALL");
end

end % END

%% ------------------------------------------------------------------------
function tdic = getTimeDatedfs(tdic, d2t, dict, opts)

if opts.ukbrap
    tdic.Link = tdic.linkout;
else
    if ~any(colnames(tdic) == "Link")
        tdic.Link = "https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=" + tdic.FieldID;
    end
end

[f1, f2] = ismember(tdic.FieldID, d2t.df);
tdic.dateField(:) = strings(height(tdic), 1);
tdic.dateField(f1) = d2t.date(f2(f1));

% add non-participant date-fields, e.g. hesin_oper, gp_clinical etc.
% note that there may be different date valuetypes per each non-participant
% entity, which should be defined carefully. So, here (as of 26SEP2024), if
% only on Date valuetype is peresent, that field is picked. Otherwise an
% error will be generated.
if opts.ukbrap

    other_ent_idx = ~ismember(tdic.entity, "participant");
    if any(other_ent_idx)
        tdic2 = tdic(other_ent_idx, :);
        uent = unique(tdic2.entity);
        tdic3 = cell(numel(uent), 1);
        for k = 1:numel(uent)
            idx = tdic2.entity == uent(k);
            tmp_df = tdic2(idx, :);
    
            % find date-fields for this entity in dictionary file
            tmp_dt = dict(dict.entity == uent(k) & dict.ValueType == "Date", :);
            if height(tmp_dt) > 1
                error("multiple date-fields for entity %s was found!", uent(k))
            elseif isempty(tmp_dt)
                tdic3{k} = tmp_df;
                continue % no date-fields in this entity
            end

            tmp_df.dateField(:) = tmp_dt.name;
            tdic3{k} = tmp_df;
            clear tmp_df tmp_dt
        end
        tdic3 = vertcat(tdic3{:});
        tdic(other_ent_idx, :) = tdic3;
    end

end

% for fields data fields missing from d2t we check UKBB website to
% find the corresponding date
if any(~f1)
    tdic2 = tdic(~f1, :);
    ulinks = unique(tdic2.Link);

    for k = 1:numel(ulinks)
        if ismissing(ulinks(k)), continue; end

        wb = webread(ulinks(k), weboptions("Timeout", 5e4));
        wb = splitlines(string(wb));
        idx1 = find(wb.contains("Related Data-Field"));
        idx2 = find(wb == "</div>");
        idx2= idx2(find(idx2 > idx1, 1));
        wb = wb(idx1:idx2);
        wb = htmlTree(wb);
        wb = wb.extractHTMLText;
        wb(wb == "") = [];
        if wb(1).startsWith("0") % no related data-fields
            continue
        end
        wb(1) = [];
        wb = split(wb, whitespacePattern(2));

        % find date if any
        idx = contains(lower(wb(:, 2)), "date") & contains(lower(wb(:, 3)), "date") & contains(lower(wb(:, 3)), "current field");
        
        if any(idx)
            fix = tdic2.Link == ulinks(k);
            tdic2.dateField(fix) = wb(idx, 1); % date data-field associated with current data field
        end

    end

    tdic(~f1, :) = tdic2;
end


% add time associated data fields information
df = tdic.dateField;

if opts.ukbrap
    tdic.tmp = tdic.dateField;
    tdic.aa = tdic.array; tdic.ii = tdic.instance;
    tdic = fillmissing(tdic, "constant", -1, "DataVariables", ["aa", "ii"]);
    df = unique(df);
    N = numel(df);
else
    N = height(tdic);
end

tdic.dateField = cell(height(tdic), 1);
for k = 1:N
    if df(k) == "" || ismissing(df(k)), continue; end
    tmp = dict(dict.FieldID == df(k), :);

    % append data transformation funcs
    tmp.fnc = strings(height(tmp), 1);
    tmp.fnc(ismember(lower(tmp.ValueType), ["date", "time"])) = "datetime";
    tmp.fnc(ismember(lower(tmp.ValueType), ["continuous", "integer"])) = "double";

    if opts.ukbrap
        tmp.aa = tmp.array; tmp.ii = tmp.instance;
        tmp = fillmissing(tmp, "constant", -1, "DataVariables", ["aa", "ii"]);
        for j = 1:height(tmp)
            fidx = tdic.tmp == df(k) & tdic.aa == tmp.aa(j) & tdic.ii == tmp.ii(j);
            if ~any(fidx)
                % such as 41202 which has no array, but 41262 has 80
                % arrays in RAP dictionary. In this case we should merge
                % arrays similar to 41202 (| delimited)
                fidx = tdic.tmp == df(k);
                assert(sum(fidx) == 1, "something went wrong with mapping fields to data-fields!")
                tdic.dateField{fidx} = tmp(:, setdiff(colnames(tmp), ["aa", "ii"]));

                fprintf("\tfield:%s --> (correponds to) --> many date field:%s\n", ...
                    tdic.FieldID(fidx) + " " + tdic.Field(fidx), ...
                    tmp.FieldID(1) + " " + tmp.Field(1))

                break
            end
            tdic.dateField{fidx} = tmp(j, setdiff(colnames(tmp), ["aa", "ii"]));

            % disp
            fprintf("\tfield:%s --> (correponds to) --> date field:%s\n", ...
                tdic.FieldID(fidx) + " " + tdic.Field(fidx), ...
                tmp.FieldID(j) + " " + tmp.Field(j))
        end

    else
        tdic.dateField{k} = tmp;
    
        % disp
        fprintf("\tfield:%s --> (correponds to) --> date field:%s\n", ...
            tdic.FieldID(k) + " " + tdic.Field(k), ...
            tmp.FieldID + " " + tmp.Field)
    end
end

if opts.ukbrap
    tdic(:, ["tmp", "aa", "ii"]) = [];
end

end % END

%% ------------------------------------------------------------------------
function [df, opts] = rapChunkMerger(ds, opts)
% each chunk of UKB-RAP data has its own eid with possible different sizes.
% So all of them 

ds = cellfun(@(x)struct2table(x), ds, uni=false);

% compress tables if necessary
% merge duplicated eids
ct = 1;
for k = 1:numel(ds)
    if numel(unique(ds{k}.eid)) ~= numel(ds{k}.eid)
        if opts.surv % don't merge codes, since it would disrupt dt/df relations
            opts.schunk{ct} = ds{k};
            ds{k} = [];
            ct = ct + 1;
        else
            ds{k} = RAPcompressBasket(ds{k});
        end
    end
end
ds(cellfun(@isempty, ds)) = [];

df = ds{1}; ds{1} = [];
df = rmmissing(df, 2, MinNumMissing=height(df));

for k = 2:numel(ds) 
    ds{k} = rmmissing(ds{k}, 2, MinNumMissing=height(ds{k}));
    if all(ismember(ds{k}.eid, df.eid)) && all(ismember(df.eid, ds{k}.eid))
        df = rapChunkMerger_tab(df, ds{k});
    else
        df = outerjoin(df, ds{k}, Keys="eid", MergeKeys=true);
    end
    ds{k} = [];
end
df = rmmissing(df, 2, MinNumMissing=height(df));
df = table2struct(df, "ToScalar", true);

end % END

%% ------------------------------------------------------------------------
function df = rapChunkMerger_tab(df, ds)
% a subfunction for rapChunkMerger, faster than outerjoin for simpler cases

[~, idx] = ismember(df.eid, ds.eid);
ds.eid = [];
df = [df, ds(idx, :)];

end

%% ------------------------------------------------------------------------
function df = rapPruneDelimiter(df, codes, threads)
% some UKB-RAP fields are | or comma delimited. This function removes other
% codes and keeps only those in opts.categorical.

if isempty(codes), return; end
codes = string(codes); 

if ~isnan(threads) && isempty(gcp("nocreate"))
    parpool("Processes", threads);
end

if isnan(threads)
    dopar = false;
else
    dopar = true;
end

% build a single regex that matches any one of your codes, but only if
% it sits between start/| or |/end (so we don’t accidentally match substrings)
% tokPattern = strjoin(codes, "|");         % e.g. "D125|D509|D649|…"
% pattern = "(?<=^|\|)(" + tokPattern + ")(?=$|\|)";

% 1) escape delimiters for a regexp character class ---------------
dclass = regexptranslate('escape', "|,");   % -> '\|;'

% 2) escape the codes themselves and join them with "|" ------------
codes_esc = regexptranslate('escape', codes);
alternation = strjoin(codes_esc, "|");        % 'K75|K94\.1'

% 3) build the full pattern ---------------------------------------
pattern = "(?<=^|[" + dclass + "])" + alternation + "(?=[" + dclass + "]|$)";

if dopar
    tmp = cell(width(df), 1);
    parfor k = 1:width(df)
        col = df.(k);
        if ~isstring(col) || ~any(contains(col, ["|", ","])),  continue;  end
    
        % only look at non‐missing, non‐empty, pipe‐containing rows
        mask = ~ismissing(col) & col~="" & contains(col, ["|", ","]);
        data = col(mask);
    
        % REGEXP with 'match' returns a cell array {Nx1} of string arrays:
        matches = regexp(data, pattern, "match");
    
        out = strings(size(data));
        has = ~cellfun(@isempty, matches);
        matches = matches(has);

        % join the found tokens with commas, leave others missing
        if iscell(matches)
            matches = cellfun(@(x)x.join(","), matches);
        else
            matches = join(matches, ",");
        end
    
        % join the found tokens with commas, leave others missing
        out(has) = matches;
        out(~has) = missing;
    
        % write back into the table
        col(mask) = out;
        tmp{k} = col;
    end
    
    change_idx = find(~cellfun(@isempty, tmp));
    for k = 1:numel(change_idx)
        df.(change_idx(k)) = tmp{change_idx(k)};
    end
    clear tmp

else

    for k = 1:width(df)
        col = df.(k);
        if ~isstring(col) || ~any(contains(col, ["|", ","])),  continue;  end
    
        % only look at non‐missing, non‐empty, pipe‐containing rows
        mask = ~ismissing(col) & col~="" & contains(col, ["|", ","]);
        data = col(mask);
    
        % REGEXP with 'match' returns a cell array {Nx1} of string arrays:
        matches = regexp(data, pattern, "match");
    
        out = strings(size(data));
        has = ~cellfun(@isempty, matches);
        matches = matches(has);

        % join the found tokens with commas, leave others missing
        if iscell(matches)
            matches = cellfun(@(x)x.join(","), matches);
        else
            matches = join(matches, ",");
        end

        out(has) = matches;
        out(~has) = missing;
    
        % write back into the table
        col(mask) = out;
        df.(k)  = col;
    end

end


% if isempty(codes), return; end
% 
% for k = 1:width(df)
%     if isstring(df.(k)) && any(df.(k).contains("|"))
%         tmp = df.(k);
%         idx = ~(ismissing(tmp) | tmp == "");
%         tmp = arrayfun(@(x)x.split("|"), tmp(idx), uni=false);
%         % tmp = cellfun(@(x) join(codes(ismember(codes, x)), ","), tmp);
%         tmp = cellfun(@(x) rapPruneDelimiter_join(x, codes), tmp);
%         df.(k)(idx) = tmp;
%     end
% end

end % END

%% ------------------------------------------------------------------------
function str = rapPruneDelimiter_join(x, codes)
% to avoid errors when codes is scalar
    idx = ismember(codes, x);
    if any(idx)
        str = join(codes(idx), ",");
    else
        str = missing;
    end

end % END

%% ------------------------------------------------------------------------
function df = RAPcompressBasket(df)
% some RAP datasets are full 'entity' (see dx_extract_dataset.m function),
% e.g. hesin_oper, which contain multiple records per each eid. To avoid
% extra memory/computational time, we merge codes in a similar way as to
% what is done in DNANexus (| delimited). Note that this is not appropriate
% with 'surv' flag.

df = fillmissing(df, "constant", "", "DataVariables", @isstring);
df = groupsummary(df, "eid", @(x)join(unique(x), "|"));
df.GroupCount = [];
df.Properties.VariableNames = erase(colnames(df), textBoundary("start") + "fun1_");

% remove trailing/starting |
str_cols = varfun(@isstring, df);
for k = 1:width(df)
    if str_cols.(k) 
        if any(df.(k).contains(textBoundary("start") + "|"))
            df.(k) = erase(df.(k), textBoundary("start") + "|");
        end
    end
end

end % END

%% ------------------------------------------------------------------------
function [df, eid, matchidx] = filterChunkData(df, eid, opts)
% filter matching string vectors to large inpatient data

if ~isnan(opts.threads) && isempty(gcp("nocreate"))
    parpool("Processes", opts.threads);
end

if isnan(opts.threads)
    dopar = false;
else
    dopar = true;
end

matchidx = ~ismissing(df);
nCols = size(df, 2);
cats = opts.categorical;

% @23SEP2024: take care of | delimiter
hasPipe = false(nCols, 1);
for k = 1:nCols
    
    % @29MAY2025: some cases where pipe can be "," are in vitamins where on RAP
    % it's treated as integer, while the format is as [double, ... double], so
    % starting/trailing [] should be removed first
    idx_comma_sep = startsWith(df.(k), "[");
    if any(idx_comma_sep)
        df.(k)(idx_comma_sep) = extractBetween(df.(k)(idx_comma_sep), textBoundary("start") + "[", "]" + textBoundary("end"));
    end

    hasPipe(k) = any(contains(df.(k), ["|", ","]));

end
hasPipe_cols = colnames(df);
hasPipe_cols(~hasPipe) = [];

% first loop to reduce the size of df 
if ~dopar
    
    % matchidx = ~ismissing(df);
    for k = 1:nCols
    
        % if any(contains(df.(k)(matchidx(:, k)), "|"))
        if hasPipe(k)
            fidx = matchPipe(df.(k)(matchidx(:, k)), cats);
        else
            fidx = ismember(df.(k)(matchidx(:, k)), cats);
        end
    
        matchidx(matchidx(:, k), k) = fidx;
    end

    matchidx = any(matchidx, 2);

else % parallel
    
    % df = tall(df);
    % matchidx = ~ismissing(df);
    % for k = 1:nCols
    %     if hasPipe(k)
    %         fidx = contains(df.(k)(matchidx(:, k)), patt);
    %     else
    %         fidx = ismember(df.(k)(matchidx(:, k)), cats);
    %     end
    % 
    %     matchidx(matchidx(:, k), k) = fidx;
    % end
    % 
    % matchidx = any(matchidx, 2);

    tempRes = cell(1, nCols);  

    parfor k = 1:nCols
        currentMask = matchidx(:, k); % Logical mask (read-only slice)
        outCol = false(size(currentMask));
        dataSubset = df.(k)(currentMask);

        % if any(contains(dataSubset, "|"))
        if hasPipe(k)
            fidx = matchPipe(df.(k)(matchidx(:, k)), cats);
        else
            fidx = ismember(dataSubset, cats);
        end
        outCol(currentMask) = fidx;
        tempRes{k} = outCol;
    end

    % Reassemble the results into matchidx.
    matchidx = false(size(matchidx));
    for k = 1:nCols
        matchidx(:, k) = tempRes{k};
    end

    matchidx = any(matchidx, 2);

end

df(~matchidx, :) = [];

if istall(df)
    matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
    df = gather(df);
    matchidx = gather(matchidx);
    matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);
end

eid(~matchidx) = [];

% second round loop
df = rmmissing(df, 2, 'MinNumMissing', height(df));
nCols = size(df, 2);

% @23SEP2024: take care of | delimiter
hasPipe = false(nCols, 1);
hasPipe(ismember(colnames(df), hasPipe_cols)) = true;
% for k = 1:nCols
%     hasPipe(k) = any(contains(df.(k), "|"));
% end

if ~dopar

    matchidx = ~ismissing(df);
    for k = 1:nCols
        if hasPipe(k)
            fidx = matchPipe(df.(k)(matchidx(:, k)), cats);
        else
            fidx = ismember(df.(k)(matchidx(:, k)), cats);
        end
    
        matchidx(matchidx(:, k), k) = fidx;
    end

else % parallel

    % if ~istall(df), df = tall(df); end
    % matchidx = ~ismissing(df);
    % for k = 1:nCols
    %     if hasPipe(k)
    %         fidx = contains(df.(k)(matchidx(:, k)), patt);
    %     else
    %         fidx = ismember(df.(k)(matchidx(:, k)), cats);
    %     end
    % 
    %     matchidx(matchidx(:, k), k) = fidx;
    % end
    
    matchidx = ~ismissing(df);
    tempRes = cell(1, nCols);  

    parfor k = 1:nCols
        currentMask = matchidx(:, k); % Logical mask (read-only slice)
        outCol = false(size(currentMask));
        dataSubset = df.(k)(currentMask);

        % if any(contains(dataSubset, "|"))
        if hasPipe(k)
            fidx = matchPipe(df.(k)(matchidx(:, k)), cats);
        else
            fidx = ismember(dataSubset, cats);
        end
        outCol(currentMask) = fidx;
        tempRes{k} = outCol;
    end

    % Reassemble the results into matchidx.
    matchidx = false(size(matchidx));
    for k = 1:nCols
        matchidx(:, k) = tempRes{k};
    end

end

if istall(df)
    matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
    df = gather(df);
    matchidx = gather(matchidx);
    matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);
end

df = fillmissing(df, 'constant', "", 'MissingLocations', ~matchidx);
df = standardizeMissing(df, "");

matchidx = any(matchidx, 2);
df(~matchidx, :) = [];

df = rmmissing(df, 2, 'MinNumMissing', height(df));

end % END

%%
function fidx = matchPipe(dataSubset, cats)
fidx = false(size(dataSubset));

% Identify entries that contain the delimiter.
hasPipe = contains(dataSubset, ["|", ","]);

% For entries with no '|', a direct ismember is sufficient.
fidx(~hasPipe) = ismember(dataSubset(~hasPipe), cats);

if any(hasPipe)

    % 1) escape delimiters for a regexp character class ---------------
    dclass = regexptranslate('escape', "|,");   % -> '\|;'

    % 2) escape the codes themselves and join them with "|" ------------
    codes_esc = regexptranslate('escape', cats);
    alternation = strjoin(codes_esc, "|");        % 'K75|K94\.1'

    % 3) build the full pattern ---------------------------------------
    pattern = "(^|[" + dclass + "])(" + alternation + ")([" + dclass + "]|$)";

    % Vectorized split:
    % The result is an array where each row corresponds to an entry and each
    % column to a token. The number of columns will be the maximum tokens found.   
    dataSubset = dataSubset(hasPipe);
    matchCells = regexp(dataSubset, pattern, 'once');
    if ~iscell(matchCells), matchCells = {matchCells}; end % if matchCells is scalar
    fidx(hasPipe) = ~cellfun(@isempty, matchCells);
end

end % END