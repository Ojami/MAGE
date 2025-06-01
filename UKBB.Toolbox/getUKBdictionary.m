function [data_dict, codings] = getUKBdictionary(opts)
% getUKBdictionary fetches data dictionaries/encodings from UK
% Biobank website
% (https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide)
% 
% Oveis Jamialahmadi, University of Gothenburg, March 2019.

% @26APR2024: data_coding/codings tables have been removed by UK Biobank
% and replaced by schema files. 'dicturl', 'encurl' are deprecated now.
% 
% @19SEP2024: 'ukbrap' flag was added to fetch/parse dictionaries from
% UKB-RAP. Note that dxtoolkit is required. 

arguments
    % opts.dicturl {mustBeTextScalar} = "https://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.tsv"
    % opts.encurl {mustBeTextScalar} = "https://biobank.ctsu.ox.ac.uk/~bbdatan/Codings.tsv";
    opts.updateonly (1,1) logical = false % only updates dictionaries/encodings
    opts.updatefreq (1,1) double = 30 % update frequency of dictionaries/encodings in days
    opts.dir {mustBeTextScalar} = fullfile(fileparts(which('getUKBdictionary.m')), "UKB_DICTIONARY") % output dir

    %19SEP2024
    opts.ukbrap (1,1) logical = true % use dictionaries from UKB-RAP instead
end

if opts.ukbrap
    opts.dname = "data_dict_codings_rap.mat";
else
    opts.dname = "data_dict_codings.mat";
end

if ~isfolder(opts.dir)
    mkdir(opts.dir)
    fetchflag = true;
elseif ~isfile(fullfile(opts.dir, opts.dname))
    fetchflag = true;
else
    % check if needed to update database
    checkdate = load(fullfile(opts.dir, opts.dname), 'datev').datev;
    diffdays = caldays(between(checkdate, datetime('now'), 'Days'));
    if diffdays >= opts.updatefreq
        fetchflag = true;
    else
        fetchflag = false;
    end
end

if fetchflag || opts.updateonly
    
    if opts.ukbrap
        [~, project_id] = runbash("dx env | grep project- | awk -F '\t' '{print $2}'");
        [~, dataset_name] = runbash("dx ls " + project_id(end));
        dataset_name(~dataset_name.endsWith(".dataset")) = [];
        dataset_name = natsort(dataset_name);
        dataset_name = dataset_name(end); % latest

        [~, dataset_id] = runbash("dx describe " + project_id + ":" + dataset_name);
        dataset_id(~dataset_id.startsWith("ID")) = [];
        dataset_id = dataset_id.extractAfter("ID").strip;

        % fetch dictionaries
        dfiles = getfilenames(opts.dir, "csv", "fullpath", true).csv;
        if ~isempty(dfiles)
            for k = 1:numel(dfiles), delete(dfiles(k)); end
        end
        runbash("dx extract_dataset " + project_id + ":" + dataset_id + ...
                " -ddd --delimiter ',' --output '" + ...
                makeWSLpath(opts.dir) + "'", wait=true, verbose=true);
        
        dfiles = getfilenames(opts.dir, "csv").csv;
        dfiles(~dfiles.startsWith(dataset_name)) = [];
        dfiles = fullfile(opts.dir, dfiles);

        data_dict = dfiles(dfiles.endsWith(".data_dictionary.csv"));
        data_dict = readtable(data_dict, TextType="string");
        data_dict = rmmissing(data_dict, 2, MinNumMissing=height(data_dict));

        % create FieldID
        data_dict.FieldID = data_dict.name;
        idx = contains(data_dict.name, textBoundary("start") + ...
            "p" + digitsPattern + (textBoundary("end") | "_"));
        data_dict.FieldID(idx) = extractBetween(data_dict.FieldID(idx), ...
            textBoundary("start") + "p", ...
            (textBoundary("end") | "_"));

        data_dict = renamevars(data_dict, [ "title", "units", "coding_name", "type"], ...
            ["Field", "Units", "Coding", "ValueType"]);
        data_dict.Coding = extractAfter(data_dict.Coding, "data_coding_");
        data_dict.Coding(ismissing(data_dict.Coding)) = "";
        data_dict.ValueType(data_dict.ValueType == "date") = "Date";
        data_dict.ValueType(data_dict.ValueType == "datetime") = "Time";
        data_dict.ValueType(data_dict.ValueType == "float") = "Continuous";
        data_dict.ValueType(data_dict.ValueType == "integer") = "Integer";
        data_dict.ValueType(data_dict.ValueType == "string") = "Text";

        % extract instances/arrays
        data_dict.instance(:) = nan; data_dict.array(:) = nan;
        idx = data_dict.name.contains("_i" + digitsPattern + (textBoundary("end")|"_"));
        data_dict.instance(idx) = extractBetween(data_dict.name(idx), "_i", textBoundary("end")|"_");
        
        idx = data_dict.name.contains("_a" + digitsPattern + (textBoundary("end")|"_"));
        data_dict.array(idx) = extractBetween(data_dict.name(idx), "_a", textBoundary("end")|"_");

        
        % read codings
        codings = dfiles(dfiles.endsWith(".codings.csv"));
        op = detectImportOptions(codings, TextType="string");
        op.VariableTypes = repmat("string", 1, numel(op.VariableTypes));
        codings = readtable(codings, op);
        codings = rmmissing(codings, 2, MinNumMissing=height(codings));
        codings = renamevars(codings, ["coding_name", "meaning", "code"], ["Coding", "Meaning", "Value"]);
        codings.Coding = extractAfter(codings.Coding, "data_coding_");
        codings.Value = string(codings.Value);


    else

        data_dict = readtable("https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=1", weboptions=weboptions("Timeout", 5e4), FileType="text", TextType="string");
        cd = [5:8, 11, 12, 20];
        codings = cell(numel(cd), 1);
        for k = 1:numel(cd)
            codings{k} = readtable("https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=" + cd(k), ...
                weboptions=weboptions("Timeout", 5e4), FileType="text", TextType="string", ...
                ReadVariableNames=true, NumHeaderLines=0);
            codings{k} = convertvars(codings{k}, ["value", "encoding_id"], @string);
            codings{k} = codings{k}(:, ["encoding_id", "value", "meaning"]);
    
        end
        codings = vertcat(codings{:});
        codings = renamevars(codings, ["encoding_id", "meaning", "value"], ["Coding", "Meaning", "Value"]);
    
        data_dict = renamevars(data_dict, ["field_id", "title", "main_category", "units", "encoding_id", "value_type"], ...
            ["FieldID", "Field", "Category", "Units", "Coding", "ValueType"]);
        data_dict = convertvars(data_dict, ["FieldID", "Category", "Coding", "ValueType"], @string);
        data_dict.Coding(data_dict.Coding == "0") = "";
    
        % change data types
        dd.t = ["Categorical", "Categorical", "Categorical", "Continuous", "Date",...
            "Integer", "Text", "Time"];
        dd.c = string([22, 21, 101, 31, 51, 11, 41, 61]);
        [f1, f2] = ismember(data_dict.ValueType, dd.c);
        data_dict.ValueType(f1) = dd.t(f2(f1));
    end


    datev = datetime('now');

    save(fullfile(opts.dir, opts.dname), 'data_dict', 'codings', 'datev')
   
    if opts.updateonly % only update database
        return
    end
else
    data_dict = load(fullfile(opts.dir, opts.dname), 'data_dict').data_dict;
    codings = load(fullfile(opts.dir, opts.dname), 'codings').codings;
end

end % END
