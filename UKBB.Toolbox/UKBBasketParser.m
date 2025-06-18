function UKBBasketParser(infile, opts)
% Uses inputUKBfile (output file from ukbconv, in csv or txt format) and
% stores all variables in a new folder for future use.
% 
% Oveis Jamialahmadi, University of Gothenburg, Sept 2019.
% @04/07/2022: major modifications were made.
% 
% @20SEP2024: major changes were made to parse output from
%             dx_extract_dataset (UKB-RAP).
% @25SEP2024: 'entity' flag was added. If true, 'infile' is treated as a
%             separate dataset, and hence, 'vsize' is not applied. This
%             flag is used for UKB entity tables on UKB-RAP, e.g.
%             hesin_oper or gp_clinical datasets, where each dataset should
%             be saved together in a separate chunk. Name of this chunk
%             will be exactly as the name of infile (Note: this flag is
%             intended to be used along with dx_extract_dataset which
%             handels non-"participant" entities internally using this
%             flag).

arguments
    infile {mustBeFile}
    % opts.tall (1,1) logical = false % DEPRECATED
    opts.parallel (1,1) logical = true % parfor
    opts.vsize (1,1) double = 30 % number of variables per each file
    opts.threads (1,1) double = 30 % thread for parallelizing

    %@20SEP2024
    opts.ukbrap (1,1) logical = false % files in infile are from ukbrap?
    opts.path {mustBeTextScalar} % output path. Will create if does not exist.
    opts.verbose (1,1) logical = true

    % is the input file an entity file? e.g. hesin_oper
    opts.entity (1,1) logical = false 
    opts.parquet (1,1) logical = false % only if entity is true
end

if opts.parallel && isempty(gcp('nocreate'))
    % distcomp.feature( 'LocalUseMpiexec', false );
    parpool('Processes', opts.threads); % Use parallel computations
end

[~, inputUKBname, ext] = fileparts(infile);
if strcmp(ext, ".csv")
    sep = ",";
elseif strcmp(ext, ".txt")
    sep = char(9);
else
    error('UKBBasketParser:only txt or csv files are supported!')
end

if isfield(opts, "path")
    if ~isfolder(opts.path), mkdir(opts.path); end
else
    if ~isfolder("UKBFileParser_" + inputUKBname)
         mkdir("UKBFileParser_" + inputUKBname)
    end
    opts.path = fullfile(pwd, "UKBFileParser_" + inputUKBname);
end

rawfile = tabularTextDatastore(infile, 'TextType', 'string', ...
    'NumHeaderLines', 0, 'ReadVariableNames', true, 'Delimiter', sep, ...
    'VariableNamingRule', 'preserve');

%@25SEP2024: only replace %f columns with %q
% rawfile.SelectedFormats = repmat({'%q'}, numel(rawfile.SelectedFormats), 1);
rawfile.SelectedFormats = replace(rawfile.SelectedFormats, '%f', '%q');

if opts.ukbrap
    %@20SEP2024: support for data from UKB-RAP: output of
    %dx_extract_dataset

    % check if variables come from extract_dataset: entity.eid
    eid_idx = contains(rawfile.VariableNames, ".eid");
    if any(eid_idx)
        rawfile.VariableNames = extractAfter(rawfile.VariableNames, ".");
    end

    % find eid
    eid_idx = ismember(rawfile.VariableNames, ["xeid", "eid"]);
    assert(any(eid_idx), "eid column is missing!")
    rawfile.VariableNames{eid_idx} = 'eid';
    rawfile.SelectedFormats{eid_idx} = '%f';

else
    rawfile.VariableNames = "x" + replace(string(rawfile.VariableNames), ["-", "."], "_");
    rawfile.VariableNames = replace(rawfile.VariableNames, "xeid", "eid");
end

ukbvars = rawfile.VariableNames; % Get variables

if opts.entity
    seq = [0, numel(ukbvars)];
else
    seq = unique([0:opts.vsize:numel(ukbvars), numel(ukbvars)]); % Range of saving variables together
end
variableSaver = ({});

variableMapper = strings(numel(ukbvars) - 1, 2);

opts.cnum = 0; % maximum chunk number of existing variableMapper: for ukbrap
if opts.ukbrap
    if isfile(fullfile(opts.path, 'variableMapper.mat'))
        opts.vmap = string(load(fullfile(opts.path, 'variableMapper.mat')).variableMapper);
        opts.vdate = load(fullfile(opts.path, 'variableMapper.mat')).datev;
        opts.cnum = max(double(opts.vmap(:, 1)));
    end
end

if opts.entity
    ukbvars = string(ukbvars);
    idx_range = 1:numel(ukbvars);
    variableMapper(idx_range, 1) = inputUKBname;
    variableMapper(idx_range, 2) = ukbvars;

else

    start_cnt = 1;
    for k = 2:numel(seq)
        chunk_number = opts.cnum + k - 1;
        
        varchunk = string(ukbvars(seq(k-1)+1:seq(k)));
    
        if opts.ukbrap % eid should be present in each chunk
            varchunk = union(varchunk, "eid");
        else
            varchunk(ismember(varchunk, "eid")) = [];
        end
    
        variableSaver{k-1, 1} = varchunk;
        variableSaver{k-1, 2} = chunk_number;
        
        varrange = start_cnt:numel(varchunk)+start_cnt-1;
        variableMapper(varrange, 1) = repmat(chunk_number, numel(varchunk), 1);
        variableMapper(varrange, 2) = varchunk;
        start_cnt = numel(varchunk) + start_cnt;
    end
end

if opts.ukbrap
    % date per each chunk of data fetched from UKB-RAP
    datev = repmat(datetime('today'), size(variableMapper, 1), 1);
else
    datev = datetime('today');
end

if opts.cnum ~= 0 % update variableMapper
    variableMapper = [opts.vmap; variableMapper];
    datev = [opts.vdate; datev];
else
    
end
save(fullfile(opts.path, 'variableMapper.mat'), 'variableMapper', 'datev')

if opts.verbose, tt = tic; end

matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);

if opts.entity
    % read the whole entity dataset
    if opts.parallel
        rawfile = tall(rawfile);
    else
        rawfile = readall(rawfile);
    end
    
    % check if any colnames are in ins_index, arr_index or level and
    % convert to double
    int_cols = intersect(colnames(rawfile), ["ins_index", "arr_index", "level"]);
    if ~isempty(int_cols)
        for k = 1:numel(int_cols)
            rawfile.(int_cols(k)) = double(rawfile.(int_cols(k)));
        end
    end

    rawfile = standardizeMissing(rawfile, "");

    if opts.parallel
        rawfile = gather(rawfile);
    end
    
    if opts.parquet
        parquetwrite(fullfile(opts.path, "UKB_" + inputUKBname + ".parquet"),...
            rawfile, VariableCompression="gzip")
    else
        rawfile = table2struct(rawfile, 'ToScalar', true);
        save(fullfile(opts.path, "UKB_" + inputUKBname + ".mat"), '-struct', 'rawfile')
    end

else

    if opts.parallel; n = gcp('nocreate').NumWorkers; else; n = 0; end
    
    targetDir = opts.path;
    isukbrap = opts.ukbrap;
    verbosity = opts.verbose;
    parfor (k = 1:size(variableSaver, 1), n)
    % for k = 1:size(variableSaver, 1) %DEBUG
    
        if verbosity
            fprintf('%d of %d\n', k, size(variableSaver, 1))
        end
    
        chunkdata = rawfile;
        chunkdata.SelectedVariableNames = variableSaver{k, 1};
        chunkdata = readall(rawfile);
        chunkdata = standardizeMissing(chunkdata, "");
        
        if ~isukbrap
            if any(ismember(variableSaver{k, 1}, 'eid'))
                UKB_eid = double(chunkdata.eid);
                saveDataChunk(fullfile(targetDir, 'UKB_eid'), UKB_eid, false)
                chunkdata.eid = [];
            end
        end
       
        chunkdata = table2struct(chunkdata, 'ToScalar', true);
        saveDataChunk(fullfile(targetDir, "UKB_" + variableSaver{k, 2} + ".mat"), chunkdata)
    end
end

matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);

if opts.verbose
    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('Done in %s (hr:min:sec)\n', tt)
end

end % END

%% subfunctions -----------------------------------------------------------
function saveDataChunk(name, chunkdata, sflag)
if nargin < 3
    sflag = true;
end

if sflag
    save(name, '-struct', 'chunkdata')
else
    UKB_eid = chunkdata;
    save(name, 'UKB_eid')
end
end
