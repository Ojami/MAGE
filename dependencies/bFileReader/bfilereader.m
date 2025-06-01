function [out, opts] = bfilereader(bfile, opts)
% bfilereader (big file reader) implemenets some Java methods to filter
% (regex and numeric filtering) big delimited files (can also be gz
% compressed) in a fast and efficient manner.
% 
% INPUTS:
%   bfile: path to input file: can be dellimited or gz compressed.
% 
% OPTIONAL:
%   readAll:    a logical. If true, both 'filter' and 'pattern'
%               options are ignored and whole file content will be read.
%               Note that this may cause out of memory error if input file
%               does not fit into memory.
% 
%   sep:        column delimiter (e.g. "\t" or " "). If left empty, it will
%               be automatically identified from the first line of the
%               input file. Note that the function first tries to check
%               common delimiters: char(9), " " and "\t". if failed, it
%               tries , and ; delimiters. If again fails, function throws
%               an error printing the first line of file.
% 
%   pattern:    string array/scalar of patterns for filtering. Pattern(s)
%               can be also Java regular expressions (regex). Unlike
%               'operator' and 'filter', 'col' can be a vector with
%               'pattern' option. Either of 'pattern' or 'filter' can be
%               used in each call to the function.
% 
%   patternCol: column name(s) or index(indices) for which pattern matching
%               will be applied. If 'pattern' is not provided, function
%               only reads the provided columns. If column names are used,
%               first row is treated as header (overrides 'header' flag).
%               Default is NaN: in this case, function used all columns for
%               either pattern matching (if 'pattern' is not empty) or
%               reading the whole file.
% 
%   multiCol:   a logical for matching each pattern to the corresponding
%               collumn in 'patternCol'. If set to true, 'patternCol' and
%               'pattern' must have the same size. This is useful when
%               there are multiple patterns to be matched over different
%               columns.
% 
%   filter:     the vector of cutoff values to be used for filtering
%               combined with 'operator' and 'filterCol' options. Default
%               is NaN: no numeric filtering to be applied.
% 
%   filterCol:  the vector of column name(s) or index(indices)
%               corresponding to values in 'filter'.
%
%   operator:   the vector of operators correponding to values in 'filter'
%               option. It can be a scalar operator if same operation
%               should be applied to all 'filter' values.
% 
%   extractCol: the vector of column names/indices to be extracted.
%               Default is NaN: return all columns after filtering/pattern
%               matching.
% 
%   header:     a logical to check whether first row contains header (true)
%               or file has no header (default: false). This flag is
%               used for 'return' as "table".
% 
%   skip:       a scalar double: number of lines to skip (default: 0). When
%               'header' is true, then first row is treated as header
%               (variable) names, then bfilereader automatically skips the
%               first (header) row and excludes it from processing. In this
%               case, the number of skipping lines is 'skip' + 1 (header
%               row). So, it is impotant that value used for 'skip' option
%               merely shows the number of comment/info lines in the
%               beginning of the input file, and not the header row (see
%               examples for more details).
% 
%   return:     the output format of returned rows, which can be either:
%                   - raw:        char vector output from the dependent 
%                                 Java class.
%                   - string:     a string array split into columns.
%                   - rawTable:   a table with returned rows and columns as
%                                 string.
%                   - table:      tries to further convert columns of
%                                 numeric data to double (default).
% 
%   parallel:   a logical (default: false). If set to true, function
%               uses multiple processor cores. Note that, this doesn't
%               necessarily mean an increase in the performance, as memory
%               overhead caused by parallel streaming may lower the
%               performance.
% 
%   summary:    to dipslay the first few lines of the file (default: off).
%               if set to "only", no processing will be done and these few
%               lines will be returned as output. This may be useful when
%               file is big enough that user wants to check the file
%               content before any decision on how to parse the file.
%   comment:    same as MATLAB readtable CommentStyle. bfilereader skips
%               lines starting with this comment style. Note, this option
%               overrides 'skip' option. 
% 
%   verbose:    to show warning/messages (on) or call the function
%               silently (off). Default is "timeOnly" showing the elapsed
%               time only.
% 
% OUTPUT: 
%   out:        read/filtered file content with the format set in
%               'return' option.
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, University of Gothenburg, July 2021.
% 
% Important notes: 
%                   - MATLAB R2019b or newer is needed!
%
%                   - bfilereader relies Java class bFileReaderDep.class.
%                   While bfilereader adds this dependent class to dynamic
%                   path using javaaddpath, you should consider adding this
%                   class to static path for better/more stable
%                   performance. For more help see doc javaclasspath.
% 
%                   - Java 8 is required (version -java).For more
%                   information, see: https://se.mathworks.com/help/compiler_sdk/java/configure-your-java-environment.html
% 
% @03/02/2022: a bug was fixed with single column in 'extractCol'.
% @03/11/2022: 'comment' options was added, and skips lines having this
%               pattern.
% @03/19/2022: a bug in 'sep' argument was fixed for 1 column files. 

arguments
    bfile {mustBeFile} % input file 
    opts.readAll (1, 1) logical = false; % reads whole file content.
    opts.sep {mustBeTextScalar} % if empty, delimiter is automatically identified.
    
    opts.pattern {mustBeText, mustBeVector} = "" % a string array/scalar containing the patterns to be used for filtering.
    opts.patternCol {mustBeVector} = NaN; % column names/indices for pattern search. Default is NaN: all columns.
    opts.multiCol (1, 1) logical = false; % match each pattern in 'pattern' to each column in 'patternCol', i.e. pattern(1) to col(1) & pattern(2) to col(2) ...
    
    opts.filter double {mustBeVector} = NaN; % the cutoff for which the "operator" filter will be applied.
    opts.filterCol {mustBeVector} = 1; % % column names/indices for filter values of numeric type.
    opts.operator {mustBeVector, mustBeMember(opts.operator, [">", ">=", "<", "<=", "==", ""])} = ""; % only for numeric filtering
    
    opts.extractCol {mustBeVector} = NaN; % column names/indices to be extracted. Default: NaN: return all columns after filtering/pattern matching. 
    
    opts.header (1, 1) logical = false; % is first row the header?
    opts.skip (1, 1) double = 0; % number of lines to skip. Default is 0: parse all rows in the input file (see docs above).
    opts.return {mustBeMember(opts.return, ["table", "rawTable", "raw", "string"])} = "table";
    opts.parallel (1, 1) logical = false; % use parallel (true) or sequential (false) stream.
    opts.summary {mustBeMember(opts.summary, ["on", "off", "only", "firstline", "linecount"])} = "off"; % "on": show first 6 lines of file, "only": don't process files after displaying summary, "firstline": only read the first line of file (after 'skip' lines).
    opts.verbose {mustBeMember(opts.verbose, ["on", "off", "timeOnly"])} = "off"; % "on": show warnings/messages, "timeOnly": only show elapsed time.
    opts.comment {mustBeTextScalar}

    % under development options
    opts.joinCol (1, 1) logical = false % joins (_)columns in 'patternCol' into a single column prior to pattern match. 
    opts.chunk (1, 1) double = 0 % mixture of parallel and serial processing. File's split into chunks and each chunk is processed in parallel (parallel flag must be set to true, default: not used)
    opts.extractRow {mustBeVector} = NaN % row indices to be extract. Only applicable to method getColumn.
end

%% check inputs -----------------------------------------------------------

%@11APR2025: a bug in new bfilereader with parallel/chunk
if opts.parallel
    disp("bfilereader::parallel was set to false due to a bug!")
    opts.parallel = false;
end

if opts.chunk ~= 0
    disp("bfilereader::chunk was set to false due to a bug!")
    opts.chunk = 0;
end

% check 'patternCol'/'filterCol': if col names are provided, 'header' is set to true
opts = assertCol(opts, "filterCol");
if isscalar(opts.patternCol) && isnumeric(opts.patternCol) && isnan(opts.patternCol) % i.e. all columns
    [opts.patternColAll, opts.patternColNumeric] = deal(true);
else
    opts = assertCol(opts, "patternCol");
    opts.patternColAll = false;
end

% check optional columns to be extracted
if ~(isnumeric(opts.extractCol) && any(isnan(opts.extractCol)))
    opts = assertCol(opts, "extractCol"); % else: retrun all cols
    opts.extractColAll = false;
else
    opts.extractColAll = true;
end

% if either of patternColNumeric or filterColColNumeric is false, 'header'
% flag shoud be set to true: file indeed has a header (variable names) row.
if any(~[opts.filterColNumeric, opts.patternColNumeric])
    opts.header = true;
end

% check if 'filter' and 'operator' have been set properly.
opts.operator = string(opts.operator);
opts.pattern = string(opts.pattern);
if ~opts.readAll % 'readAll' overrides these options
    if numel(opts.operator) > 1 && any(opts.operator == "")
        error('bfilereader: ''operator'' cannot contain empty operators!')
    elseif numel(opts.filter) > 1 && any(isnan(opts.filter))
        error('bfilereader: ''filter'' cannot contain NaN values!')
    elseif all(~isnan(opts.filter)) && any(opts.operator == "")
        error('bfilereader:''operator'' must be set with ''filter'' option!')
    elseif all(opts.operator ~= "") && any(isnan(opts.filter))
        error('bfilereader:''filter'' must be set with ''operator'' option!') 
    elseif numel(opts.pattern) > 1 && any(opts.pattern == "" | ismissing(opts.pattern)) % first condition: to skip default option
        error('bfilereader:''pattern'' cannot contain empty/missing cells!')
    end
end
%% prepare other details --------------------------------------------------
% get home (bFilReaderDep class) and working directories
% input file must have full path, it suffices to check its presence only in
% pwd
if exist(fullfile(pwd, bfile), 'file')
    bfile = fullfile(pwd, bfile);
end

hd = fileparts(which('bfilereader.m'));
% check if bFileReaderDep.class is present in bfilereader.m directory
if ~isfile(fullfile(hd, 'bFileReaderDep.class'))
    error('bfilereader: bFileReaderDep.class must be present in bfilereader.m directory!')
end
% check if bFileReaderDep.class has already added to dynamic/static path.
jpath = javaclasspath('-all');
if ~any(ismember(jpath, fullfile(hd))) % not found in jpath
    javaaddpath(fullfile(hd)); % see notes under 'important notes'.
end

% decide what to do: readAll, readCol (only specific columns), numeric filtering, or pattern match
opts.dict = table([">"; ">="; "<"; "<="; "=="], ...
        ["gt"; "ge"; "lt"; "le"; "eq"]);
if opts.readAll
    opts.method = "readAll";
elseif any(opts.operator == "") && all(opts.pattern == "") 
    opts.method = "readCol";
elseif all(opts.operator ~= "") && ~all(opts.pattern == "")
    opts.method = "match & filter";
elseif any(opts.operator == "")
    opts.method = "match";
else
    opts.method = "filter";
end


%% begin file proccessing -------------------------------------------------
reader = bFileReaderDep;

% check comment style
if isfield(opts, 'comment')
    start = 0;
    while true
        checkline = string(reader.readHeader(bfile, java.lang.Integer(start)));
        if ~startsWith(checkline, opts.comment) && checkline ~= ""
             % no reason to continue with start == 0: comment can't be even
             % found in first line
            if start ~= 0
               opts.skip = start; % comment lines
            end
            break
        end
        start = start + 1;
    end
end

% get file header and count rows ------------------------------------------
opts.line = string(reader.readHeader(bfile, java.lang.Integer(opts.skip))); % first line (for header)
if strcmp(opts.verbose, "on") || strcmp(opts.summary, "linecount")
    opts.lineCount = double(reader.lineCount(bfile)); % count number of rows
    out.linecount = opts.lineCount;
    out.firstline = opts.line;
    return
else
    opts.lineCount = 1;
end

% check file and fetch needed info: delimiter, header and column indices --
opts = getbFileInfo(opts);

% get summary (only first 6 lines) ----------------------------------------
if ~strcmp(opts.summary, "off")
    if strcmp(opts.summary, "firstline") % don't read lines after this.
        out = split(opts.line, opts.sep);
        return
    end
    
    fileSummary = strings(5, 1);
    if opts.header
        skip  = double(opts.skip) - 1;
    else
        skip = double(opts.skip);
    end
    
    for i = 1:numel(fileSummary)
        try
            fileSummary(i) = string(reader.readHeader(bfile, java.lang.Integer(skip+i)));
        catch
            break % exceeds line number
        end
    end
    fileSummary(fileSummary == "") = [];
    
    if strcmp(opts.return, 'raw')
        fileSummary = [opts.line;fileSummary];
    elseif strcmp(opts.return, 'string')
        fileSummary = split([opts.line;fileSummary], opts.sep);
    else
        if opts.header
            if strcmp(opts.verbose, "on"); disp('file first 5 rows:'); end
            fileSummary = array2table(split(fileSummary, opts.sep, 2), ...
                'VariableNames', split(opts.line, opts.sep));
        else
            if strcmp(opts.verbose, "on"); disp('file first 6 rows:'); end
            fileSummary = array2table(split([opts.line;fileSummary], opts.sep, 2));
            fileSummary.Properties.VariableNames = "Var" + (1:size(fileSummary, 2));
        end
    end
    if strcmp(opts.verbose, "on"); disp(fileSummary); fprintf('\n'); end
    if strcmp(opts.summary, "only")
        out = fileSummary;
        return
    end
end

% if method is "readCol" with all columns to be read, change the method to
% "readAll" this means both patternColAll and extractColAll are true amd
% 'extractRow' is nan. This also means, "readCol" can be invoked by setting
% 'patternCol', 'extractCol' or 'extractRow' options.
if strcmp(opts.method, "readCol") && opts.patternColAll && ...
        opts.extractColAll && isempty(opts.extractRow(1))
    opts.method = "readAll"; % which is faster than "readCol" with all cols
elseif opts.patternColAll && ~opts.extractColAll
    opts.patternCol = opts.extractCol; % use 'extractCol' option for "readCol"
end

if ~strcmp(opts.verbose, "off"); tic; end
out = ''; % if crashes, return this empty char
try
    switch opts.method
        case "readAll"
            out = reader.readAll(bfile, opts.sep, opts.skip, opts.parFlag);
        case "readCol" % --------------------------------------------------
            if numel(opts.patternCol) > 1 && strcmp(opts.parFlag, "s") && isempty(opts.extractRow(1))
                out = reader.getColumnBuffer(bfile, opts.sep, opts.patternCol, opts.skip, opts.parFlag);
            else
                out = reader.getColumn(bfile, opts.sep, opts.patternCol, opts.skip, opts.parFlag, opts.extractRow, java.lang.Integer(opts.chunk));
            end
        case "filter" % ---------------------------------------------------
            out = reader.filterCol(bfile, opts.sep, opts.filterCol, opts.operator, opts.filter, opts.skip, opts.parFlag, opts.extractCol);
        case "match" % ----------------------------------------------------
            if opts.multiCol
                out = reader.multiCompareToCols(bfile, opts.skip, opts.sep, opts.pattern, opts.patternCol, opts.parFlag, opts.extractCol);
            elseif opts.joinCol
                out = reader.compareToJoinedCols(bfile, opts.pattern, opts.sep, opts.patternCol, opts.skip, opts.parFlag, opts.extractCol, java.lang.Integer(opts.chunk));
            else
                if numel(opts.patternCol) == numel(opts.headerVars) % search over all columns
                    out = reader.compare(bfile, opts.pattern, opts.sep, opts.skip, opts.parFlag, opts.extractCol);
                else
                    out = reader.compareToCols(bfile, opts.pattern, opts.sep, opts.patternCol, opts.skip, opts.parFlag, opts.extractCol);
                end
            end
        case "match & filter" % -------------------------------------------
            if opts.multiCol
                out = reader.multiCompareFilterCol(bfile, opts.skip, opts.sep, opts.pattern,...
                    opts.patternCol, opts.operator, opts.filterCol, opts.filter, opts.parFlag, opts.extractCol);
            else
                out = reader.compareFilterCol(bfile, opts.skip, opts.sep, opts.pattern,...
                    opts.patternCol, opts.operator, opts.filterCol, opts.filter, opts.parFlag, opts.extractCol);
            end
    end
catch ME
    disp(ME.message)
    error('bfilereader:internal Java error')
end

clear reader
java.lang.System.gc()


%% make output MATLAB friendly --------------------------------------------
out = javaStr2MATLAB(out, opts);

if ~strcmp(opts.verbose, "off")
    elapseT = toc;
    fprintf('Elapsed time is %.3f seconds\n', elapseT)
end

end % END

%% subfunctions -----------------------------------------------------------
function opts = assertCol(opts, name)
% checks if 'patternCol' and 'filterCol' are of numeric or string type.
if isnumeric(opts.(name)) % column index(indices)
    if any(isnan(opts.(name)) | opts.(name) < 1)
        error('bfilereader:''%s'' cannot contain NaN or < 1 values!', name)
    end
    opts.(name + "Numeric") = true;
else % col names
    opts.(name + "Numeric") = false;
    opts.(name) = string(opts.(name));
    if any(opts.(name) == "" | ismissing(opts.(name)))
        error('bfilereader:''%s'' cannot contain empty/missing cells!', name)
    end
end
end

%% ------------------------------------------------------------------------
function opts = getbFileInfo(opts)
% prepares/gets file delimiter, column indices and header (first row)
% variable names (either if 'header' is true ot 'col' is a string array).
lineCount = opts.lineCount;
line = opts.line;

% check 'sep' delimiter ---------------------------------------------------
% @03/19/2022: commented out, since may return unexpected results
% especially with 1 column datasets.
% if opts.sep ~= "" % detect file delimiter
%     lineCols = line.split(opts.sep);
%     if numel(lineCols) == 1 % wrong delimiter
%         if strcmp(opts.verbose, "on")
%             fprintf('WARNING: wrong delimiter: %s\n', string(opts.sep))
%             fprintf('will try to find it automatically!\n')
%         end
%         opts.sep = "";
%     end
% end

if ~isfield(opts, 'sep')
    % see https://se.mathworks.com/help/matlab/ref/whitespacepattern.html
    sepList = [char(9), char(32), char(160), char(8239), char(8199)...
        , " ", "\t", ",", ";", "|", "||", "/", "//"];
    for i = 1:numel(sepList)
        lineCols = line.split(sepList(i));
        if numel(lineCols) > 1 % found delimiter
            opts.sep = sepList(i);
            break
        end
    end
else
    lineCols = line.split(opts.sep);
end

if opts.patternColAll % if 'patternCol' is nan, use all columns
    opts.headerVars = lineCols;
    opts.patternCol = 1:numel(lineCols);
    opts.patternColNumeric = true;
end

opts.headerVars = lineCols; % header variable names
if strcmp(opts.method, 'readCol')
     % keep only header names for selected columns
    if opts.patternColNumeric
        opts.headerVars = lineCols(opts.patternCol);
    else
        opts.headerVars = lineCols(ismember(lineCols, opts.patternCol));
    end
end

% check 'col'--------------------------------------------------------------
if ~opts.readAll 
    % check 'patternCol'/'filterCol' option
    opts.filterCol = unique(opts.filterCol, 'stable');
    opts.patternCol = unique(opts.patternCol, 'stable');
    opts = assertColLen(opts, 'filterCol', lineCols);
    opts = assertColLen(opts, 'patternCol', lineCols);
    
    % check for 'multiCol' option
    if opts.multiCol && (numel(opts.patternCol) ~= numel(opts.pattern))
        error('bfilereader:getbFileInfo: with ''multiCol'', ''patternCol'' and ''pattern'' must have the same size!')
    end
    
    % 'filter' and 'filterCol' must have the same length.
    if all(~isnan(opts.filter)) 
        if numel(opts.filter) ~= numel(opts.filterCol)
            error('bfilereader:getbFileInfo: ''filter'' and ''filterCol'' must have the same size!')
        elseif isscalar(opts.operator) % apply the same operator to all filter cols.
            opts.operator = repmat(opts.operator, numel(opts.filter), 1);
        elseif numel(opts.operator) ~= numel(opts.filter)
            error('bfilereader:getbFileInfo: ''filter'' and ''operator'' must have the same size!')
        end
    end
    
    [~, idx] = ismember(opts.operator, opts.dict.(1));
    opts.operator = opts.dict.(2)(idx(idx > 0));
    
    % convert col to Java integer class
    opts.patternCol = makeJavaInt(opts.patternCol, opts.method);
    opts.filterCol = makeJavaInt(opts.filterCol, opts.method);
    if ~all(isnan(opts.extractRow))
        % + 1 as it contains line numbers 
        opts.extractRow = makeJavaInt(opts.extractRow + 1, "none"); % convert it to Integer array anyway
    else
        opts.extractRow = javaArray('java.lang.Integer', 1); % read all rows
    end
    
    if ~opts.extractColAll
        opts = assertColLen(opts, 'extractCol', lineCols);
        % keep only header names for optional columns to be returned
        opts.headerVars = lineCols(opts.extractCol);
        
        % @03/02/2022 -> opts.method was replaced with "none" to allow fo
        % single column to be read
        opts.extractCol = makeJavaInt(opts.extractCol, "none"); 
    else
        opts.extractCol = javaArray('java.lang.Integer', 1); % return all columns
    end
end

% make 'skip' Java friendly
if opts.header
    opts.skip = opts.skip + 1; % don't process/parse the first line (see docs)
    lineCount = lineCount - 1; % first line
end
opts.skip = java.lang.Integer(opts.skip);

% parallel flag
if opts.parallel
    if opts.chunk > 0 % only is supported for joined pattern search
        opts.parFlag = "cp";
    else
        opts.parFlag = "sp";
    end
else
    opts.parFlag = "s";
end

if strcmp(opts.verbose, "on")
    fprintf('file has %d lines and %d columns\n', lineCount, numel(lineCols))
end
end

%% ------------------------------------------------------------------------
function opts = assertColLen(opts, name, lineCols)
if opts.(name + "Numeric")
    if max(opts.(name)) > numel(lineCols)
        fprintf('file column length:       %d\n', numel(lineCols))
        error('bfilereader:getbFileInfo:input column indices exceedes file columns!')
    end
else
    if ~all(ismember(opts.(name), lineCols))
        fprintf('file column names: %s\n', join(lineCols, "|"))
        error('bfilereader:getbFileInfo:input column names cannot be found within file header!')
    end
    opts.(name) = find(ismember(lineCols, opts.(name)));
end
end

%% ------------------------------------------------------------------------
function col = makeJavaInt(col, method)
if numel(col) > 1 || ~contains(method, 'read')
    javacol = javaArray('java.lang.Integer', numel(col));
    for i = 1:numel(col)
        javacol(i) = java.lang.Integer(col(i)-1); % MATLAB idx starts from 1.
    end
else
    javacol = java.lang.Integer(col-1);
end
col = javacol;
end

%% ------------------------------------------------------------------------
function out = javaStr2MATLAB(out, opts)

if isempty(out) || strcmp(opts.return, "raw") % nothing found?
    return
end

try
    out = reshape(split(convertCharsToStrings(out), opts.sep), numel(opts.headerVars), []).';
catch
    out = convertCharsToStrings(out);
    out = split(out, opts.sep);
    out(out == "") = [];
    out = reshape(out, numel(opts.headerVars), []).';
end

if contains(lower(opts.return), "table")
    out = array2table(out);
    if opts.header % if first row of input file contains header names
        out.Properties.VariableNames = opts.headerVars;
    else
        out.Properties.VariableNames = "Var" + (1:size(out, 2));
    end
end

if strcmp(opts.return, "table")
    % check which columns can be further converted to numeric
    emptyCells = ["", ".", "..", '""']; % ignore these
    if size(out, 1) < 200
        N = size(out, 1);
    else
        N = 200; % only check first 200 rows 
    end
    
    for i = 1:size(out, 2)
        check = out.(i)(1: N);
        check(ismember(check, emptyCells) | ismissing(check)) = [];
        if ~any(isnan(double(check))) % column may contain numeric data type
            out.(i) = double(out.(i));
        end
    end
end

end
% full END