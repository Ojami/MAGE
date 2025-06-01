function res = tabix(gzfile, region, opts)
% tabix is part of hstlib (currently version 1.14)
% to compile and make see htslib page
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, Jan 2022.
% 
% @25/01/2022: 'offline' option was added: "force" to download files before
%               query fetch. "false" to use online data. "lazy" to download
%               only if data already exist locally (in 'cachedir')
%               'downloadiffailed' flag was added: true, to download data
%               if online fetching was failed; false, otherwise (error)
% @07/09/2022: some bugs were fixed.
% @25JAN2023: some bugs were fixed.

arguments
    gzfile {mustBeTextScalar}
    region {mustBeText, mustBeVector}
    opts.reader {mustBeMember(opts.reader, ["readtable", "bfilereader", "best"])} = "readtable" % best: if gcp isn't empty, uses tall, otherwise uses bfilereader
    opts.downloadindex (1,1) logical = false % -D flag
    opts.offline {mustBeMember(opts.offline, ["force", "false", "lazy"])} = "force" % lazy: use offline data only if they already axist 
    opts.downloadiffailed (1,1) logical = true % download files if server suppresses remote access
    opts.verbose (1,1) logical = true
    opts.cachedir {mustBeTextScalar}
end

if ispc
    wslpref = "wsl ";
else
    wslpref = "";
end

% check if tabix has been installed
[~, msg] = system(wslpref + "/bin/bash -ilc 'which tabix'");
if isempty(msg)
    error('tabix cannot be found. See: http://www.htslib.org/download/\n')
end

name = getRandomName("tabix", 5);

if ~isfield(opts, 'cachedir')
    opts.cachedir = fullfile(pwd, "tabix.cache");
end

if ~isfolder(opts.cachedir) && (opts.downloadiffailed || opts.offline == "force")
    mkdir(opts.cachedir)
end

if opts.offline == "lazy" % check existence of files in cachedir
    tmpgzfile = regexprep(gzfile, "^ftp://", "https://");
    [~, savename, saveext] = fileparts(tmpgzfile);
    if isfile(fullfile(opts.cachedir, savename + saveext))
        opts.offline = "force";
    end
end

[hdr, res] = deal([]);

if ~isempty(hdr) || opts.offline == "force"

    hdr = splitlines(string(hdr));
    hdr(contains(hdr, "file not found)")) = [];
    hdr(hdr == "") = [];

    if any(contains(hdr, "onnection reset by peer")) || opts.offline == "force"
        if opts.downloadiffailed || opts.offline == "force" % download locally
            gzfile = regexprep(gzfile, "^ftp://", "https://");
            [~, savename, saveext] = fileparts(gzfile);
            if ~isfile(fullfile(opts.cachedir, savename + saveext))
                websave(fullfile(opts.cachedir, savename + saveext), gzfile);
                websave(fullfile(opts.cachedir, savename + saveext + ".tbi"), gzfile + ".tbi");
            end
            gzfile = fullfile(opts.cachedir, savename + saveext);
            hdr = bfilereader(gzfile, 'summary', 'firstline');
            gzfile = makeWSLpath(gzfile);

        else
            error('%s\n', hdr{:})
        end
    end
end


if isempty(hdr)
    rcode(1, 1) = "enum <- FALSE";
    rcode(2, 1) = "for (i in 1:70){";
    rcode(3, 1) = "tryCatch( { hdr <- colnames(readr::read_tsv('" + gzfile + "', n_max = 1))}";
    rcode(4, 1) = "       , error = function(e) {enum <<- TRUE})";
    rcode(5, 1) = "if(!enum){break}}";
    rcode(6, 1) = "write(hdr, file = '" + name + "')";
    writematrix(rcode, name + ".r", 'FileType', 'text', 'QuoteStrings', false)
    MATLAB2Rconnector(name + ".r", 'delr', true);
    try
        hdr = readmatrix(name, 'FileType', 'text', 'OutputType', 'string');
    catch
        fprintf('tabix:failed fetch hdr from %s\n', name + ".r")
        if isfile(name + ".r"), delete(name + ".r"); end 
    end
    delete(name)
end

if isempty(hdr), return, end

if numel(region) > 1 % multiple regions
    region = squeeze(split(region, [":", "-", "_"]));
    writematrix(region, name + ".bed", 'QuoteStrings', false,...
        'FileType', 'text', 'Delimiter', '\t')
    multiFlag = true;
else
    multiFlag = false;
end

% @29MAY2023: -D flag has been removed.
% if ~opts.downloadindex
%     cmd = wslpref + "/bin/bash -ilc ""tabix -D ";
% else
    cmd = wslpref + "/bin/bash -ilc ""tabix ";
% end

if multiFlag
    [~, ~] = system(cmd + gzfile + " -R " + name + ".bed" + " > " + name + '"');
    delete(name + ".bed")
else
    [~, ~] = system(cmd + gzfile + " " + region + " > " + name + '"');
end


% check if empty
if dir(name).bytes < 1
    res = []; 
    delete(name)
    return
end

if opts.reader ~= "bfilereader" && ~isempty(gcp('nocreate'))
    [~, ~, ext] = fileparts(name);

    % silent gather: https://se.mathworks.com/matlabcentral/answers/511958-getting-a-command-like-gather-to-run-silently
    matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
    res = gather(tall(tabularTextDatastore(name, 'TextType' , 'string', 'FileExtensions', ext)));
    matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);
elseif opts.reader == "bfilereader"
    res = bfilereader(name, 'parallel', true, 'verbose', 'off');
else
    res = readtable(name, 'FileType', 'text', 'TextType', 'string');
end
delete(name)
if numel(hdr) ~= width(res)
    disp(hdr)
    fprintf('hdr size: %d - table size: %d\n', numel(hdr), width(res))
    error('tabix:column names mismatch table width!')
end
res.Properties.VariableNames = hdr;

end 