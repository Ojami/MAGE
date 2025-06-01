function filterline(index, file, opts)
% a wrapper for https://github.com/miku/filterline

arguments
    index {mustBeTextScalar, mustBeFile} % line numbers (1-based) per each file
    file {mustBeTextScalar, mustBeFile} % file that filtering should be performed on
    opts.keep (1,1) logical = false % keep the lines (true) or remove (false, default)?
    opts.output {mustBeTextScalar} % if left empty inputfile + "_filtered" suffix will be added
    opts.verbose (1,1) logical = true
end

% locate the binary file
exe = string(regexprep(which("filterline.m"), ".m$", "")); 
assert(isfile(exe), "cannot find filterline binary file!")
exe = "'" + makeWSLpath(exe) + "'";

% check the keep flag
if opts.keep
    opts.keep = "";
else
    opts.keep = " -v ";
end

% get full path
[~, index] = fileattrib(index);
index = string(index.Name);

[~, file] = fileattrib(file);
file = string(file.Name);

if ~isfield(opts, "output")
    [pth, name, ext] = fileparts(file);
    opts.output = fullfile(pth, name + "_filtered" + ext);
else
    [pth, name, ext] = fileparts(opts.output);
    if string(pth) == "", pth = pwd; end
    opts.output = fullfile(pth, name + ext);
end

index = "'" + makeWSLpath(index) + "'";
file = "'" + makeWSLpath(file) + "'";
opts.output = "'" + makeWSLpath(opts.output) + "'";

cmd = string;
cmd(1) = "sort -n " + index + " -o " + index;
if opts.output == file % in-place
    cmd(2) = exe + opts.keep + " " + index + " " + file + " > " + ...
        "temp.txt && mv temp.txt " + opts.output;
else
    cmd(2) = exe + opts.keep + " " + index + " " + file + " > " + opts.output;
end

runbash(cmd, "parallel", false, "verbose", opts.verbose);

end % END
