function R2MATLABWrapper(fnc, pkg, opts)

% writes a MATLAB wrappers function using R function 'fnc' in the package
% 'pkg'.
% @19APR2023: some bugs were fixed.

arguments
    fnc {mustBeTextScalar}
    pkg {mustBeTextScalar}
    opts.funcname {mustBeTextScalar} % MATLAB wrapper function
    opts.nargout (1,1) double = 0 % output arguments
    opts.figureOutput (1,1) logical = true % R function's output is a figure. This adds input arguments to the wrapper function needed for figure output.
end

if ~isfield(opts, 'funcname')
    opts.funcname = fnc + "Wrapper";
end

fargs = getRfuncArgs(fnc, pkg);
fargs(fargs.arg == "...", :) = [];
if isempty(fargs)
    fprintf('cannot find arguments in for this func! Are you sure it''s not a class?\n')
    return
end
fargs.fis = replace(fargs.arg, ".", "_");
fargs.dot = fargs.arg ~= fargs.fis;

% write MATLAB wrapper body------------------------------------------------
inargs = fargs.fis(ismissing(fargs.val)); % function's required input arguments
fn(1) = "function ";
if opts.nargout > 0
    fn(1) = fn(1) + "[" + ("out" + (1:opts.nargout)) + "] = ";
end

if isempty(inargs)
    fn(1) = fn(1) + opts.funcname + "(opts)";
else
    fn(1) = fn(1) + opts.funcname + "(" + inargs.join(", ") + ", opts)";
end
fn(2) = "% Written by R2MATLABWrapper on " + string(datetime("today")); % descriptions

% arguments block
fn(4) = "arguments";
ct = numel(fn) + 1;
for k = ct:height(fargs)+ct-1
    if any(fargs.fis(k-ct+1) == inargs) % required arguments
        fn(k) = compose("\t%s", fargs.fis(k-ct+1));
    else % optional
        fn(k) = compose("\t%s", parseOptionalRarg(fargs(k-ct+1, :)));
    end
end

if opts.nargout > 0
    fn(end + 2) = compose("\t%s", "% non-native arguments for output");
    fn(end + 1) = compose("\t%s", 'opts.mat_out {mustBeTextScalar} = "' + fnc + '_mat" % name of output file');
end

if opts.figureOutput % add additional arguments for figure export
    fn(end + 2) = compose("\t%s", "% non-native arguments for figure export");
    fn(end + 1) = compose("\t%s", 'opts.matfig_out {mustBeTextScalar} = "' + fnc + '_fig" % name of output image');
    fn(end + 1) = compose("\t%s", "opts.matfig_res (1,1) double = 300 % image resolution");
    fn(end + 1) = compose("\t%s", "opts.matfig_width (1,1) double = 8 % image width in inch");
    fn(end + 1) = compose("\t%s", "opts.matfig_height (1,1) double = 7 % image height in inch");
end

ct = numel(fn) + 1;
fn(ct) = "end";

% prepare inputs
if opts.nargout > 0
    fn(end + 2) = "[" + join("out" + (1:opts.nargout), ",") + "] = deal([]);"; % initialize outputs
end

% check outputs (either vargout or figure)
if opts.nargout > 0
    fn(end + 2) = 'wd = fileparts(opts.mat_out);';
elseif opts.figureOutput
    fn(end + 2) = 'wd = fileparts(opts.matfig_out);';
end

if opts.nargout > 0 || opts.figureOutput
    fn(end + 1) = 'if isempty(wd) || wd == ""';
    fn(end + 1) = compose("\t%s", "wd = pwd;");

    if opts.nargout > 0
        fn(end + 1) = compose("\t%s", "opts.mat_out = fullfile(wd, opts.mat_out);");
    end

    if opts.figureOutput
        fn(end + 1) = compose("\t%s", "opts.matfig_out = fullfile(wd, opts.matfig_out);");
    end

    fn(end + 1) = "end";
else
    fn(end + 1) = "wd = pwd;";
end

if opts.figureOutput
    fn(end + 2) = 'opts.matfig_out = regexprep(opts.matfig_out, ".png$", "");';
    fn(end + 1) = 'opts.matfig_out = opts.matfig_out + ".png";';    
    fn(end + 1) = 'opts.matfig_out = replace(opts.matfig_out, filesep, "/"); % R path'; 
end

if opts.nargout > 0
    fn(end + 2) = 'opts.mat_out = replace(opts.mat_out, filesep, "/"); % R path'; 
end

fn(end + 2) = 'file = fullfile(wd, getRandomName("' + fnc + '_R", 5)); % random tmp file name';

% function's body should be filled out by user based on 'inargs'
fn(end + 2) = "% R script is defined in string 'r'; modify it as appropriate.";
fn(end + 1) = 'r(1) = "require(' + pkg + ')";';
fn(end + 1) = "% initialize and prepare required arguments to be passed to R here, e.g. r(2) = ...";

if opts.figureOutput
    fn(end + 2) = 'popts = opts;';
    fn(end + 1) = 'popts = rmfield(popts, setdiff(fieldnames(opts), ["matfig_" + ["out", "width", "height", "res"], "mat_out"]));';
    fn(end + 1) = 'fis = string(fieldnames(popts));';
    fn(end + 1) = 'for k = 1:numel(fis)';
    fn(end + 1) = compose("\t%s", 'popts.(regexprep(fis(k), "^matfig_", "")) = popts.(fis(k));');
    fn(end + 1) = "end";
    fn(end + 1) = 'popts = rmfield(popts, fis);';
    fn(end + 1) = 'r(10) = "png(" + struct2rfun(popts, skip="out") + ", units = ''in'')";';
end

fn(end + 1) = 'opts = rmfield(opts, "matfig_" + ["out", "res", "width", "height"]);';
fn(end + 1) = 'underFis = ["' + join(fargs.fis(fargs.dot) + "," + fargs.arg(fargs.dot), '", "') + '"];';

if ~isempty(inargs)
    fn(end + 1) = 'r(11) = "out = ' + string(fnc) + '(' + inargs.join(",") + ...
        ', " + struct2rfun(opts, replace=underFis) + ")";';
else
    fn(end + 1) = 'r(11) = "out = ' + string(fnc) + '(' + ...
        '" + struct2rfun(opts, replace=underFis) + ")";';
end

if ~any(fargs.dot)
    fn(end) = erase(fn(end), ", replace=underFis");
end
    
if opts.figureOutput
    fn(end + 1) = 'r(12) = "print(out)";';
    fn(end + 1) = 'r(13) = "dev.off()";';
end

fn(end + 2) = "% modify the 'r' script to handle the 'out' variable";

fn(end + 1) = 'MATLAB2Rconnector(file + ".r", code=r, delr=true);';
fn(end + 2) = "% read output(s) of R script 'r' here";
fn(end + 2) = "end % END";

opts.funcname = regexprep(opts.funcname, ".m$", "");
writelines(fn', opts.funcname + ".m")
fprintf('Wrapper function %s is ready!\n', opts.funcname + ".m")

end % END

%% subfunctions ===========================================================
function out = parseOptionalRarg(in)
% checks R function's data type and returns an equivalent MATLAB data type

out = "opts." + in.fis;

val = in.val;
if ~any(isnumeric(val)) && startsWith(in.val, "c(") % vector
    val = regexprep(val, ["^c(", ")$"], "");
    val = val.split(",");
    val = strtrim(erase(val, '"'));
    vec = true;
else
    vec = false;
end

if all(~isnan(double(val))) % numeric
    if vec, term = "double"; else, term = "(1,1) double"; end
    if vec, val = "[" + val.join(",") + "]"; end

elseif all(ismember(val, ["T", "F", "TRUE", "FALSE"])) % logical
    if vec, term = "logical"; else, term = "(1,1) logical"; end
    idx = ismember(val, ["T", "TRUE"]);
    val(idx) = "true";
    val(~idx) = "false";

    if vec
        val = "[" + val.join(",") + "]";
    end

else % string
     if all(val ~= "NULL")
         if vec
            % defVal = '"' + val(1) + '"';
            val = '["' + val.join('", "') + '"]';
         else
             val = '"' + val + '"';
         end
     else
        val = missing;
     end

    if ~vec
        term = "{mustBeTextScalar}";
    else
        term = "{mustBeMember(" + out + ", " + val + ")}";
        % val = defVal;
    end

end

if ~ismissing(val)
    out = join([out, term, "=", val], " ");
end
end