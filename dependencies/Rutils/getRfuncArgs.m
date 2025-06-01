function [args, fnc] = getRfuncArgs(fnc, pkg, opts)
% read input arguments with their default values from an R function in a
% package.

arguments
    fnc {mustBeTextScalar} % function name
    pkg {mustBeTextScalar} % package name
    opts.readfunc (1,1) logical = false % reads function content
end

r(1) = "out = formals(" + pkg + "::" + fnc + ")";
r(2) = "outn = names(out)";
r(3) = "out = as.character(out)";
r(4) = "names(out) = outn";
r(5) = "write.table(out, 'RfncArgs_4matlab.txt', col.names = F, sep='||', quote = F)";

if opts.readfunc
    r(6) = "fncContent = deparse(" + pkg + ":::" + fnc + ")";
    r(7) = "write(fncContent, 'RfncContent_4matlab.txt')";
end

MATLAB2Rconnector("getRfuncArgs.r", code=r, delr=true);
args = readlines("RfncArgs_4matlab.txt");
args(args == "") = [];
try
    args = array2table(args.split("||"), VariableNames=["arg", "val"]);
    args.val(args.val == "") = missing;
catch % 
end
delete("RfncArgs_4matlab.txt")

if opts.readfunc
    fnc = readlines("RfncContent_4matlab.txt");
    delete("RfncContent_4matlab.txt")
else
    fnc = [];
end

end % END