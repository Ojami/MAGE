function cmd = getPycmd(func, str, opts)
% contructs a python command to be run using pyrun.
% Oveis Jamialahmadi, University of Gothenburg, 15JULY2023

arguments
    func {mustBeTextScalar} % python func
    str {mustBeA(str, "struct")} % struct of arguments
    opts.fi {mustBeText, mustBeVector} % use only these fields
    opts.out {mustBeTextScalar} % output variable
    opts.in {mustBeText, mustBeVector} % a subset of 'fi' or str fields that are inputs in 'pyrun'. In this case, field name is used as field value (e.g. df=df)
end

cmd = func + "(";

if isfield(opts, "out")
    cmd = opts.out + " = " + cmd;
end

if ~isfield(opts, "fi")
    fi = string(fieldnames(str));
else
    fi = string(opts.fi);
end

if isfield(opts, "in")
    % modify input fields and use their field name instead of value
    infi = intersect(fi, opts.in);
    for k = 1:numel(infi)
        str.(infi(k)) = string(infi(k));
    end
else
    opts.in = "";
end

for k = 1:numel(fi)
    if isfield(str, fi(k))
        if isstring(str.(fi(k))) && ~ismember(fi(k), opts.in)
            if str.(fi(k)).lower == "none"
                cmd = cmd + fi(k) + "=" + str.(fi(k)) + ",";
            else
                cmd = cmd + fi(k) + "='" + str.(fi(k)) + "',";
            end
        
        elseif iscell(str.(fi(k)))
            cmd = cmd + fi(k) + "=" + str.(fi(k)) + ",";

        else
            if islogical(str.(fi(k)))
                str.(fi(k)) = regexprep(string(str.(fi(k))), '^\w{1}', '${upper($0)}');
            end
            cmd = cmd + fi(k) + "=" + str.(fi(k)) + ",";

        end
    end
end

cmd = regexprep(cmd, ",$", "") + ")";

end % END