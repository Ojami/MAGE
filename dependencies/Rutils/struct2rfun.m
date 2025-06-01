function out = struct2rfun(in, opts)
% converts input struct 'in' to an R function input.
% @22MAY2023: a bug was fixed
% @05SEP2023: a bug was fixed
% @04NOV2023: now supports data.frames. If some fields of 'in' are of
% 'table' data class, will convert them to inline data.frames and appends
% to 'out' string.

arguments
    in {mustBeA(in, "struct")}
    opts.par (1,1) logical = false % enclose in parentheses
    opts.underlineTodot {mustBeText} % e.g. argument of max.size in R, and max_size in MATLAB struct. This option is deprecated and 'replace' should be used instead.
    opts.replace {mustBeText} % similar to 'underlineTodot'. Comma separated string vector of argument names to be replaced, e.g. ["old1,new1", "old2,new2"], the order doesn't matter. 
    opts.skip {mustBeText} % e.g just insert the value and not the argument, e.g. v1 instead of v1=v1
end

% find dataframe fields and parse those fields separately
istab = structfun(@(x)istable(x), in);
if any(istab)
    fis_rem = string(fieldnames(in));
    fis_rem = fis_rem(istab);
    
    df = cell(numel(fis_rem), 1);
    for k = 1:numel(fis_rem)
        df{k} = table2Rdf(in.(fis_rem(k)), name="mat_df" + k);
    end
    df = vertcat(df{:});

    for k = 1:numel(fis_rem)
        in.(fis_rem(k)) = {"mat_df" + k};
    end
end


out = string;
fis = string(fieldnames(in));

fisr = fis;
if isfield(opts, 'underlineTodot') && ~isfield(opts, 'replace')
    idx = ismember(fisr, opts.underlineTodot);
    fisr(idx) = replace(fisr(idx), "_", ".");
elseif isfield(opts, 'replace')
    opts.replace = string(opts.replace);
    if isrow(opts.replace), opts.replace = opts.replace'; end

    if isscalar(opts.replace)
        opts.replace = squeeze(split(opts.replace, ",", 2));
    else
        opts.replace = squeeze(split(opts.replace, ","));
    end
    [idx, idx2] = ismember(opts.replace, fisr);
    oldCol = any(idx, 1);
    idx = idx(:, oldCol); idx2 = idx2(:, oldCol);
    fisr(idx2(idx)) = opts.replace(idx, ~oldCol);
end

for i = 1:numel(fis)
    tmp = in.(fis(i));
    if islogical(tmp)
        tmp = string(tmp).upper;
    elseif isstring(tmp)
        tmp = replace(tmp, filesep, "/");
        if all(tmp ~= "NULL")
            tmp = '"' + tmp + '"';
        end

        if numel(tmp) > 1 % vector
            tmp = "c(" + tmp.join(",") + ")";
        end
    elseif iscellstr(tmp) || isstruct(tmp) || iscell(tmp) % df, list or other non-string arguments
        if isstruct(tmp)
            tmp = fis(i);
        else
            tmp = string(tmp);
        end
    elseif isnan(tmp)
        tmp = "NULL";        
    end

    if numel(tmp) > 1 % logical or numeric vector
        tmp = "c(" + join(string(tmp), ",") + ")";
    end
    
    if isfield(opts, 'skip') && any(fis(i) == opts.skip)
        out(i) = tmp;
    else
        out(i) = fisr(i) + "=" +  tmp;
    end

end

out = join(out, ",");
if opts.par
    out = "(" + out + ")";
end

if any(istab)
    out = [df; out];
end

end % END

