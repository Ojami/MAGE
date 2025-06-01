function [df, fname] = table2Rdf(tab, opts)

% converts MATLAB table to inline (within R script) definition of an R
% data.frame. Oveis Jamialahmadi, University of Gothenbug, 11 Nov 2023.
% 
% @19MAY2024: now supports writing data in the table to a temp file and
% read it into R. Usefull for big tables that cannot be defined inline.
arguments
    tab {mustBeA(tab, "table")}
    opts.name {mustBeTextScalar} = "df"; % name of data.frame
    opts.n (1,1) double = 20 % to break each string and create a new line (to avoid R wrapping errors)
    
    opts.write (1,1) logical = false
    opts.dt {mustBeMember(opts.dt, ["parquet", "txt"])} = "parquet" % output file type
    opts.dir {mustBeTextScalar} = pwd % dir to write/read the df to/from
end

opts.name = string(opts.name);
cols = colnames(tab);

% check row names
if isempty(tab.Row)
    norow = true;
else
    norow = false;
    row_names = string(tab.Row);
end

fname = missing;

% 19MAY2024: writes the table to a file to be read into R
if opts.write
    if ~isfolder(opts.dir), mkdir(opts.dir); end
    fname = getRandomName("df");

    if opts.dt == "parquet"
        fname = fullfile(opts.dir, fname + ".parquet");
        parquetwrite(fname, tab)
    else
        fname = fullfile(opts.dir, fname + ".txt");
        if isempty(gcp('nocreate'))
            writetable(tab, fname, "Delimiter", "\t", "WriteRowNames", false)
        else
            fastWriteTable(tab, output=fname, delimiter="\t")
        end
    end

    df = string;
    df(1) = "library(dplyr)";
    df(2) = "tmp_file = '" + fname.replace(filesep, "/") + "'";

    if opts.dt == "parquet"
        df(3) = opts.name + "= as.data.frame(arrow::read_parquet(tmp_file))";
    else
        df(3) = opts.name + "= data.table::fread(tmp_file, sep='\t')";
        df(numel(df) + 1) = "names(" + opts.name + ') = c("' + join(cols, '","') + '")';
        df(numel(df) + 1) = "if (file.exists(tmp_file)) {";
        df(numel(df) + 1) = "file.remove(tmp_file)}";
        df(numel(df) + 1) = opts.name + " <- " + opts.name + " %>% mutate_all(~ifelse(is.nan(.), NA, .))";
    end

    % convert to appropriate variable types
    ct = numel(df) + 2;
    for k = 1:width(tab)
        if islogical(tab.(k))
            df(ct) = opts.name + " <- " + opts.name + " %>% mutate_at" + ...
                "('" + cols(k) + "', as.logical)";
            ct = ct + 1;
        elseif iscategorical(tab.(k))
            if isordinal(tab.(k))
                ordered = "T";
            else
                ordered = "F";
            end
    
            cats = string(categories(tab.(k)));
            cats = 'c("' + join(cats, '", "') + '")';
            df(ct) = opts.name + "$" + cols(k) + "=" + ...
                    "factor(" + opts.name + "$" + cols(k) + ", order=" + ordered + ", levels=" + cats + ")";
            ct = ct + 1;
        end
    end

    df = df';
    return
end

df_out = cell(width(tab), 1);
for k = 1:width(tab)

    v = tab.(k);
    if isstring(v) || iscellstr(v)
        if iscell(v)
            v(ismissing(v)) = {missing};
            v = string(v);
        end
        
        v(ismissing(v)) = "NA";
        v = string(v);
        v(v == "") = "NA";
        v = createRvec(v, opts.n, true); % 'c("' + join(v, '", "') + '")';
        v = replace(v, '"NA"', "NA");
    elseif islogical(v(~ismissing(v)))
        v = createRvec(upper(string(v)), opts.n, false);
    elseif iscategorical(v)
        if isordinal(v)
            ordered = "T";
        else
            ordered = "F";
        end

        cats = string(categories(v));
        v = string(v);
        v(ismissing(v) | v == "") = "NA";
        v = createRvec(v, opts.n, true); % 'c("' + join(v, '", "') + '")';
        v = replace(v, '"NA"', "NA");

        cats = 'c("' + join(cats, '", "') + '")';
        cats = replace(cats, '"NA"', "NA");
        
        if numel(v) > 1
            v = ["factor(" + v(1); v(2:end); ", order=" + ordered + ", levels=" + cats + ")"];
        else
            v = "factor(" + v + ", order=" + ordered + ", levels=" + cats + ")";
        end
       
    elseif isnumeric(v)
        v = string(v);
        v(ismissing(v)) = "NA";
        v = createRvec(v, opts.n, false);
    else
        error("unknown data class!")
    end
    
    if k == 1
        if numel(v) > 1
            df = [opts.name + " = data.frame(a" + k + "=" + v(1); v(2:end)];
        else
            df = opts.name + " = data.frame(a" + k + "=" + v;
        end
    else
        if numel(v) > 1
            df = [", a" + k + "=" + v(1); v(2:end)];
        else
            df = ", a" + k + "=" + v;
        end
    end

    if k == width(tab)
        if ~norow
            rnames = createRvec(row_names, opts.n, true);
            if numel(rnames) > 1
                df_out{k} = [df; ', row.names = ' + rnames(1); rnames(2:end-1); rnames(end) + ")"]; %c("' + join(row_names, '","') + '")';
            else
                df_out{k} = [df; ', row.names = ' + rnames + ")"];
            end
        else
            df_out{k} = [df; ")"];
        end
        
    else
        df_out{k} = df;
    end

end

df_out = vertcat(df_out{:});
df = [df_out; "names(" + opts.name + ') = c("' + join(cols, '","') + '")'];

end % END

