function df = table2pandasDF(tab)
% converts MATLAB table to Pandas data frame
% Oveis Jamialahmadi, University of Gothenburg, 14JULY2023.

arguments
    tab {mustBeA(tab, "table")}
end

pd = py.importlib.import_module("pandas");
cols = colnames(tab);

df = struct;
cat_cols = false(width(tab), 1);
for k = 1:width(tab)
    df.("x" + k) = tab.(k)';
    if isstring(df.("x" + k))
        df.("x" + k) = cellstr(df.("x" + k));
    elseif islogical(df.("x" + k)(~ismissing(df.("x" + k))))
        df.("x" + k) = py.numpy.array(df.("x" + k));
        df.("x" + k) = py.numpy.ravel(df.("x" + k));
    elseif iscategorical(df.("x" + k))
        df.("x" + k) = cellstr(string(df.("x" + k)));
        cat_cols(k) = true;
    end
end
df = py.dict(df);
df = pd.DataFrame(df);
df.columns = cellstr(cols);

if any(cat_cols)
    cat_cols = cellstr(cols(cat_cols));
    for k = 1:numel(cat_cols)
        cat_type = pd.api.types.CategoricalDtype(categories=categories(tab.(cat_cols{k}))', ordered=true);
        df = pyrun("df[cat_cols] = df[cat_cols].astype(cat_type)", "df", df=df, cat_cols=py.list(cat_cols(k)), cat_type=cat_type);
    end
end

end % END