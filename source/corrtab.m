function [rho, p] = corrtab(tab, opts)
% correlation for table input
% Oveis Jamialahmadi, June 2023.

arguments
    tab {mustBeA(tab, "table")}
    opts.sigOnly (1,1) logical = false % only return significant (alpha = 0.05) correlations
    opts.plot (1,1) logical = false % heatmap of correlation coefficients
    opts.name {mustBeTextScalar} = "corrTab"
    opts.height (1,1) double = 9 % for plot
    opts.width (1,1) double = 9 % for plot

    % corr arguments
    opts.Type {mustBeMember(opts.Type, ["Pearson", "Kendall", "Spearman"])} = "Pearson"
    opts.Rows {mustBeMember(opts.Rows, ["all", "complete", "pairwise"])} = "all"
    opts.Tail {mustBeMember(opts.Tail, ["both", "left", "right"])} = "both"
end

doubleIdx = cellfun(@(x)string(x) == "double", varfun(@class, tab, OutputFormat="cell"));
tab = tab(:, doubleIdx);

optsin = opts;
optsin = rmfield(optsin, ["plot", "sigOnly", "name", "height", "width"]);
optsin = namedargs2cell(optsin);
[rho, p] = corr(tab{:, :}, optsin{:});

if opts.sigOnly
    rho(p > 0.05) = nan;
end
rho = array2table(rho, VariableNames=colnames(tab), RowNames=colnames(tab));
p = array2table(p, VariableNames=colnames(tab), RowNames=colnames(tab));

if opts.plot
    pheatmap(rho{:, :}, "colnames", colnames(tab), "rownames", colnames(tab), ...
    "cluster_rows", ~opts.sigOnly, "cluster_cols", ~opts.sigOnly, "out", ...
    opts.name, "width", opts.width, "height", opts.height, "display_numbers", true)
end

end % END