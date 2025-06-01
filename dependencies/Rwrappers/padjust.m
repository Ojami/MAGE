function padj = padjust(p, opts)
% A wrapper around R p.adjust function.
% @16MAY2023

arguments
	p {mustBeVector, mustBeNonempty, mustBeNonNan}
	opts.method {mustBeTextScalar, mustBeMember(opts.method, ["holm", ...
        "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" ])} = "BH"
	opts.n {mustBeTextScalar} 
end

if ~isfield(opts, "n"), opts.n = numel(p); end
file = getRandomName("p.adjust_R", 5); % random tmp file name

save(file + ".mat", "p")
r(1) = "p = R.matlab::readMat('" + file + ".mat" + "')$p";
r(2) = "out = stats::p.adjust(p, method='" + opts.method + "', n=" + opts.n + ")";
r(3) = "data.table::fwrite(as.matrix(out), '" + file + ".txt', row.names = F, " + ...
    "col.names = F, quote = F)";

MATLAB2Rconnector(file + ".r", code=r, delr=true);
padj = readmatrix(file + ".txt");

delete(file + ".mat")
delete(file + ".txt")

end % END
