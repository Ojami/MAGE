function res = Rpam(data, opts)
% Written by R2MATLABWrapper on 19-Mar-2024
% @19MAR2024: uses PAM (kmediods in MATLAB) in R to find optimal number of
% clusters (if 'k' is nan:default), and after clustering, evaluates
% different measures using fpc and clValid R packages:
% https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/
%   
%   Internal validation measures:
%   1-Silhouette coefficient: for observations: ~1 best, ~0 between
%   clusters, <0 wrong cluster
%   2-Dunn index: If the data set contains compact and well-separated
%   clusters, the diameter of the clusters is expected to be small and the
%   distance between the clusters is expected to be large. Thus, Dunn index
%   should be maximized.
% 
%   External validation measures: External clustering validation, can be
%   used to select suitable clustering algorithm for a given data set.
%   1-Rand index and Meila2s VI: provides a measure for assessing the
%   similarity between two partitions, adjusted for chance varies from -1
%   (no agreement) to 1 (perfect agreement).
% 
%   Stability measures: testing how sensitive it is to perturbations in the
%   input data.
%   1-In clValid package this means removing each column one at
%   a time and re-rnning the clustering. There are several measures
%   included, such as average proportion of non-overlap (APN), the average
%   distance (AD),the average distance between means (ADM), and the figure
%   of merit (FOM), all of which should be minimised.
%   2-Stability of clusters can be also computed using fpc::clusterboot()
%   function, however the perturbations in the input data is a bit
%   different here: it’s done by resampling the data in a chosen way (e.g.
%   bootstraping).
%       -Clusterwise Jaccard bootstrap mean should be maximised
%       -number of dissolved clusters should be minimised
%       -number of recovered clusters should be maximised and as close to
%       the number of pre-defined bootstraps as possible
% 
% Notes: 
%   1-Note that the stability of a cluster is assessed, but stability is
%   not the only important validity criterion - clusters obtained by very
%   inflexible clustering methods may be stable but not valid, as discussed
%   in Hennig (2007). https://search.r-project.org/CRAN/refmans/fpc/html/clusterboot.html
%   

arguments
	data {mustBeA(data, "table")}
	
	opts.diss {mustBeTextScalar} % = "inherits(x, "dist")"
	opts.metric {mustBeMember(opts.metric, ["euclidean", "manhattan"])} = "euclidean"
	opts.medoids {mustBeTextScalar} % = "if (is.numeric(nstart)) "random""
	opts.nstart {mustBeTextScalar} % = "if (variant == "faster") 1 else NA"
	opts.stand (1,1) logical = false
	opts.cluster_only (1,1) logical = false
	opts.do_swap (1,1) logical = true
	opts.keep_diss {mustBeTextScalar} % = "!diss && !cluster.only && n < 100"
	opts.keep_data {mustBeTextScalar} % = "!diss && !cluster.only"
	opts.variant {mustBeMember(opts.variant, ["original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster"])}
	opts.pamonce (1,1) logical = false
	opts.trace_lev (1,1) double = 0

	% non-native arguments
	opts.out {mustBeTextScalar} = "Rpam" % name of output files
	opts.res (1,1) double = 300 % image resolution
	opts.width (1,1) double = 8 % image width in inch
	opts.height (1,1) double = 7 % image height in inch

    opts.k (1,1) double = nan % number of clusters
    opts.seed (1,1) double = 123 % for reproducibility
    opts.nboot (1,1) double = 2e3 % for fpc:clusterboot
end

wd = fileparts(opts.out);
if isempty(wd) || wd == ""
	wd = pwd;
	opts.out = fullfile(wd, opts.out);
end
if ~isfolder(wd), mkdir(wd); end

opts.out = regexprep(opts.out, [".pdf$", ".png$"], "");
opts.out = replace(opts.out, filesep, "/"); % R path

file = fullfile(wd, getRandomName("pam_R", 5)); % random tmp file name

% R script is defined in string 'r'; modify it as appropriate.
r(1) = "pckg = c(""cluster"", ""factoextra"", ""fpc"", ""clValid"")";
r(2) = "lapply(pckg, function(x) suppressMessages(require(x, character.only = TRUE)))";

% get row and col names
data = rmmissing(data);
rownames = data.Row;
if isempty(rownames)
    rownames = "y" + (1:height(data));
end
rnames = cellstr(rownames);
cnames = cellstr(colnames(data));
data = normalize(data{:, :}); % normalize (same as scale in R)
save(file + ".mat", "data", "rnames", "cnames")

r(3) = "tmp = R.matlab::readMat('" + replace(file + ".mat", filesep, "/") + "')";
r(4) = "data = tmp$data";
r(5) = "rownames(data) = unlist(tmp$rnames)";
r(6) = "colnames(data) = unlist(tmp$cnames)";
r(7) = "set.seed(" + opts.seed + ", kind = 'Mersenne-Twister')";

% check number of clusters
if isnan(opts.k)
    r(8) = "sfig = fviz_nbclust(data, pam, method ='silhouette')+theme_minimal()"; % silhouette plot
    r(9) = "pdf('" + opts.out + "_silhouette.pdf', width = " + opts.width + ...
        ", height = " + opts.height + ")";
    r(10) = "print(sfig)";
    r(11) = "dev.off()";

    r(12) = "png('" + opts.out + "_silhouette.png', width = " + opts.width + ...
        ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
    r(13) = "print(sfig)";
    r(14) = "dev.off()";

    r(15) = "nc = pamk(data)$nc";
else
    r(15) = "nc = " + opts.k;
end

optsin = opts;
optsin = rmfield(optsin,  ["out", "width", "height", "res", "k", "seed", "nboot"]);
underFis = ["cluster_only,cluster.only", "do_swap,do.swap", "keep_diss,keep.diss", "keep_data,keep.data", "trace_lev,trace.lev"];
r(17) = "out = pam(data,nc, " + struct2rfun(optsin, replace=underFis) + ")";
r(18) = "write.table(out$clustering, file='" + opts.out + "_cluster.txt', quote = F, sep = '||', col.names = F, row.names = T)";

% visualize
r(19) = "cfig = fviz_cluster(out,ellipse.type ='convex',repel =TRUE, ggtheme =theme_bw(),main='')"; % silhouette plot
r(20) = "pdf('" + opts.out + "_cluster.pdf', width = " + opts.width + ...
    ", height = " + opts.height + ")";
r(21) = "print(cfig)";
r(22) = "dev.off()";

r(23) = "png('" + opts.out + "_cluster.png', width = " + opts.width + ...
    ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
r(24) = "print(cfig)";
r(25) = "dev.off()";

% stability measures
r(26) = "ss = fpc::clusterboot(data, B=" + opts.nboot + ...
    ", bootmethod = 'boot', clustermethod = claraCBI, " + ...
    "k=nc, metric='" + opts.metric + ...
    "', seed=set.seed(" + opts.seed + ", kind = 'Mersenne-Twister'))";
r(27) = "write.table(t(ss$bootresult), '" + opts.out + "_stability.txt', row.names = F, col.names = F, sep='||')";

% validity measures
r(28) = "valid_test <- clValid::clValid(data, union(nc, 2:10), " + ...
    "clMethods = c('kmeans', 'pam', 'hierarchical'), validation = c('internal', 'stability'))";
r(29) = "ctab1 = as.data.frame(valid_test@measures)";
r(30) = "ctab1$Measure = rownames(ctab1)";
r(31) = "write.table(ctab1, '" + opts.out + "_validity1.txt', quote = F, row.names = F, col.names = T, sep = '||')";
r(32) = "ctab2 = optimalScores(valid_test)";
r(33) = "ctab2$Measure = rownames(ctab2)";
r(34) = "write.table(ctab2, '" + opts.out + "_validity2.txt', quote = F, row.names = F, col.names = T, sep = '||')";

% modify the 'r' script to handle the 'out' variable
MATLAB2Rconnector(file + ".r", code=r, delr=true);

% read output(s) of R script 'r' here
delete(file + ".mat")

% read clusters
ctab = readtable(opts.out + "_cluster.txt", "VariableNamingRule", "preserve", ...
    "TextType", "string", "Delimiter", "||");
ctab.Properties.VariableNames = ["var", "cluster"];

% read stability info (Jaccard index from bootstrapping step)
stab = readtable(opts.out + "_stability.txt", "VariableNamingRule", "preserve", ...
    "TextType", "string", "Delimiter", "||");
stab.Properties.VariableNames = "cluster" + (1:width(stab));
stab_mean = varfun(@mean, stab);
stab_sd = varfun(@std, stab);
stab_mean = stab_mean{:, :}; stab_sd = stab_sd{:, :};
overall_jaccard = compose("%.2f (%.2f)", mean(stab_mean), std(stab_mean));
cluster_jaccard = compose("%.2f (%.2f)", stab_mean', stab_sd');

% overall validity
vtab1 = readtable(opts.out + "_validity1.txt", "VariableNamingRule", "preserve", ...
    "TextType", "string", "Delimiter", "||");
vtab1 = movevars(vtab1, "Measure", Before=1);

% optimal validity
vtab2 = readtable(opts.out + "_validity2.txt", "VariableNamingRule", "preserve", ...
    "TextType", "string", "Delimiter", "||");

delete(opts.out + "_stability.txt")
delete(opts.out + "_validity1.txt")
delete(opts.out + "_validity2.txt")
delete(opts.out + "_cluster.txt")

res = struct;
res.cluster = ctab;
res.stability.boot = stab;
res.stability.overall_jaccard = overall_jaccard;
res.stability.cluster_jaccard = cluster_jaccard;
res.validity.overall = vtab1;
res.validity.optimal = vtab2;

end % END
