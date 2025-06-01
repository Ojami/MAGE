function EnhancedVolcano(tab, opts)

% UNDER DEVELOPMENT: not all options have been wrapped.
% a wrapper around EnhancedVolcano R package.
% Oveis Jamialahmadi, University of Gothenburg, January 2023.

arguments
    tab {mustBeA(tab, "table")} % a table of DE genes
    opts.out {mustBeTextScalar} = "EnhancedVolcano" % name of output image
    opts.res (1,1) double = 300 % image resolution
    opts.width (1,1) double = 8 % image width in inch
    opts.height (1,1) double = 7 % image height in inch
    opts.Nlabs (1,1) double = 15 % max number of labels to show
    
    opts.lab {mustBeTextScalar} = "gene_name" % column with labels
    opts.FCcutoff (1,1) double = 1
    opts.pCutoff (1,1) double = 0.05
    opts.x {mustBeTextScalar} = "log2FoldChange"
    opts.y {mustBeTextScalar} = "pvalue"
    opts.title {mustBeTextScalar} = "NULL"
    opts.subtitle {mustBeTextScalar} = "NULL"
    opts.caption {mustBeTextScalar} = "NULL"
    opts.colAlpha (1,1) double = 0.8
    opts.pointSize (1,1) double = 3.0
    opts.boxedLabels (1,1) logical = false
    opts.drawConnectors (1,1) logical = true
    opts.max_overlaps (1,1) double = 50

    opts.selectLab {mustBeText, mustBeVector} % selected labels. If not provided, will be selected automatically based on FC and P.
end

wd = fileparts(opts.out);
if isempty(wd) || wd == ""
    wd = pwd;
    opts.out = fullfile(wd, opts.out);
end

opts.out = regexprep(opts.out, ".png$", "");
opts.out = opts.out + ".png";

opts.out = replace(opts.out, filesep, "/"); % R path 
file = fullfile(wd, getRandomName("EnhancedVolcano", 5)); % random tmp file name

idx = ismissing(tab.(opts.lab));
tab.(opts.lab)(idx) = " ";

% identify top genes ------------------------------------------------------
idxmin = tab.(opts.y) < realmin;
tab.(opts.y)(idxmin) = realmin;
idx = abs(tab.(opts.x)) >= opts.FCcutoff & tab.(opts.y) <= opts.pCutoff;
tab = tab(:, [opts.x, opts.y, opts.lab]);
labs = tab(idx, :);
labs.metric = abs(labs.(opts.x)).*-log10(labs.(opts.y));
labs = sortrows(labs, "metric", "descend");
opts.Nlabs = min(height(labs), opts.Nlabs);

if ~isfield(opts, "selectLab")
    opts.selectLab = labs.(opts.lab)(1:opts.Nlabs);
end
    
opts.selectLab(opts.selectLab == " ") = [];

% prepare input table -----------------------------------------------------
labs = cellstr(tab.(opts.lab));
slabs = cellstr(opts.selectLab);
tab = table2array(tab(:, [opts.x, opts.y]));
cols = cellstr([opts.x, opts.y, opts.lab]);
save(file + ".mat", "tab", "labs", "slabs", "cols")
opts.selectLab = 'slabs';

rcode(1) = "tmp = R.matlab::readMat('" + replace(file + ".mat", filesep, "/") + "')";
rcode(2) = "cols = unlist(tmp$cols, use.names=F)";
rcode(3) = "slabs = unlist(tmp$slabs, use.names=F)";
rcode(4) = "lab = unlist(tmp$labs, use.names=F)";
rcode(5) = "tmp = data.frame(tmp$tab, lab)";
rcode(6) = "colnames(tmp) = cols";
rcode(9) = "png('" + opts.out + "', width = " + opts.width + ...
    ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";

rmfis = ["lab", "Nlabs", "out", "res", "width", "height"];
opts = rmfield(opts, rmfis);
rcode(8) =  "evplot=EnhancedVolcano::EnhancedVolcano(tmp, lab=lab, " + ...
    struct2rfun(opts, underlineTodot="max_overlaps") + ")";
rcode(10) = "print(evplot)";
rcode(11) = "dev.off()";
MATLAB2Rconnector(file + ".r", code=rcode, delr=true);
delete(file + ".mat")

end % END

