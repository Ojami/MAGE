function pheatmap(tab, opts)

% UNDER DEVELOPMENT: not all options have been wrapped.
% generates heatmap using R pheatmap package.
% Oveis Jamialahmadi, University of Gothenburg, December 2022.
% 
% @19JUNE2023: more arguments were added.
% @05SEP2023: general struct2rfun is now used instead of internal one.


arguments
    tab double % a double matrix
    opts.out {mustBeTextScalar} = "pheatmap" % name of output image
    opts.res (1,1) double = 300 % image resolution
    opts.width (1,1) double = 7 % image width in inch
    opts.height (1,1) double = 8 % image height in inch

    opts.color {mustBeMember(opts.color, ["Blues", "BuGn", "BuPu", "GnBu", ...
        "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", ...
        "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd",...
        "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", ...
        "RdYlGn", "Spectral", "Accent", "Dark2", "Paired", "Pastel1", ...
        "Pastel2", "Set1", "Set2", "Set3", "magma", "inferno", "plasma", ...
        "viridis", "cividis", "rocket", "mako", "turbo"])}
    opts.ncolor (1,1) double = 100
    opts.scale {mustBeMember(opts.scale, ["none", "row", "column"])} = "none"
    opts.rownames {mustBeText, mustBeVector} % row names 
    opts.colnames {mustBeText, mustBeVector} % column names
    opts.annotation_row {mustBeA(opts.annotation_row, "table")} % data frame that specifies the annotations shown on left side of the heatmap
    opts.annotation_col {mustBeA(opts.annotation_col, "table")} % similar to annotation_row, but for columns.
    opts.cluster_rows (1,1) logical = true
    opts.cluster_cols (1,1) logical = true
    opts.cutree_rows (1,1) double % number of clusters the rows are divided into, based on the hierarchical clustering (using cutree)
    opts.cutree_cols (1,1) double
    opts.legend (1,1) logical = true
    opts.fontname {mustBeTextScalar} = "Garamond"
    opts.format {mustBeMember(opts.format, ["png", "svg", "pdf"])} = "png"
    opts.bold (1,1) logical = false
    opts.na_col {mustBeTextScalar} = "white"
    opts.display_numbers (1,1) = false % logical determining if the numeric values are also printed to the cells. 
    opts.number_format {mustBeTextScalar} = "%.2f"
    opts.number_color {mustBeTextScalar}
    opts.cellwidth (1,1) double 
    opts.cellheight (1,1) double
    opts.gaps_row double {mustBeVector}
    opts.gaps_col double {mustBeVector}
    opts.show_rownames (1,1) logical
    opts.show_colnames (1,1) logical
    opts.fontsize (1,1) double 
    opts.fontsize_row (1,1) double % ontsize for rownames (Default: fontsize)
    opts.fontsize_col (1,1) double % fontsize for colnames (Default: fontsize)
end

if isfield(opts, "color")
    if any(opts.color == ["magma", "inferno", "plasma", "viridis", ...
            "cividis", "rocket", "mako", "turbo"])
        opts.color = {"viridis::" + opts.color + "(" + opts.ncolor + ")"};
    else
        opts.color = {"RColorBrewer::brewer.pal(" + opts.ncolor + ", '" + opts.color + "')"};
    end
    opts = rmfield(opts, "ncolor");
end

fontname = opts.fontname;
opts = rmfield(opts, 'fontname');

wd = fileparts(opts.out);
if isempty(wd) || wd == ""
    wd = pwd;
    opts.out = fullfile(wd, opts.out);
end

opts.out = regexprep(opts.out, "."+opts.format+"$", "");
opts.out = opts.out + "."+opts.format;

opts.out = replace(opts.out, filesep, "/"); % R path 
file = fullfile(wd, getRandomName("pheatmap", 5)); % random tmp file name

% prepare input table -----------------------------------------------------
if isfield(opts, 'colnames')
    if numel(opts.colnames) ~= size(tab, 2)
        error("pheatmap:size mismatch between column names and input table!")
    end
    opts.colnames = cellstr(opts.colnames);
    if ~isfield(opts, "show_colnames")
        opts.show_colnames = true;
    end
else
    opts.colnames = {'NA'};
    opts.show_colnames = false;
end

if isfield(opts, 'rownames')
    if numel(opts.rownames) ~= size(tab, 1)
        error("pheatmap:size mismatch between row names and input table!")
    end
    opts.rownames = cellstr(opts.rownames);
    if ~isfield(opts, "show_rownames")
        opts.show_rownames = true;
    end
else
    opts.show_rownames = false;
    opts.rownames = {'NA'};
end

rownames = opts.rownames;
colnames = opts.colnames;
opts = rmfield(opts, ["rownames", "colnames"]);
save(file + ".mat", "tab", "rownames", "colnames")

rcode(1) = "tmp = R.matlab::readMat('" + replace(file + ".mat", filesep, "/") + "')";
rcode(2) = "tab = tmp$tab";
rcode(3) = "if (length(setdiff(unlist(tmp$rownames), 'NA')) > 0){";
rcode(4) = "rownames(tab) = unlist(tmp$rownames)}";
rcode(5) = "if (length(setdiff(unlist(tmp$colnames), 'NA')) > 0){";
rcode(6) = "colnames(tab) = unlist(tmp$colnames)}";

rcode(7) = "library(tidyverse)";
rcode(8) = "make_bold_names <- function(mat, rc_fun, rc_names=NULL) {";
rcode(9) = "  bold_names <- rc_fun(mat)";
rcode(10) = "  if (!is.null(rc_names))";
rcode(11) = "    ids <- rc_names %>% match(rc_fun(mat))";
rcode(12) = "  else";
rcode(13) = "    ids = 1:length(rc_fun(mat))";
rcode(14) = "  ids %>%";
rcode(15) = "    walk(";
rcode(16) = "      function(i)";
rcode(17) = "        bold_names[i] <<-";
rcode(18) = "        bquote(bold(.(rc_fun(mat)[i]))) %>%";
rcode(19) = "        as.expression()";
rcode(20) = "    )";
rcode(21) = "  bold_names}";

if opts.display_numbers
    rcode(23) = "idx = is.na(tab) | is.nan(tab)";
    rcode(24) = "disp_mat = matrix(sprintf('" + opts.number_format + "', tab), nrow = nrow(tab))";
    rcode(25) = 'disp_mat[idx] = ""';
    opts.display_numbers = {'disp_mat'};
end

if opts.format == "svg"
    rcode(30) = "svg('" + opts.out + "', width = " + opts.width + ...
        ", height = " + opts.height + ")";
elseif opts.format == "pdf"
    rcode(30) = "pdf('" + opts.out + "', width = " + opts.width + ...
        ", height = " + opts.height + ")";
else
    rcode(30) = opts.format + "('" + opts.out + "', width = " + opts.width + ...
        ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
end

optsin = opts;
optsin = rmfield(optsin, ["height", "res", "width", "out"]);

rtmp = struct2rfun(optsin);
rtmp(end) =  "pheatmap::pheatmap(tab, " + rtmp(end);
rcode = [rcode'; rtmp];
n = numel(rcode);
if opts.bold
    rcode(n) = rcode(n) + ", labels_row = make_bold_names(tab, rownames), labels_col = make_bold_names(tab, colnames))";
else
    rcode(n) = rcode(n) + ")";
end
rcode(n + 1) = "dev.off()";

rcode = addRfont(rcode, fontFamily=fontname);
MATLAB2Rconnector(file + ".r", code=rcode);
delete(file + ".mat")

end % END

% %% subfunctions ===========================================================
% function out = struct2rfun(in, opts)
% % converts input struct 'in' to an R function input.
% arguments
%     in {mustBeA(in, "struct")}
%     opts.par (1,1) logical = false % enclose in parentheses
% end
% 
% out = string;
% fis = string(fieldnames(in));
% for i = 1:numel(fis)
%     tmp = in.(fis(i));
%     if islogical(tmp)
%         tmp = string(tmp).upper;
%     elseif isstring(tmp)
%         tmp = replace(tmp, filesep, "/");
%         tmp = '"' + tmp + '"';
%     elseif iscellstr(tmp) % df or other non-string arguments
%         tmp = string(tmp);
%     end
%     out(i) = fis(i) + "=" +  tmp;
% end
% 
% out = join(out, ",");
% if opts.par
%     out = "(" + out + ")";
% end
% 
% end % END