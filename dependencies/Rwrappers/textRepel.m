function res = textRepel(xx, yy, labels, opts)
% https://stackoverflow.com/questions/45065567/getting-coordinates-for-the-label-locations-from-ggrepel/45067419#45067419
% https://github.com/slowkow/ggrepel/issues/24
% Oveis Jamialahmadi, University of Gothenburg, 18 APR 2024.

arguments
    xx double {mustBeVector} 
    yy double {mustBeVector}
    labels {mustBeText, mustBeVector}
    opts.width (1,1) double = 1 % npc unit
    opts.height (1,1) double = 1 % npc unit

    opts.MarkerSize (1,1) double = 2
    opts.TextSize (1,1) double = 2
    opts.BoxPadding (1,1) double = 0.1
    opts.maxOverlap (1,1) double = 4 % increases this sequentially until no warning is found
    opts.maxIter (1,1) double = 1000
    opts.force_pull (1,1) double = 1 % Force of attraction between a text label and its corresponding data point. Defaults to 1.
end

% check if xx, yy, and labels have the same length
assert(numel(xx) == numel(yy) & numel(xx) == numel(labels), "input points and labels must have the same size!")
labels = cellstr(labels);

% R script
file = fullfile(pwd, getRandomName("textRepel", 5)); % random tmp file name
save(file + ".mat", "xx", "yy", "labels")

rcode = string;
rcode(1) = "library(grid)";
rcode(2) = "library(ggrepel)";
rcode(3) = "library(ggplot2)";
rcode(4) = "library(tidyr)";
rcode(5) = "library(dplyr)";

rcode(7) = "get_repel_coords <- function(.data, g_base, width=" + opts.width + ", height=" + opts.height + ", ...) {";
% rcode(numel(rcode)+1) = "    tmp_file <- tempfile(fileext = '.png')";
% rcode(numel(rcode)+1) = "    png(tmp_file)";
rcode(numel(rcode)+1) = "    grid.newpage()";
rcode(numel(rcode)+1) = "    pushViewport(viewport(width = width, height = height))";
rcode(numel(rcode)+1) = '    g <- g_base + geom_text_repel(aes(x, y), label = ".", data = .data, max.overlaps = Inf, ...)';
rcode(numel(rcode)+1) = "    panel_params <- ggplot_build(g)$layout$panel_params[[1]]";
rcode(numel(rcode)+1) = "    xrg <- panel_params$x.range";
rcode(numel(rcode)+1) = "    yrg <- panel_params$y.range";

rcode(numel(rcode)+1) = '    textrepeltree <- ggplotGrob(g) %>% grid.force(draw = F) %>% getGrob("textrepeltree", grep = T)';
rcode(numel(rcode)+1) = '    children <- childNames(textrepeltree) %>% grep("textrepelgrob", ., value = T)';

rcode(numel(rcode)+1) = "    get_xy <- function(n) {";
rcode(numel(rcode)+1) = "        grob <- getGrob(textrepeltree, n)";
rcode(numel(rcode)+1) = '        data.frame(x.repel = xrg[1] + diff(xrg) * convertX(grob$x, "native", valueOnly = T), y.repel = yrg[1] + diff(yrg) * convertY(grob$y, "native", valueOnly = T))';
rcode(numel(rcode)+1) = "    }";
rcode(numel(rcode)+1) = "    lapply(children, get_xy) %>% bind_rows %>% cbind(.data)";
% rcode(numel(rcode)+1) = "    dev.off()  # Close the PNG device";
% rcode(numel(rcode)+1) = "    file.remove(tmp_file)  # Remove the temporary file";
% rcode(numel(rcode)+1) = "    return(res)";
rcode(numel(rcode)+1) = "}";

% read data
rcode(numel(rcode)+1) = "tmp = R.matlab::readMat('" + replace(file + ".mat", filesep, "/") + "')";
rcode(numel(rcode)+1) = "df <- data.frame(x = tmp$xx, y = tmp$yy, z = unlist(tmp$labels))";
% rcode(numel(rcode)+1) = "pdf(NULL)";
rcode(numel(rcode)+1) = "mo <- " + opts.maxOverlap;
rcode(numel(rcode)+1) = "custom_env <- new.env()";
rcode(numel(rcode)+1) = 'assign("last.warning", NULL, envir = custom_env)';

rcode(numel(rcode)+1) = "repeat {";
rcode(numel(rcode)+1) = "  tryCatch({";
rcode(numel(rcode)+1) = "    p = ggplot(df, aes(x = x, y = y, label = z)) + " + ...
    "geom_point(size = " + opts.MarkerSize + ") + " + ...
     "geom_text_repel(size = " + opts.TextSize + ...
     ", box.padding = unit(" + opts.BoxPadding + ", 'line'), " + ...
     " max.overlaps = mo, max.iter = " + ...
     opts.maxIter + ", force_pull = + " + opts.force_pull + ") + theme_minimal() +  theme(legend.position = 'none')";
  
rcode(numel(rcode)+1) = "    print(p)";
rcode(numel(rcode)+1) = "    break";
rcode(numel(rcode)+1) = "  }, warning = function(w) {";
rcode(numel(rcode)+1) = '    if (grepl("Consider increasing max.overlaps", w$message)){';
rcode(numel(rcode)+1) = "      mo <<- mo + 1}";
rcode(numel(rcode)+1) = '     assign("last.warning", w, envir = custom_env)';
rcode(numel(rcode)+1) = "  })";
rcode(numel(rcode)+1) = "}";

rcode(numel(rcode)+1) = "res = get_repel_coords(df, p)";
rcode(numel(rcode)+1) = "write.table(res, '" + replace(file, filesep, "/") + ".txt', row.names = F, col.names = T, sep = '||', quote = F)";
rcode(numel(rcode)+1) = "dev.off()";

MATLAB2Rconnector(file + ".r", code=rcode);
res = readtable(file + ".txt", "TextType", "string", ...
    "VariableNamingRule", "preserve", "NumHeaderLines", 0, ...
    "Delimiter", "||");
res = renamevars(res, ["x.repel", "y.repel"], ["x_repel", "y_repel"]);
delete(file + ".mat")
delete(file + ".txt")
if isfile("Rplots.pdf"), delete("Rplots.pdf"); end

end % END