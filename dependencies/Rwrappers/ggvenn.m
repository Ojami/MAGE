function ggvenn(sets, opts)
% a wrapper around ggVennDiagram R pacakge
% Oveis Jamialahmadi, 6 Feb 2023, University of Gothenburg.

arguments
    sets {mustBeA(sets, 'struct')} % struct with each field corresponding to a different set
    opts.label {mustBeMember(opts.label, ["both", "count", "percent", "none"])} = "count"
    opts.label_alpha (1,1) double = 0
    opts.label_geom {mustBeMember(opts.label_geom, ["label", "text"])} = "text"
    opts.label_color {mustBeTextScalar} = "Black"
    opts.label_size (1,1) double
    opts.set_size (1,1) double
    opts.category_names {mustBeText} % default is names(sets) in R
    opts.x_expand (1,1) double = 0.2 % to expand x-axis when category_names are long

    opts.out {mustBeTextScalar} = "ggvenn" % name of output image
    opts.res (1,1) double = 300 % image resolution
    opts.width (1,1) double = 8 % image width in inch
    opts.height (1,1) double = 7 % image height in inch
end

if opts.label_geom == "text"
    opts.label_alpha = 0.5;
end

wd = fileparts(opts.out);
if isempty(wd) || wd == ""
    wd = pwd;
    opts.out = fullfile(wd, opts.out);
end

opts.out = regexprep(opts.out, ".png$", "");
opts.out = opts.out + ".png";

opts.out = replace(opts.out, filesep, "/"); % R path 
file = fullfile(wd, getRandomName("ggvenn", 5)); % random tmp file name

% convert to cellstr all fields
fis = string(fieldnames(sets));
for k = 1:numel(fis), sets.(fis(k)) = cellstr(string(sets.(fis(k)))); end
save(file + ".mat", "-struct", "sets")

r(1) = "lapply(c('ggVennDiagram', 'ggplot2'), require, character.only = TRUE)";
r(2) = "genes = R.matlab::readMat('" + replace(file, filesep, "/") + ".mat')";
r(3) = "genes = lapply(genes, unlist)";
r(4) = "png('" + opts.out + "', width = " + opts.width + ...
    ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
rmfis = ["out", "res", "width", "height"];
opts = rmfield(opts, rmfis);
r(5) = "p1 = ggVennDiagram::ggVennDiagram(genes, " + ...
    struct2rfun(opts, underlineTodot="category_names") + ") + " + ...
    "scale_x_continuous(expand = expansion(mult = " + opts.x_expand + "))";
r(6) = "print(p1)";
r(7) = "dev.off()";
MATLAB2Rconnector(file + ".r", code=r, delr=true);
delete(file + ".mat")

% r(10) = "ggplot() + geom_sf(aes(fill = id), data = venn_region(data)," + ...
%     "show.legend = FALSE) +  geom_sf(color='grey', size = 3, data = " + ...
%     "venn_setedge(data), show.legend = FALSE) + " + ...
%     "geom_sf_text(aes(label = name), fontface = 'bold', " + ...
%     "data = venn_setlabel(data)) + geom_sf_text(aes(label = name)," + ...
%     "data = venn_region(data)) + theme_void()";

end % END