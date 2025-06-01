function upsetWrapper(data, opts)
% Written by R2MATLABWrapper on 07-Feb-2023

arguments
	data {mustBeA(data, 'struct')}
	opts.nsets (1,1) double = 5
	opts.nintersects (1,1) double = 40
	opts.sets
	opts.keep_order (1,1) logical = false
	opts.set_metadata
	opts.intersections
	opts.matrix_color {mustBeTextScalar} = "gray23"
	opts.main_bar_color {mustBeTextScalar} = "gray23"
	opts.mainbar_y_label {mustBeTextScalar} = "Intersection Size"
	opts.mainbar_y_max
	opts.sets_bar_color {mustBeTextScalar} = "gray23"
	opts.sets_x_label {mustBeTextScalar} = "Set Size"
	opts.point_size (1,1) double = 2.2
	opts.line_size (1,1) double = 0.7
	opts.mb_ratio double = [0.7,0.3]
	opts.expression
	opts.att_pos
	opts.att_color {mustBeTextScalar} = "main.bar.color"
	opts.order_by {mustBeMember(opts.order_by, ["freq", "degree"])} = ["freq", "degree"]
	opts.decreasing logical = [true,false]
	opts.show_numbers {mustBeTextScalar} = "yes"
	opts.number_angles (1,1) double = 0
	opts.group_by {mustBeMember(opts.group_by, ["sets", "degree"])} = "degree"
	opts.cutoff
	opts.queries
	opts.query_legend {mustBeTextScalar} = "none"
	opts.shade_color {mustBeTextScalar} = "gray88"
	opts.shade_alpha (1,1) double = 0.25
	opts.matrix_dot_alpha (1,1) double = 0.5
	opts.empty_intersections
	opts.color_pal (1,1) double = 1
	opts.boxplot_summary
	opts.attribute_plots
	opts.scale_intersections {mustBeTextScalar} = "identity"
	opts.scale_sets {mustBeTextScalar} = "identity"
	opts.text_scale (1,1) double = 1
	opts.set_size_angles (1,1) double = 0
	opts.set_size_show (1,1) logical = false
	opts.set_size_numbers_size
	opts.set_size_scale_max

	% non-native arguments for figure export
	opts.matfig_out {mustBeTextScalar} = "upset_fig" % name of output image
	opts.matfig_res (1,1) double = 300 % image resolution
	opts.matfig_width (1,1) double = 8 % image width in inch
	opts.matfig_height (1,1) double = 7 % image height in inch
    opts.format {mustBeMember(opts.format, ["png", "pdf"])} = "png"
end

wd = fileparts(opts.matfig_out);
if isempty(wd) || wd == ""
	wd = pwd;
	opts.matfig_out = fullfile(wd, opts.matfig_out);
end

opts.matfig_out = regexprep(opts.matfig_out, opts.format + ".$", "");
opts.matfig_out = opts.matfig_out + "." + opts.format;
opts.matfig_out = replace(opts.matfig_out, filesep, "/"); % R path

file = fullfile(wd, getRandomName("upset_R", 5)); % random tmp file name

% convert to cellstr all fields
fis = string(fieldnames(data));
for k = 1:numel(fis), data.(fis(k)) = cellstr(string(data.(fis(k)))); end
save(file + ".mat", "-struct", "data")

% R script is defined in string 'r'; modify it as appropriate.
r(1) = "require(UpSetR)";
r(2) = "data = R.matlab::readMat('" + replace(file, filesep, "/") + ".mat')";
r(3) = "data = lapply(data, unlist)";

popts = opts;
popts = rmfield(popts, setdiff(fieldnames(opts), "matfig_" + ["out", "width", "height", "res"]));
fis = string(fieldnames(popts));
for k = 1:numel(fis)
	popts.(regexprep(fis(k), "^matfig_", "")) = popts.(fis(k));
end
popts = rmfield(popts, fis);
if opts.format ~= "png"
    popts = rmfield(popts, setdiff(fieldnames(popts), ["out", "width", "height"]));
    r(10) = opts.format + "(" + struct2rfun(popts, skip="out") + ", onefile=FALSE)";
else
    r(10) = "png(" + struct2rfun(popts, skip="out") + ", units = 'in')";
end
opts = rmfield(opts, ["format", "matfig_" + ["out", "res", "width", "height"]]);
underFis = ["keep_order,keep.order", "set_metadata,set.metadata", "matrix_color,matrix.color", "main_bar_color,main.bar.color", "mainbar_y_label,mainbar.y.label", "mainbar_y_max,mainbar.y.max", "sets_bar_color,sets.bar.color", "sets_x_label,sets.x.label", "point_size,point.size", "line_size,line.size", "mb_ratio,mb.ratio", "att_pos,att.pos", "att_color,att.color", "order_by,order.by", "show_numbers,show.numbers", "number_angles,number.angles", "group_by,group.by", "query_legend,query.legend", "shade_color,shade.color", "shade_alpha,shade.alpha", "matrix_dot_alpha,matrix.dot.alpha", "empty_intersections,empty.intersections", "color_pal,color.pal", "boxplot_summary,boxplot.summary", "attribute_plots,attribute.plots", "scale_intersections,scale.intersections", "scale_sets,scale.sets", "text_scale,text.scale", "set_size_angles,set_size.angles", "set_size_show,set_size.show", "set_size_numbers_size,set_size.numbers_size", "set_size_scale_max,set_size.scale_max"];
r(11) = "out = upset(fromList(data), " + struct2rfun(opts, replace=underFis) + ")";
r(12) = "print(out)";
r(13) = "dev.off()";

% modify the 'r' script to handle the 'out' variable
MATLAB2Rconnector(file + ".r", code=r, delr=true);

delete(file + ".mat")

end % END
