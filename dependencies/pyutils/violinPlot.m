function violinPlot(tab, opts)
% A simple seaborn.violinplot wrapper:
% https://seaborn.pydata.org/generated/seaborn.violinplot.html
% Oveis Jamialahmadi, University of Gothenburg, 14JULY2023.

arguments
    tab {mustBeA(tab, "table")}
    opts.x {mustBeTextScalar}
    opts.y {mustBeTextScalar}
    opts.hue {mustBeTextScalar} % grouping variable
    opts.split (1,1) logical = false % to split in half
    opts.dodge (1,1) logical = true
    opts.orient {mustBeMember(opts.orient, ["v", "h"])} = "v"
    opts.cut (1,1) double = 0
    opts.saturation (1,1) double = 1
    opts.inner {mustBeMember(opts.inner, ["box", "quartile", "point", "stick", "None"])} = "box"

    % bot plox: only when inner == "box"
    opts.box_width (1,1) double = 0.3
    opts.box_marker {mustBeTextScalar} = "o";
    opts.box_markersize (1,1) double = 3
    opts.box_alpha (1,1) double = 0.5 % marker alpha

    % theme
    opts.palette {mustBeMember(opts.palette, ["Accent", "Accent_r", ...
        "Blues", "Blues_r", "BrBG", "BrBG_r", "BuGn", "BuGn_r", "BuPu", "BuPu_r", ...
        "CMRmap", "CMRmap_r", "Dark2", "Dark2_r", "GnBu", "GnBu_r", "Greens", "Greens_r", "Greys", "Greys_r", "OrRd", ...
         "OrRd_r", "Oranges", "Oranges_r", "PRGn", "PRGn_r", "Paired", "Paired_r", "Pastel1", ...
         "Pastel1_r", "Pastel2", "Pastel2_r", "PiYG", "PiYG_r", "PuBu", "PuBuGn", "PuBuGn_r", ...
         "PuBu_r", "PuOr", "PuOr_r", "PuRd", "PuRd_r", "Purples", "Purples_r", "RdBu", "RdBu_r", ...
         "RdGy", "RdGy_r", "RdPu", "RdPu_r", "RdYlBu", "RdYlBu_r", "RdYlGn", "RdYlGn_r", "Reds", ...
         "Reds_r", "Set1", "Set1_r", "Set2", "Set2_r", "Set3", "Set3_r", "Spectral", "Spectral_r", ...
         "Wistia", "Wistia_r", "YlGn", "YlGnBu", "YlGnBu_r", "YlGn_r", "YlOrBr", "YlOrBr_r", "YlOrRd",... 
         "YlOrRd_r", "afmhot", "afmhot_r", "autumn", "autumn_r", "binary", "binary_r", "bone", ...
         "bone_r", "brg", "brg_r", "bwr", "bwr_r", "cividis", "cividis_r", ...
         "cool", "cool_r", "coolwarm", "coolwarm_r", "copper", "copper_r",...
         "cubehelix", "cubehelix_r", "flag", "flag_r", "gist_earth",...
         "gist_earth_r", "gist_gray", "gist_gray_r", "gist_heat", "gist_heat_r", "gist_ncar", "gist_ncar_r",...
         "gist_rainbow", "gist_rainbow_r", "gist_stern", "gist_stern_r", "gist_yarg", ...
         "gist_yarg_r", "gnuplot", "gnuplot2", "gnuplot2_r", "gnuplot_r", "gray", "gray_r",...
         "hot", "hot_r", "hsv", "hsv_r", "icefire", "icefire_r", "inferno", ...
         "inferno_r", "magma", "magma_r", "mako", "mako_r", ...
         "nipy_spectral", "nipy_spectral_r", "ocean", "ocean_r", "pink", "pink_r",...
         "plasma", "plasma_r", "prism", "prism_r", "rainbow", "rainbow_r",...
         "rocket", "rocket_r", "seismic", "seismic_r", "spring", "spring_r",...
         "summer", "summer_r", "tab10", "tab10_r", "tab20", "tab20_r", "tab20b",...
         "tab20b_r", "tab20c", "tab20c_r", "terrain", "terrain_r", "twilight",...
         "twilight_r", "twilight_shifted", "twilight_shifted_r", "viridis",...
         "viridis_r", "vlag", "vlag_r", "winter", "winter_r", "None"])} = "None"
    opts.style {mustBeMember(opts.style, ["darkgrid", "whitegrid", "dark", "white", "ticks"])} = "darkgrid"
    opts.font {mustBeTextScalar} = "Garamond"
    opts.n_colors (1,1) double % number of colors for palette
    
    opts.alpha (1,1) double = 1 % transparency
    opts.x_fontsize (1,1) double = 12
    opts.x_label {mustBeTextScalar}
    opts.y_fontsize (1,1) double = 12
    opts.y_label {mustBeTextScalar}
    opts.tick_fontsize (1,1) double = 9
    opts.legend_fontsize (1,1) double = 9
    opts.title_fontsize (1,1) double = 17
    opts.linewidth (1,1) double 
    opts.title {mustBeTextScalar}    
    opts.output {mustBeTextScalar} = "violin_plot"
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "svg", "pdf"])} = "png" % format of output plots
    opts.dpi (1, 1) double = 400 % resolution of output plot
end

% output
opts.output = string(opts.output);
pth = fileparts(opts.output);
if pth == "", pth = pwd; end
if ~isfolder(pth), mkdir(pth); end

if ~isfield(opts, "x_label")
    opts.x_label = opts.x;
end

if ~isfield(opts, "y_label")
    opts.y_label = opts.y;
end

% terminate(pyenv)
df = table2pandasDF(tab);
p = pyenv;
pyTCL = fullfile(p.Home, 'tcl', 'tcl8.6');
pyTK = fullfile(p.Home, 'tcl', 'tk8.6');
setenv('TCL_LIBRARY', pyTCL);
setenv('TK_LIBRARY', pyTK);
sns = py.importlib.import_module("seaborn");
plt = py.importlib.import_module("matplotlib.pyplot");
font_manager = py.importlib.import_module("matplotlib.font_manager");

% close all figures
plt.close("all");
plt.figure();
sns.reset_defaults();

font_path = pyrun("font_path = font_manager.findfont(font_manager.FontProperties(family=font_name))", ...
    "font_path", font_name=opts.font, font_manager=font_manager);
font_prop = pyrun("font_prop = font_manager.FontProperties(fname=font_path)",...
    "font_prop", font_manager=font_manager, font_path=font_path);

% set theme
cmd = getPycmd("sns.set_theme", opts, fi=["palette", "style", "font"]);
pyrun(cmd, sns=sns);
cmd = getPycmd("sns.color_palette", opts, out="pt", fi=["palette", "n_colors"]);
pt = pyrun(cmd, "pt", sns=sns, data=df);
opts.palette = pt;

fi = ["x", "y", "hue", "split", "dodge", "orient", "cut", "linewidth", ...
    "saturation", "inner", "data", "palette"];
inopts = opts;
inopts.data = df;
if opts.inner == "box", inopts.inner = "None"; end
cmd = getPycmd("sns.violinplot", inopts, out="ax", fi=fi, in=["data","palette"]);
ax = pyrun(cmd, "ax", sns=sns, data=df, palette=pt);

% set alpha transparency
pyrun("for patch in ax.collections: patch.set_alpha(" + opts.alpha + ")", ax=ax);

% draw box plots 
if opts.inner == "box"
    inopts.ax = "ax";
    inopts.width = opts.box_width;
    inopts.boxprops = {"{'zorder': 2}"};
    inopts.flierprops= {"{'marker': '" + opts.box_marker + ...
        "', 'markersize': " + opts.box_markersize + ", 'alpha':" +...
        opts.box_alpha + "}"};
    fi = setdiff(fi, ["inner", "cut", "split"]);
    fi = union(fi, ["width", "ax", "boxprops", "flierprops"]);
    cmd = getPycmd("sns.boxplot", inopts, out="ax", fi=fi, in=["ax", "data", "palette"]);
    ax = pyrun(cmd, "ax", sns=sns, data=df, palette=pt);
end

% set font size and name (not necessary now since we use set_theme)
if isfield(opts, "title")
    % ax.axes.set_title(opts.title, fontsize=opts.title_fontsize, fontproperties=font_prop);
    ax.set_title(opts.title, fontproperties=font_prop, fontsize=opts.title_fontsize);
end
ax.set_xlabel(opts.x_label, fontproperties=font_prop, fontsize=opts.x_fontsize);
ax.set_ylabel(opts.y_label, fontproperties=font_prop, fontsize=opts.y_fontsize);
ax.tick_params(labelsize=opts.tick_fontsize);
plt.setp(ax.get_xticklabels(), fontname=opts.font);
plt.setp(ax.get_yticklabels(), fontname=opts.font);

if isfield(opts, "hue")
    if opts.inner == "box"
        outvars = ax.get_legend_handles_labels();
        pyrun("ax.legend(handles[:2], labels[:2], title='" + opts.hue + "')", ax=ax, labels=outvars{2}, handles=outvars{1});
    end

    sns.move_legend(ax, "upper left", bbox_to_anchor=py.tuple([1, 1]), ...
        fontsize=opts.legend_fontsize, title_fontsize=opts.legend_fontsize+1);
    % plt.legend(title=opts.hue, fontsize=opts.legend_fontsize, ...
    %     loc="best", bbox_to_anchor=py.tuple([1.05, 1]));
    
end

plt.tight_layout();
plt.savefig(opts.output.replace("\", "/") + "." + opts.format, dpi=opts.dpi)

end % END