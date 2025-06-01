function rgb_matrix = get_colormap(name, n)

% uses seaborn package to get colormaps not implemented in MATLAB.
% Oveis Jamialahmadi, University of Gothenburg, 09 SEPT 2024. 

arguments
    name {mustBeMember(name, ["magma", "inferno", "plasma", "viridis", ...
        "cividis", "twilight", "twilight_shifted", "turbo", "Blues", ...
        "BrBG", "BuGn", "BuPu", "CMRmap", "GnBu", "Greens", "Greys", ...
        "OrRd", "Oranges", "PRGn", "PiYG", "PuBu", "PuBuGn", "PuOr", ...
        "PuRd", "Purples", "RdBu", "RdGy", "RdPu", "RdYlBu", "RdYlGn", ...
        "Reds", "Spectral", "Wistia", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", ...
        "afmhot", "autumn", "binary", "bone", "brg", "bwr", "cool", ...
        "coolwarm", "copper", "cubehelix", "flag", "gist_earth", ...
        "gist_gray", "gist_heat", "gist_ncar", "gist_rainbow", ...
        "gist_stern", "gist_yarg", "gnuplot", "gnuplot2", "gray", "hot", ...
        "hsv", "jet", "nipy_spectral", "ocean", "pink", "prism", "rainbow", ...
        "seismic", "spring", "summer", "terrain", "winter", "Accent", "Dark2", ...
        "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3", "tab10", ...
        "tab20", "tab20b", "tab20c", "rocket", "mako", "icefire", "vlag", ...
        "flare", "crest", "colorblind", "deep", "muted"])} = "magma"
    n = 1
end

% % % list of names:
% plt = py.importlib.import_module("matplotlib.pyplot");
% colormap_names = plt.colormaps;
% colormap_names = colormap_names.dictionary.keys;
% colormap_names(colormap_names.endsWith("_r")) = []; % delete reverse

warning off
% plt = py.importlib.import_module("matplotlib.pyplot");
% try
%     colormap = plt.get_cmap(name);
% catch
% 
%     colormap = py.seaborn.color_palette(name, as_cmap=true);
% end
% values = py.numpy.linspace(0, 1, py.int(n));
% rgb_matrix = colormap(values);
% rgb_matrix = double(rgb_matrix);
% rgb_matrix = rgb_matrix(:, 1:3);

try
    % Try to get a matplotlib colormap (callable)
    plt = py.importlib.import_module("matplotlib.pyplot");
    colormap = plt.get_cmap(name);
    values = py.numpy.linspace(0, 1, py.int(n));
    rgb_matrix = colormap(values);
catch
    % If that fails, use seaborn's color_palette.
    % First, make sure to import seaborn
    seaborn = py.importlib.import_module("seaborn");
    % Specify the number of colors using the n_colors parameter.
    colormap_list = seaborn.color_palette(name, pyargs('n_colors', py.int(n)));
    % Convert the list to a numpy array.
    rgb_matrix = py.numpy.array(colormap_list);
end

% Convert the Python array (or list) to a MATLAB numeric matrix (dropping alpha if present).
rgb_matrix = double(rgb_matrix);

if size(rgb_matrix,2) > 3
    rgb_matrix = rgb_matrix(:, 1:3);
end

warning on

% colormap = py.seaborn.color_palette(name, as_cmap=true);
% 
% % Generate RGB values for the specified number of data points
% rgb_values = pyrun("rgb_values = colormap(range(num_colors))", "rgb_values", colormap=colormap, num_colors=py.int(n));
% 
% % Convert to a matrix (list of RGB tuples)
% rgb_matrix = pyrun("rgb_matrix = [list(color) for color in rgb_values]", "rgb_matrix", rgb_values=rgb_values);
% rgb_matrix = cell(rgb_matrix);
% rgb_matrix = cellfun(@(x)double(x), rgb_matrix, uni=false);
% rgb_matrix = vertcat(rgb_matrix{:});
% rgb_matrix = rgb_matrix(:, 1:3).;

end