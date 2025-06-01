function t = fig2tile(files, opts)
% fig2tile copies axes from input fig files (MATLAB figures with fig
% suffix) to a single tiledlayout.
% Oveis Jamialahmadi, University of Gothenburg, August 2022.
% 
% TODO: 'Position' and 'TightInset' should be used to fill the horizontal
% space betwee tiles when input figures are 'square'.

arguments
    files {mustBeFile, mustBeVector}
    opts.GridSize (1, 2) double % [m, n] grid size of tiledlayout, by default this is not used, and instead 'flow' decides the grid size.
    opts.optimGridSize (1,1) logical = false % optimizes the grid size and overrides 'GridSize' option. It uses a ga optmization with nonlinear constaints (UNDER CONSTRUCTION!)
    opts.title {mustBeText, mustBeVector} % titles to be shown on each tile. Will replace the existing titles
    opts.titleFontSize (1,1) double % font size of titles in 'title' option. If left empty, will be axes font size + 2
    opts.TileSpacing {mustBeMember(opts.TileSpacing, ["compact", "none", "loose", "tight"])} = "compact"
    opts.Padding {mustBeMember(opts.Padding, ["compact", "loose", "tight"])} = "loose"
    opts.WindowState {mustBeTextScalar, mustBeMember(opts.WindowState, ["maximized", "normal"])} = "maximized"
    opts.output {mustBeTextScalar} = "fig2tile_output"
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "tif" % format of output plots
    opts.resolution (1, 1) double = 600 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    opts.legLocation {mustBeMember(opts.legLocation, ["north", "south", ...
        "east", "west", "northeast", "northwest", "southeast", "southwest", ...
        "northoutside", "southoutside", "eastoutside", "westoutside", ...
        "northeastoutside", "northwestoutside", "southeastoutside", ...
        "southwestoutside", "best", "bestoutside", "none", "layout"])} % if there exists a legend 
end

if any(~endsWith(files, ".fig"))
    error('fig2tile:all inputs must be MATLAB figures!')
end

if isfield(opts, 'title')
    opts.title = string(opts.title);
    if numel(opts.title) ~= numel(files)
        error("fig2title:'title' size should match number of figures!")
    end
end

close all force
f = figure("WindowState", opts.WindowState);
t = tiledlayout(f, 'flow', "TileSpacing", opts.TileSpacing);
if isfield(opts, 'GridSize')
    if opts.optimGridSize
        opts.GridSize = optimizeGrids(numel(files));
    end
    t.GridSize = opts.GridSize;
end

for i = 1:numel(files)
    setTitle = true;
    fig = openfig(files(i), "invisible");
    hax = copyobj(fig.Children, t);
    for j = 1:numel(hax)
        % check legend
        if endsWith(class(hax(j)), ".Legend")
            if isfield(opts, 'legLocation')
                hax_loc = opts.legLocation;
            else
                hax_loc = hax(j).Location;
            end
        else
            hax_loc = "";
        end
        hax(j).Layout.Tile = i;

        if hax_loc ~= ""
            drawnow;
            hax(j).Location = hax_loc;
        end

        if isfield(opts, 'title')
            if (~strcmp(hax(j).Title.String, '') || j == numel(hax)) && setTitle
                hax(j).Title.String = opts.title(i);
                if ~isfield(opts, 'titleFontSize')
                    opts.titleFontSize = hax(j).FontSize + 2;
                end
                hax(j).Title.FontSize = opts.titleFontSize;
                setTitle = false;
            end
            
        end
    end
%     set(fig.Children, "Parent", t)
    close(fig)
end

if opts.save
    exportgraphics(t.Parent, opts.output + "." + opts.format, 'Resolution', opts.resolution)
end

if opts.savefig % doesn't work with "scatter"
    savefig(t.Parent, opts.output + ".fig", 'compact')
end

end % END

%% subfunctions ===========================================================
function gsize = optimizeGrids(n)
opts = optimoptions('ga', 'Display', 'off');
opts.ConstraintTolerance=1e-9;
opts.FunctionTolerance=1e-9;
% 
A = [1 -1];
b = 1;
gsize = ga(@(x)(x(1)+x(2)),2,A,n,[],[],[1,1],[],@(x)nonLinConst(x,n),opts);
gsize = ceil(round(gsize, 2, 'significant'));
end

function [c,ceq] = nonLinConst(x, a)

c = a - x(1)*x(2);
ceq = [];
end
