function effDotPlotter(tab, opts)
% dot plot of effects for a set of loci and traits:
%         T1 T2 T3
%  locus1 O--O--O
%  locus2 O--O--O
%  locus3 O--O--O
% Oveis Jamialahmadi, University of Gothenburg, 22 Jan 2024.

arguments
    tab {mustBeA(tab, "table")} % table of summary stats with the following cols:
    opts.pCol {mustBeTextScalar} = "P" % p-value of the association
    opts.xCol {mustBeTextScalar} = "Trait" % shown on the x-axis
    opts.yCol {mustBeTextScalar} = "Locus" % shown on the y-axis
    opts.bCol {mustBeTextScalar} = "Beta" % effect column.
    opts.bLabel {mustBeTextScalar} = "Effect" % Colorbar label of bCol variable

    opts.szRange (1,2) double = [30, 110]
    opts.WindowState {mustBeTextScalar, mustBeMember(opts.WindowState, ...
        ["maximized", "normal"])} = "maximized" 
    opts.WindowPosition (1,4) double = nan(4, 1) % windows position in pixels, overrides WindowState
    opts.visible {mustBeTextScalar, mustBeMember(opts.visible, ["on", "off"])} = "on" % use off only for debugging/rendering issues
    opts.pbaspect (1,3) double = [1, 4, 1] % plot box aspect ratio
    

    % options for lines
    % constant horizontal lines connecting dots
    opts.linecolor {mustBeTextScalar, mustBeMember(opts.linecolor, ["Coral", "OrangeRed",...
        "Gold", "Yellow", "Khaki", "DarkKhaki", "IndianRed", "Salmon", ...
        "Red", "FireBrick", "DarkRed", "Pink", "HotPink", "DeepPink", ...
        "Orange", "Violet", "Magenta", "BlueViolet", "DarkOrchid", ...
        "Purple", "Indigo", "SlateBlue", "MediumSlateBlue", "GreenYellow",...
        "LawnGreen", "Lime", "LimeGreen", "SpringGreen", "MediumSeaGreen", ...
        "SeaGreen", "Green", "YellowGreen", "Olive", "DarkCyan", ...
        "Teal", "Cyan", "SteelBlue", "SkyBlue", "DeepSkyBlue", ...
        "DodgerBlue", "RoyalBlue", "Blue", "MediumBlue", "Navy", ...
        "Bisque", "Wheat", "BurlyWood", "Tan", "Goldenrod", "Peru", ...
        "Chocolate", "Brown", "Maroon", "White", "HoneyDew", "Azure", ...
        "WhiteSmoke", "Beige", "Ivory", "Linen", "LightGray", "Silver", ...
        "DarkGray", "Gray", "SlateGray", "DarkSlateGray", "Black", ...
        "BlueLight", "BrutalBlue", "CrownYellow", "YassQueen", ...
        "SisterSister"])} = "DarkSlateGray"
    opts.lineWidth (1,1) double = 1.1

    % options for adjusting p-values for multiple testing
    opts.padjustScheme {mustBeTextScalar, mustBeMember(opts.padjustScheme, ...
        ["overall", "pertrait"])} = "pertrait" % if "pertrait", applies 'padjust' to each trait separately
    opts.padjust {mustBeTextScalar, mustBeMember(opts.padjust, ["holm", ...
        "hochberg", "hommel", "bonferroni", "BH", "BY", "none" ])} = "none"
    opts.padjustAlpha = 0.05 % adjusted P-value cut-off

    % options for effect sizes
    opts.colormap {mustBeMember(opts.colormap, ["jet", "turbo", "abyss", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", "summer", ...
        "colorcube", "lines", "magma", "plasma", "inferno", "cividis", ...
        "viridis"])} = "plasma" % colormap for effect size
    opts.binEffect (1,1) logical = true % effects will be divided in two groups of risk inscreasing/decreasing.

    % aesthetic
    opts.xOrder {mustBeText, mustBeVector} % string vector for the order of xCol values to be shown
    opts.markerAlpah (1,1) double = 0.9
    opts.markerEdgeAlpah (1,1) double = 0.6
    opts.fontsize (1,1) double = 11
    opts.Xfontsize (1,1) double % if not used 'fontsize' option is used
    opts.Yfontsize (1,1) double % if not used 'fontsize' option is used
    opts.fontname {mustBeTextScalar, mustBeMember(opts.fontname, ...
        ["Bauhaus", "Futura", "Bodoni", "Garamond", "Corbel", ...
        "Corbel Light", "Arial", "Tahoma"])} = "Arial"
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "png" % format of output plots
    opts.resolution (1, 1) double = 300 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    opts.xlabel {mustBeTextScalar} = ""
    opts.ylabel {mustBeTextScalar} = ""
    opts.xlabelAngel (1,1) double = 45
    opts.output {mustBeTextScalar} = "effDotPlot"

    opts.excludeNonSig (1,1) logical = false % to exclude non-significant associations (treated like missing: shown as X).
    opts.title {mustBeTextScalar} = "" % plot title 
    opts.rmmissing (1,1) logical = true % to remove all missngness by default
end

% create a new table 
df = struct;
df.Trait = tab.(opts.xCol);
df.P = tab.(opts.pCol);
df.Beta = tab.(opts.bCol);
df.Locus = tab.(opts.yCol);
df = struct2table(df);
if opts.rmmissing, df = rmmissing(df); end

if isfield(opts, "xOrder") % remove traits not requested
    df(~ismember(df.Trait, opts.xOrder), :) = [];
end

% binarize effect
if opts.binEffect
    idx_pos = df.Beta > 0;
    idx_neg = df.Beta < 0;
    df.Beta(idx_pos) = 1; df.Beta(idx_neg) = -1;

    try
        cmap = get_colormap(opts.colormap, 256);
    catch
        cmap = feval(opts.colormap);
    end
    idx = round([0.1, 0.9].*size(cmap, 1)); % keep 2 colors only
    opts.cmap = cmap(idx, :);
end

% adjust p-values for multiple testing
df.marker(:) = "O"; 
if opts.padjust ~= "none"
    
    if opts.padjustScheme == "pertrait" % adjust p in each trait category
        ut = unique(df.Trait);
        for k = 1:numel(ut)
            idx = df.Trait == ut(k);
            padj = padjust(df.P(idx), method=opts.padjust);
            df.P(idx) = padj;
            clear padj
        end

    else % overll adjustment
        df.P = padjust(df.P, method=opts.padjust);
    end
   
end

% non-significant p-adj will be shown with a different marker
dash_idx = df.P > opts.padjustAlpha;

if opts.excludeNonSig
    df(dash_idx, :) = [];
else
    df.marker(dash_idx) = "square";
end

% add markerAlpha to df table
df.markerAlpah(:) = opts.markerAlpah;
df.markerAlpah(df.marker == "square") = 0.2;
df.Beta(df.marker == "square") = 0;

% Since some p-value cannot go under realmin, modify 0 p-values 
f_zero = df.P <= realmin;
df.P(f_zero) = realmin*eps;
df.log10p = -log10(df.P); % for legend

% set size factors for P-values
df.P = rescale(df.log10p, opts.szRange(1), opts.szRange(2)); % for dots

% use percentiles
% df.P = stratify(df.log10p, "q", (1:100)./100);
% df.P = rescale(df.P, opts.szRange(1), opts.szRange(2)); % for dots


% size of plot
yn = numel(unique(df.Locus));
xn = numel(unique(df.Trait));

% draw now
close all force
if all(~isempty(opts.WindowPosition) & ~ismissing(opts.WindowPosition))
    fig = figure(Position=opts.WindowPosition, Units="pixels", Visible=opts.visible);
elseif isfield(opts, "WindowState")
    fig = figure("WindowState", opts.WindowState, Visible=opts.visible);
else
    fig = figure(Visible=opts.visible);
end

ax = gca(fig);
% axis(ax, "equal");
ax.XLim = [0.5, xn+0.5];
ax.YLim = [0.5, yn+0.5];
ax.XTick = 1:xn; ax.YTick = 1:yn;

% make y-axis long enough
drawnow;
pos = ax.Position;
ax.Position(4) = 0.98 - pos(2);

pbaspect(ax, opts.pbaspect); % aspect ratio of plot

% draw constant horizontal lines connecting dots --------------------------
plt = palette;
opts.linecolor = hex2rgb(plt{string(opts.linecolor), "code"});
xx = repmat(1:xn, yn, 1);
yy = repmat(1:yn, xn, 1);
line(ax, xx', yy, Color=opts.linecolor, LineWidth=opts.lineWidth);


% scatter plot of dots with size proportion of -log10 P (rescaled)
hold(ax, "on")
ut = natsort(unique(df.Trait, "Stable"));
ul = natsort(unique(df.Locus, "Stable"));

if isfield(opts, "xOrder")
    opts.xOrder(~ismember(opts.xOrder, ut)) = [];
    ut = opts.xOrder;
end

for k1 = 1:yn
    
    for k2 = 1:xn
    
        idx = df.Trait == ut(k2) & df.Locus == ul(k1);
        viz = df(idx, :);
        if isempty(viz) % no summary stat for this loci/trait pair
            scatter(ax, k2, k1, mean(opts.szRange), ...
                [0, 0, 0], "filled", ...
                "MarkerEdgeColor", "black", "Marker", "x", ...
                "LineWidth", 2);
        else

            if opts.binEffect
                if viz.marker == "square" || ismissing(viz.Beta) % no color (non-significant)
                    tmp_color = [0, 0, 0];
                elseif viz.Beta < 0
                    tmp_color = opts.cmap(1, :);
                else
                    tmp_color = opts.cmap(2, :);
                end
                
                scatter(ax, k2, k1, viz.P, "filled", ...
                    "MarkerEdgeColor", "black", "Marker", viz.marker, ...
                    "MarkerEdgeAlpha", opts.markerEdgeAlpah, ...
                    "MarkerFaceColor", tmp_color, ...
                    "MarkerFaceAlpha", viz.markerAlpah);
            else % continuous effect size
                scatter(ax, k2, k1, viz.P, viz.Beta, "filled", ...
                    "MarkerEdgeColor", "black", "Marker", viz.marker, ...
                    "MarkerEdgeAlpha", opts.markerEdgeAlpah, ...
                    "MarkerFaceAlpha", viz.markerAlpah);
            end

        end

    end
end


% apply colormap
if opts.binEffect
    colormap(ax, opts.cmap);
else
    try
        opts.color = feval(opts.colormap);
    catch
        opts.color = get_colormap(opts.colormap, 256);
    end
    colormap(ax, opts.color);
end

% add ytick and xtick lables
ax.XTickLabel = ut; ax.YTickLabel = ul;
ax.XTickLabelRotation = opts.xlabelAngel;
ax.XLabel.String = opts.xlabel;
ax.YLabel.String = opts.ylabel;
ax.FontName = opts.fontname;

if isfield(opts, "Xfontsize")
    ax.XAxis.FontSize = opts.Xfontsize;
else
    ax.XAxis.FontSize = opts.fontsize;
end

if isfield(opts, "Yfontsize")
    ax.YAxis.FontSize = opts.Yfontsize;
else
    ax.YAxis.FontSize = opts.fontsize;
end

ax.TickDir = "out";
ax.TickLength = [0, 0];
ax.Box = "on";

ax.Title.String = opts.title;
ax.Title.FontName = opts.fontname;
ax.Title.FontSize = opts.fontsize + 2;

% add/adjust colorbar
c = colorbar(ax, 'Location', 'eastoutside');
% drawnow;
c.Position(4) = c.Position(4)/4;
c.Position(3) = c.Position(3)/2;
tpos = ax.tightPosition;
c.Position(1) = tpos(1) + tpos(3) + 0.001;
c.Position(2) = tpos(2) + tpos(4) - c.Position(4);
c.AxisLocation  = 'out';
c.Label.String = opts.bLabel;
% c.Label.Units = 'normalized';
% c.Label.HorizontalAlignment = 'left';
% c.Label.VerticalAlignment = 'middle';
c.FontSize = opts.fontsize - 2;
c.FontName = opts.fontname;
c.Label.FontSize = opts.fontsize - 2;
c.Label.FontName = opts.fontname;
% c.Label.Position = [1 0.5 0];

if opts.binEffect
    c.Ticks = [0.25, 0.75];
    c.TickLabels = ["-", "+"];
end

% add legend for p-values
% pvec_str = unique([min(df.log10p), mean(df.log10p), max(df.log10p)]); % for legend: actual p-values
% pvec_line = unique([min(df.P), mean(df.P), max(df.P)]);
% pvec_line = 2 * sqrt(pvec_line / pi);
% % pvec_line = sqrt(pvec_line);  % scatter to line sizes
% lg = (0);
% for k = 1:numel(pvec_line)
%     lg(k) = plot(ax, nan, nan, 'MarkerFaceColor', [0.3, 0.3, 0.3],...
%         'MarkerEdgeColor', 'k', 'LineStyle', 'none',...
%         'Marker', "O", "MarkerSize", pvec_line(k));
% end
% 
% leg = legend(lg, 'String', compose("%.1E", string(10.^-pvec_str)), ...
%     'Location', 'northeast',...
%     'Color', '#BAD8E0', "AutoUpdate", false);
% leg.Title.String = 'P-value';
% leg.TextColor = '#3F3F3F';
% leg.Position(1) = c.Position(1);
% leg.Position(2) = c.Position(2) - leg.Position(4) - 0.01;

pvec_scatter = unique([min(df.P), mean(df.P), max(df.P)]);
% Generate the text labels based on the original -log10(p) values
pvec_labels = compose("%.1E", string(10.^-unique([min(df.log10p), mean(df.log10p), max(df.log10p)]) ));

lg_width  = 0.06;
lg_height = 0.1;
lg_pos_x  = c.Position(1);
lg_pos_y  = c.Position(2) - lg_height - 0.01;

axLeg = axes('Position', [lg_pos_x, lg_pos_y, lg_width, lg_height], ...
             'Color', '#BAD8E0', 'Box', 'on'); % box creates the surrounding border
hold(axLeg, 'on');
axis(axLeg, [0 1 0 1]); 
axLeg.TickLength = [0 0];
axLeg.XTick = [];
axLeg.YTick = [];

% Define positions for the dummy markers and the associated text
x_marker = 0.2;  % x position of the markers in the custom legend
x_text   = 0.4;  % x position for the labels
y_coords = linspace(0.8, 0.2, numel(pvec_scatter)); % evenly spaced vertical positions

% Plot the dummy markers and labels
for k = 1:numel(pvec_scatter)
    % Draw a marker (scatter uses area values directly)
    scatter(axLeg, x_marker, y_coords(k), pvec_scatter(k), ...
        [0.3, 0.3, 0.3], 'filled', ...
        'MarkerEdgeColor', 'k');

    % Add a text label next to the marker
    text(axLeg, x_text, y_coords(k), pvec_labels(k), ...
        VerticalAlignment="middle", FontSize=opts.fontsize-2, ...
        Color="#3F3F3F", FontName=opts.fontname);
end

% Add a title for the custom legend
axLeg.Title.String = "P-value";
axLeg.Title.FontName = opts.fontname;
axLeg.Title.FontSize = opts.fontsize - 2;

% Move downwards to avoid the collision with colorbar
tipos = axLeg.TightInset + axLeg.Position;
axLeg.Position(2) = c.Position(2) - tipos(4) - 0.01;

% Maker sure that the legend does cover all text labels
teObj = findobj(axLeg, "type", "Text");
teExt = reshape([teObj.Extent], 4, []);
while max(sum(teExt([1, 3], :), 1)) > 0.95
    
    % extend axLeg X-position
    is_within = axLeg.Position;
    is_within = (is_within(1) + is_within(3)) < 0.99; % only if so, otherwise reduce text fontsize

    if is_within
        axLeg.Position(3) = axLeg.Position(3) + 5e-3;
    else
        for k = 1:numel(teObj)
            teObj(k).FontSize = teObj(k).FontSize - 1;
        end
    end

    drawnow;
    teExt = reshape([teObj.Extent], 4, []);

end

% % add legend
% if any(df.fade) % any square?
%     lg = (0);
%     lg(1) = plot(NaN, NaN, 'MarkerFaceColor', [0.7, 0.7, 0.7],...
%                 'LineStyle', 'none', 'Marker', "O",...
%                 'MarkerEdgeColor',  'k',...
%                 'MarkerSize', mean(opts.szRange));
%     lg(2) = plot(NaN, NaN, 'MarkerFaceColor', [0.7, 0.7, 0.7],...
%                 'LineStyle', 'none', 'Marker', "square",...
%                 'MarkerEdgeColor',  'k',...
%                 'MarkerSize', mean(opts.szRange));
% 
%     leg = legend(lg, ["Sig", "Non-sig"], 'Color', [0.93 0.93 0.93], ...
%         'AutoUpdate', 'off',...
%         'FontSize', 11, 'Location', 'northeast', ...
%         'Orientation', 'Vertical');
% 
% end

hold(ax, "off")
shade(ax, FaceAlpha=0.1, dir="xy");

if opts.save
    if any(ismember(opts.format.lower, ["pdf", "eps"]))
        exportgraphics(fig, opts.output + "." + opts.format, ...
            "ContentType", "vector")
    else
        exportgraphics(fig, opts.output + "." + opts.format, ...
            'Resolution', opts.resolution)
    end
end

if opts.savefig 
    savefig(fig, opts.output + ".fig")
end

end % END

%% custom legend TODO
% ----- Custom Legend Setup (all in points) -----
% 
% % Get the scatter marker sizes (areas in points^2) as used in your main plot
% pvec_scatter = unique([min(df.P), mean(df.P), max(df.P)]);
% % Generate the text labels based on the original -log10(p) values
% pvec_labels = compose("%.1E", string(10.^-unique([min(df.log10p), mean(df.log10p), max(df.log10p)]) ));
% 
% % Compute the marker diameters (in points) for each marker area.
% % (Remember: scatter marker size is an area, so diameter = 2*sqrt(area/π).)
% diameters = 2 * sqrt(pvec_scatter / pi);
% 
% % Define a vertical gap between markers as a fraction of the smallest diameter.
% % (Using 10% of the smallest diameter as the gap keeps it relative to marker size.)
% gap_vert = 0.1 * min(diameters);
% 
% % Define padding for the legend axes (both horizontally and vertically).
% % This will be used to ensure that both markers and texts are inside the box.
% pad_x = gap_vert;    % use the same fraction (10% of smallest diameter)
% pad_y = gap_vert;
% 
% % Compute the vertical stacking (in points) of the markers.
% nMarkers = numel(diameters);
% % Total height needed for markers (without the top and bottom padding)
% totalMarkersHeight = sum(diameters) + gap_vert * (nMarkers - 1);
% % Now, the full legend height includes the vertical padding at top and bottom:
% totalLegendHeight = totalMarkersHeight + 2 * pad_y;
% 
% % Compute the y-center positions for each marker so that markers are stacked
% % starting from pad_y at the bottom. (They will be placed at the center of their own space.)
% yCenters = zeros(nMarkers, 1);
% yCenters(1) = pad_y + diameters(1)/2;
% for i = 2:nMarkers
%     yCenters(i) = yCenters(i-1) + (diameters(i-1)/2 + gap_vert + diameters(i)/2);
% end
% % (yCenters(end) + diameters(end)/2 equals totalMarkersHeight + pad_y)
% 
% % For horizontal placement, we want the marker to be fully inside the box.
% % Let the marker be centered at x_marker.
% % A natural choice is to set x_marker so that the left edge is at pad_x:
% %   x_marker - diameters(k)/2 >= pad_x  -> simplest: x_marker = pad_x + (max(diameters)/2)
% x_marker = pad_x + max(diameters)/2;
% 
% % Now, for the text label. We want the text to appear to the right of the largest marker.
% % We choose x_text to be just after the right edge of the largest marker plus a small gap.
% % In other words, x_text = (x_marker + max(diameters)/2) + gap_hor.
% % Here we set the horizontal gap (gap_hor) relative to the marker size.
% gap_hor = 0.1 * max(diameters);  
% x_text = x_marker + max(diameters)/2 + gap_hor;
% 
% % Next, we need to determine the width required for the text.
% % We create a temporary text object (in the new units "points") with the longest label,
% % then get its extent, and then delete it.
% [~, idxLongest] = max(cellfun(@length, cellstr(pvec_labels)));
% tempText = text('Units','points','Position',[0,0],'String', pvec_labels(idxLongest), ...
%     'FontSize', opts.fontsize-2, 'FontName', opts.fontname);
% tempExtent = get(tempText, 'Extent');  % format: [x, y, width, height]
% textWidth = tempExtent(3);
% delete(tempText);
% 
% % Now compute the total legend width.
% % It must include left padding, the marker area (which, horizontally, spans max(diameters)),
% % the horizontal gap, the width of the text, and right padding.
% totalLegendWidth = pad_x + max(diameters) + gap_hor + textWidth + pad_x;
% 
% % ----- Create Custom Legend Axes -----
% 
% % Create an axes that will be our custom legend.
% % We start with a temporary normalized position (we will override it in points).
% lg_pos_norm = [c.Position(1), c.Position(2) - 0.05, 0.06, 0.12];  % temporary values
% axLeg = axes('Position', lg_pos_norm, 'Color', '#BAD8E0', 'Box', 'on');
% 
% % Switch the axes units to points so that we can set size exactly.
% axLeg.Units = 'points';
% % Set the axes limits to exactly accommodate the computed width and height.
% axLeg.XLim = [0, totalLegendWidth];
% axLeg.YLim = [0, totalLegendHeight];
% 
% % ----- Position the Legend Axes Relative to the Colorbar -----
% 
% % Get the colorbar position in points.
% oldUnits = c.Units;
% c.Units = 'points';
% cPos = c.Position;  % [x, y, width, height] in points
% c.Units = oldUnits;
% 
% % In our previous approach the legend was positioned flush with c's top,
% % but here we want the legend to appear below the colorbar.
% % So we want the legend's top (i.e. newLegend_y + totalLegendHeight) to be equal to cPos(2) minus a small gap.
% % We can use gap_hor (or pad_y) as the gap between the colorbar and legend.
% newLegend_y = cPos(2) - pad_y - totalLegendHeight;  
% % For horizontal position, align the left edge with the colorbar's x-position.
% newLegend_x = cPos(1);
% % Set the new position (in points) for axLeg.
% axLeg.Position = [newLegend_x, newLegend_y, totalLegendWidth, totalLegendHeight];
% 
% % ----- Plot the Dummy Markers and Labels -----
% 
% hold(axLeg, 'on');
% for k = 1:nMarkers
%     % Plot the marker with its area using scatter.
%     scatter(axLeg, x_marker, yCenters(k), pvec_scatter(k), [0.3, 0.3, 0.3], ...
%         'filled', 'MarkerEdgeColor', 'k');
%     % Place the text to the right of the marker.
%     text(axLeg, x_text, yCenters(k), pvec_labels(k), ...
%          'VerticalAlignment', 'middle', 'FontSize', opts.fontsize-2, ...
%          'Color', '#3F3F3F', 'FontName', opts.fontname);
% end
% 
% % Add a title to the legend at the top center.
% % We position the title within the axes so it sits inside the box.
% title(axLeg, 'P-value', 'FontName', opts.fontname, 'FontSize', opts.fontsize-2);
% 
% % Optionally, remove tick marks from the custom legend axes.
% axLeg.XTick = [];
% axLeg.YTick = [];

