function gax = effNetPlotter(tab, opts)
% similar to 'effPlotter' function, but visualizes the effects in a graph.
% Inspired by https://www.nature.com/articles/s41588-018-0047-6. 
% 
% Oveis Jamialahmadi, University of Gothenburg, 15 MAY 2023.
% 
% Tips: - for large networks with "force" layout, use WeightEffect="direct"
%         and UseGravity=true.
%       - 'nominal' option was designed only for combination of adjusted
%          and nominally significant p-values. In case no adjustment is
%          desired, set 'padjust' to "none" and 'padjustAlpha' to the
%          desired cutoff (e.g. 5e-8).
% 
% @21JAN2024: colormaps were updated.

arguments
    tab {mustBeA(tab, "table")} % results of cell-type specific analysis from LDSC. Example table can be found within the same directory: ldsc_cell_type_specific.mat
    opts.pCol {mustBeTextScalar} % = "P-adj" % p-values column
    opts.snpCol {mustBeTextScalar} % = "loci" % variants column
    opts.locusCol {mustBeTextScalar} % = "loci" % locus column (optional): to group 'snpCol' into different categories (e.g. gene/locus names)
    opts.traitCol {mustBeTextScalar} % = "Pheno" % traits column
    opts.traitCatCol {mustBeTextScalar} % = "Pheno" % trait categories column
    opts.effectCol {mustBeTextScalar} % = "β|OR"
    opts.OReffs logical {mustBeVector} % vector of odds ratio effects (same size as input 'tab'). If left empty (default), all effects are treated as regression coefficients.
    opts.effectBinarize (1,1) logical = true % binarize effect sizes (> 0 and < 0)
    
    % P adjustment
    opts.padjust {mustBeTextScalar, mustBeMember(opts.padjust, ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none" ])} = "BH"
    opts.padjustAlpha = 0.05 % adjusted P-value cut-off
    opts.nominal (1,1) double = 1 % nominal p-value cutoff to include also nominal p-values (shown as dashed edges). Default is 1, i.e. ignores nominally significant associations.
    opts.pVis {mustBeMember(opts.pVis, ["adjusted", "nominal"])} = "adjusted" % to visualize either adjusted or nomianl p-values for edges' weights (if padjust is "none", will be "nominal")

    opts.edgeFactor (1, 2) double = [0.3 3] % min, max width of edges
    opts.nodeFactor = 9 % adds this value to node size
    opts.maxNodeSize = 17 % max node size: if any node size (upon the addition of 'nodeFactor') > this value, node sizes will be capped by this value.
    opts.mapSnp (1,1) {mustBeMember(opts.mapSnp, ["jet", "turbo", "abyss", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", "summer", ...
        "colorcube", "lines", "magma", "plasma", "inferno", "cividis", ...
        "viridis"])} = "magma"
    opts.mapTrait(1,1) {mustBeMember(opts.mapTrait, ["jet", "turbo", "abyss", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", "summer", ...
        "colorcube", "lines", "magma", "plasma", "inferno", "cividis", ...
        "viridis"])} = "parula"
    opts.mapEffect (1,1) {mustBeMember(opts.mapEffect, ["jet", "turbo", "abyss", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", "summer", ...
        "colorcube", "lines", "magma", "plasma", "inferno", "cividis", ...
        "viridis"])} = "jet"
    opts.fontname {mustBeTextScalar} = "Garamond" % node font name
    opts.fontsize (1,1) double = 11 % node font size
    opts.fontsizeLegend (1,1) double = 9 % legend's font size
    opts.layout {mustBeMember(opts.layout, ["force", "circle", "layered"])} = "force"
    opts.UseGravity (1,1) logical = false
    opts.WeightEffect {mustBeMember(opts.WeightEffect, ["none", "inverse", "direct"])} = "none" % for layoyt == "force"
    opts.Iterations (1,1) double = 1e4
    opts.Direction {mustBeMember(opts.Direction, ["down", "up", "left", "right"])} = "right"
    opts.AssignLayers {mustBeMember(opts.AssignLayers, ["auto", "asap", "alap"])} = "auto"
    opts.BezierPoints (1,1) double = 60 % number of points used for Bezier curve
    opts.d (1,1) double = 0.2 % Bezier control point distance factor
    opts.ArrowHeadLength (1,1) double = 5
    opts.ArrowHeadWidth (1,1) double = 5

    opts.output {mustBeTextScalar} = "effNetPlot"
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "jpg" % format of output plots
    opts.resolution (1, 1) double = 300 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well

end

% convert ORs to log ORs
if ~isfield(opts, "OReffs")
    opts.OReffs = false(height(tab), 1);
end
tab.(opts.effectCol)(opts.OReffs) = log(tab.(opts.effectCol)(opts.OReffs));

% filter based on 'fdr' cutoff
tab.p = tab.(opts.pCol);
tab.p(tab.p < realmin*eps) = realmin*eps;
tab.q = padjust(tab.p, method=opts.padjust);
if isfield(opts, 'nominal') && opts.nominal < 1
    tab(tab.p > opts.nominal, :) = [];
    tab.nominal = tab.q < opts.padjustAlpha;
    opts.nominal = true;
else
    opts.nominal = false;
    tab(tab.q > opts.padjustAlpha, :) = [];
end

% weights for the graph's edges
if opts.pVis == "nominal"
    opts.legText = "P";
    tab.q = tab.p;
else
    opts.legText = "P^{adj}";
end

if isempty(tab)
    disp('no significant results remained after adjusting for FDR!')
    gax = [];
    return
end

% trait and variant counts 
cnt.trait = groupsummary(tab, opts.traitCol);
cnt.ct = groupsummary(tab, opts.snpCol);

if ~isfield(opts, "traitCatCol"), opts.traitCatCol = opts.traitCol; end
if ~isfield(opts, "locusCol"), opts.locusCol = opts.snpCol; end
tab.(opts.locusCol) = string(tab.(opts.locusCol));
tab.(opts.traitCatCol) = string(tab.(opts.traitCatCol));

% add colors
uloci = unique(tab.(opts.locusCol), "stable");
try
    col_loci = get_colormap(opts.mapSnp, numel(uloci));
catch
    col_loci = feval(opts.mapSnp, numel(uloci));
end
% col_loci = feval(opts.mapSnp, numel(uloci));
[~, idx] = ismember(tab.(opts.locusCol), uloci);
tab.col_loci = col_loci(idx, :);

ucats = unique(tab.(opts.traitCatCol), "stable");
try
    col_traits = get_colormap(opts.mapTrait, numel(ucats));
catch
    col_traits = feval(opts.mapTrait, numel(ucats));
end
% col_traits = feval(opts.mapTrait, numel(ucats));
[~, idx] = ismember(tab.(opts.traitCatCol), ucats);
tab.col_trait = col_traits(idx, :);

% create a digraph
G = digraph(tab.(opts.snpCol), tab.(opts.traitCol), -log10(tab.q));

% for each trait-cell pair, keep only the strongest association per each marker 
G = simplify(G, 'max');

% edge withs correpond to -log10 FDR
% LWidths = opts.edgeFactor*G.Edges.Weight/max(G.Edges.Weight);
LWidths = rescale(G.Edges.Weight, opts.edgeFactor(1), opts.edgeFactor(2));

% plot the graph
close all force
switch opts.layout
    case "force"
        ax = plot(G, LineWidth=LWidths, Layout=opts.layout, ...
            UseGravity=opts.UseGravity, Iterations=opts.Iterations, ...
            WeightEffect=opts.WeightEffect);
    case "layered"
        ax = plot(G, LineWidth=LWidths, Layout=opts.layout, ...
            Direction=opts.Direction, AssignLayers=opts.AssignLayers);
    otherwise
        ax = plot(G, LineWidth=LWidths, Layout=opts.layout);
end

% arrow size
ax.ArrowSize = LWidths.*3;
ax.ArrowPosition = 0.7;
ax.EdgeAlpha = 0.95;

% square for cell-type and circle for traits
idx = ismember(ax.NodeLabel, tab.(opts.traitCol));
ax.Marker = repmat("O", numel(idx), 1);
ax.Marker(~idx) = {'square'};

% node size based on counts
markersize = ax.MarkerSize;

[f1, f2] = ismember(ax.NodeLabel, cnt.trait.(opts.traitCol)); 
ax.MarkerSize = ones(numel(f1), 1).*markersize;
ax.MarkerSize(f1) = cnt.trait.GroupCount(f2(f1)); 

[f1, f2] = ismember(ax.NodeLabel, cnt.ct.(opts.snpCol)); 
ax.MarkerSize(f1) = cnt.ct.GroupCount(f2(f1));

ax.MarkerSize = ax.MarkerSize + opts.nodeFactor;

% cap node sizes
if any(ax.MarkerSize > opts.maxNodeSize)
    ax.MarkerSize = rescale(ax.MarkerSize, min(ax.MarkerSize), opts.maxNodeSize);
end

% node colors
[f1, f2] = ismember(ax.NodeLabel, tab.(opts.traitCol));
ax.NodeColor = repmat(tab.col_loci(1, :), numel(f1), 1);
ax.NodeColor(f1, :) = tab.col_trait(f2(f1), :);

[f1, f2] = ismember(ax.NodeLabel, tab.(opts.snpCol));
ax.NodeColor(f1, :) = tab.col_loci(f2(f1), :);

%% generate a new graph plot for better control over the figure behaviors--
fig = figure(WindowState="maximized");
h = axes(fig);
hold(h, "on")
X = ax.XData;
Y = ax.YData;

h = gscatter(h, X, Y, 1:numel(X));
ax2 = h(1).Parent;
ax2.Position(1) = 0.05; % move the axis to the left, to make extra space for legends to the right
ax2.Legend.Visible = "off";

for k = 1:numel(h)
    h(k).Marker = ax.Marker(k);
    h(k).MarkerSize = ax.MarkerSize(k);
    h(k).MarkerFaceColor = ax.NodeColor(k, :);
    h(k).MarkerEdgeColor = 'k';
    h(k).Tag = "Node";
end
hold(ax2, "off")

% add legends for nodes: traits and variant groups ------------------------
lg = (0);
tmp = unique(tab(:, [opts.traitCatCol, "col_trait"]), "rows");
hold(ax2, "on")
for k = 1:height(tmp)
    lg(k) = plot(NaN, NaN, MarkerFaceColor=tmp.col_trait(k,:),...
            LineStyle='none', Marker='O', MarkerEdgeColor='k',...
            Color=tmp.col_trait(k, :), MarkerSize=mean(ax.MarkerSize));
end
leg = legend(lg, tmp.(1), Color=[0.93 0.93 0.93], AutoUpdate='off',...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    Location='northeastoutside', Orientation='Vertical');
leg.Title.String = "Trait";
hold(ax2, "off")

tmp = unique(tab(:, [opts.locusCol, "col_loci"]), "rows");
ax3 = copyobj(ax2, gcf);
ax3.Visible = 'off';
drawnow;
for i = 1:numel(ax3.Children)
    ax3.Children(i).Visible = 'off';
end
drawnow;

hold(ax3, "on")
lg2 = (0);
for k = 1:height(tmp)
    lg2(k) = plot(ax3, NaN, NaN, MarkerFaceColor=tmp.col_loci(k,:),...
            LineStyle='none', Marker='s', MarkerEdgeColor='k',...
            Color=tmp.col_loci(k, :), MarkerSize=mean(ax.MarkerSize));
end

leg2 = legend(lg2, tmp.(1), Color=[0.93 0.93 0.93], AutoUpdate='off',...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    Location='southeastoutside', Orientation='Vertical');
leg2.Title.String = "Variant";
hold(ax3, "off")

drawnow
ax3.Position = ax2.Position;

% shift the second legend upwards
drawnow;
leg2.Position(2) = leg.Position(2) - leg2.Position(4) - 0.01;

ncol = 2;
while leg2.Position(2) < 0.45 % too low so other legens cannot fit
    drawnow;
    leg2.NumColumns = ncol;
    leg2.Position(1) = leg.Position(1);
    leg2.Position(2) = leg.Position(2) - leg2.Position(4) - 0.01;
    drawnow;
    ncol = ncol + 1;
end

% legend for marker size
nodeinfo.marker = string(ax.Marker);
nodeinfo.size = ax.MarkerSize;

% create a new axis for trait/variant marker sizes ------------------------
ax4 = copyobj(ax3, gcf);
ax4.Visible = "on";
drawnow;
for i = 1:numel(ax4.Children)
    ax4.Children(i).Visible = 'off';
end

drawnow;
ax4.Position(1) = leg2.Position(1);
ax4.Position(2) = leg2.Position(2) - ax4.Position(4) - 0.01;
ax4.Position(4) = ax4.Position(4) + ax4.Position(2) - ax3.Position(2);
ax4.Position(2) = ax3.Position(2);
ax4.Position(3) = max(leg.Position(3), leg2.Position(3)); 

% equal axis: commented, better to use unequal axes
drawnow;
ax4.Position(4) = ax4.Position(3);
afactor = ax4.PlotBoxAspectRatio(1)/ax4.PlotBoxAspectRatio(2);
ax4.Position(4) = ax4.Position(4)*afactor;
ax4.Position(2) = leg2.Position(2) - ax4.Position(4) - 0.01;

% add two annotation text boxes for traits and variants
% nodew = getMarkerBounds(h); % non-overlapping (horizontally) markers for legend
drawnow;
axtest = plot(ax4, 1, 0.5, Marker="O", MarkerSize=max(nodeinfo.size));
te =  text(ax4, 1, 0.5, "Variants per trait", FontSize=opts.fontsizeLegend, FontName=opts.fontname); % test text title
teExt = te.Extent; % how much space do we need for text title?
teExt = 1.05*teExt(4); % the space we need for the title string
nodew = getMarkerBounds(ax4.Children);
delete(axtest)
delete(te)
% nodew.w = nodew.w.*.5;
% nodew.h = nodew.h.*.5;

maxMarginX = max(nodew.w); % max width of markers aligned horizontally
ax4.XLim = [1-maxMarginX/2, 1+2.2*maxMarginX];
midX = ax4.XLim(1) + diff(ax4.XLim)/2;

ax4.YLim = [0, 1]; % initial YLim height
maxMarginY = max(nodew.h); % max height of annotation boxes should be decided based on this value

ax4.YLim = [0, 2*teExt + 2.2*maxMarginY];

pos1 = [ax4.XLim(1), 1.2*maxMarginY+teExt, diff(ax4.XLim), maxMarginY+teExt]; % upper rectangle
pos2 = [ax4.XLim(1), 0, diff(ax4.XLim), maxMarginY+teExt]; % lower rectangle
rectangle(ax4, Position=pos1);
yline(ax4, 2.2*maxMarginY+teExt);
text(ax4, midX, 2.2*maxMarginY+1.5*teExt, "Variants per trait", ...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    VerticalAlignment="middle", HorizontalAlignment="center", FontWeight="bold");

rectangle(ax4, Position=pos2);
yline(ax4, maxMarginY);
text(ax4, midX, maxMarginY+0.5*teExt, "Traits per variant", ...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    VerticalAlignment="middle", HorizontalAlignment="center", FontWeight="bold");

idx = nodeinfo.marker == "o";
msize1 = [max(nodeinfo.size(idx)), median(nodeinfo.size(idx)), min(nodeinfo.size(idx))];
msize2 = [max(nodeinfo.size(~idx)), median(nodeinfo.size(~idx)), min(nodeinfo.size(~idx))];
mlabel1 = round([max(cnt.trait.GroupCount), median(cnt.trait.GroupCount), min(cnt.trait.GroupCount)]);
mlabel2 = round([max(cnt.ct.GroupCount), median(cnt.ct.GroupCount), min(cnt.ct.GroupCount)]);

[mlabel1, idx] = unique(mlabel1, "stable");
msize1 = msize1(idx);
[mlabel2, idx] = unique(mlabel2, "stable");
msize2 = msize2(idx);

hold(ax4, "on")
for k = 1:numel(msize1)
    plot(ax4, 1 + maxMarginX*(k-1), teExt+1.7*maxMarginY, Marker="O", MarkerSize=msize1(k), Color="k");
    text(ax4, 1 + maxMarginX*(k-1), teExt+1.7*maxMarginY, string(mlabel1(k)), ...
        VerticalAlignment="middle", HorizontalAlignment="center", ...
        FontName=opts.fontname, FontSize=opts.fontsizeLegend-2);
end

for k = 1:numel(msize2)
    plot(ax4, 1 + maxMarginX*(k-1), maxMarginY/2, Marker="S", MarkerSize=msize2(k), Color="k");
    text(ax4, 1 + maxMarginX*(k-1), maxMarginY/2, string(mlabel2(k)), ...
        VerticalAlignment="middle", HorizontalAlignment="center", ...
        FontName=opts.fontname, FontSize=opts.fontsizeLegend-2);
end
hold(ax4, "off")
ax4.Visible = "off";

% another axis for -log10(FDR)---------------------------------------------
ax5 = copyobj(ax4, gcf);
ax5.Visible = "on";
drawnow;
for i = 1:numel(ax5.Children)
    ax5.Children(i).Visible = 'off';
end
ax5.Position(1) = ax4.Position(1);
ax5.Position(2) = ax3.Position(2);
ax5.Position(4) = max(1e-6, ax4.Position(2) - ax5.Position(2) - 0.015);
df = ax5.Position(4) - ax5.Position(3);
if df > 0
    ax5.Position(4) = ax5.Position(3);
    ax5.Position(2) = ax5.Position(2) + df;
end
log10fdr = ceil(-log10(tab.q));
log10fdr = quantile(unique(log10fdr), [0, 0.4, 0.8, 1]); %ceil([min(log10fdr), median(log10fdr), max(log10fdr)]);
[log10fdr, uidx] = unique(round(log10fdr), "stable");
lw = unique(round(ax.LineWidth, 3, "significant"));
log10fdrw = quantile(lw, [0, 0.4, 0.8, 1]);
log10fdrw = log10fdrw(uidx);

for k = 1:numel(log10fdr)
    yline(ax5, log10fdr(k), LineWidth=log10fdrw(k), Color="#4B0082")
end

ax5.FontName = opts.fontname;
drawnow;
axis(ax5, "padded")
ax5.YTick = log10fdr;
ax5.YAxis.Label.String = "-log_{10} " + opts.legText;

% fit the Y label
ax5.PositionConstraint = "outerposition";
df = ax5.Position(1) - ax5.OuterPosition(1);
ax5.OuterPosition(1) = ax5.Position(1);
ax5.Position(3) = ax5.Position(3) - df;
ax5.Position(4) = ax5.Position(4)/2;
if max(diff(log10fdr)) > 25
    ax5.YScale = "log";
end
ax5.Visible = "off";
ax5.YAxis.Visible = "on";

% add edges ---------------------------------------------------------------
hold(ax2, "on")

% get boundaries of node objects: this is done to avoid overlapping between
% arrow heads and node margins
nodew = getMarkerBounds(h);

edges = G.Edges;
edges.from = string(edges.EndNodes(:, 1));
edges.to = string(edges.EndNodes(:, 2));
edges.id = edges.to + ":" + edges.from;
tab.id = tab.(opts.traitCol) + ":" + tab.(opts.snpCol);
[~, idx] = ismember(edges.id, tab.id);
tab = tab(idx, :);

% nominally significant (see opts.nominal) associations are shown as dashed
% lines
edges.ls(:) = "-";
if opts.nominal
    edges.ls(~tab.nominal) = "--";
end

% add edge colors based on the effect sizes
if isfield(opts, "effectCol")
    
    edges.effect = tab.(opts.effectCol);

    if opts.effectBinarize
        edges.effect(edges.effect > 0) = 1;
        edges.effect(edges.effect < 0) = -1;
        % cmap = feval(opts.mapEffect);
        try
            cmap = get_colormap(opts.mapEffect, 256);
        catch
            cmap = feval(opts.mapEffect);
        end
        idx = round([0.1, 0.9].*size(cmap, 1)); % keep 2 colors only
        cmap = cmap(idx, :);
        edges.color = nan(height(edges), 3);
        edges.color(edges.effect < 0, :) = repmat(cmap(1, :), sum(edges.effect < 0), 1);
        edges.color(edges.effect > 0, :) = repmat(cmap(2, :), sum(edges.effect > 0), 1);
    else
        cmin = min(edges.effect);
        cmax = max(edges.effect);
        % cmap = feval(opts.mapEffect);
        try
            cmap = get_colormap(opts.mapEffect, 256);
        catch
            cmap = feval(opts.mapEffect);
        end
        idx = fix((edges.effect - cmin)/(cmax - cmin)*length(cmap)) + 1;
        edges.color = squeeze(ind2rgb(idx, cmap));
    end
    
else
    edges.color(:) = "#4B0082";
end

% Loop through the edges and draw a Bezier curve for each edge
for i = 1:height(edges)
    idx_from = ismember(ax.NodeLabel, edges.from(i));
    idx_to = ismember(ax.NodeLabel, edges.to(i));
    x1 = X(idx_from);
    y1 = Y(idx_from);
    x2 = X(idx_to);
    y2 = Y(idx_to);
    
    % Calculate the control points for the Bezier curve
    cx1 = (x1 + x2) / 2;
    cy1 = y1 + opts.d * (y2 - y1);
    cx2 = (x1 + x2) / 2;
    cy2 = y2 + opts.d * (y1 - y2);
    
    % Draw the Bezier curve
    t = linspace(0 , 1, opts.BezierPoints);
    bx = (1-t).^3*x1 + 3*(1-t).^2.*t*cx1 + 3*(1-t).*t.^2*cx2 + t.^3*x2;
    by = (1-t).^3*y1 + 3*(1-t).^2.*t*cy1 + 3*(1-t).*t.^2*cy2 + t.^3*y2;

    % avoid overlapping between arrow heads and nodes: split the curve into
    % two continuous lines
    xn = nodew.x(idx_to, :);
    yn = nodew.y(idx_to, :);

    idxX = bx >= xn(1) & bx <= xn(2);
    idxY = by >= yn(1) & by <= yn(2);
    idx = find(idxX & idxY, 1);
    
    if idx < 3, idx = find(idxX & idxY, 3); idx = idx(end); end

    bx1 = bx(1:idx-1); bx2 = bx(idx-2:end);
    by1 = by(1:idx-1); by2 = by(idx-2:end);

    lax = plot(h(1).Parent, bx1, by1, LineWidth=ax.LineWidth(i), ...
        Color=edges.color(i, :), LineStyle=edges.ls(i));
    lax = line2arrow(lax, 'HeadLength', opts.ArrowHeadLength, 'HeadWidth', opts.ArrowHeadWidth);
    lax.HeadStyle = "cback1";

    plot(h(1).Parent, bx2, by2, LineWidth=ax.LineWidth(i), ...
        Color=edges.color(i, :), LineStyle=edges.ls(i));
    
end

hold(ax2, "off")

if isfield(opts, "effectCol")
    c = colorbar(ax2, 'Location', 'north');
    colormap(ax2, cmap)
    c.Position(1) = ax2.Position(1);
    c.Position(2) = ax2.Position(2) + ax2.Position(4);
    c.Position(4) = c.Position(4)/2;
    c.Position(3) = ax2.Position(3)/10;
    c.AxisLocation  = 'out';
    c.Label.String = "Effect";
    c.Label.Units = 'normalized';
    c.Label.HorizontalAlignment = 'left';
    c.Label.VerticalAlignment = 'middle';
    c.FontSize = opts.fontsizeLegend;
    c.FontName = opts.fontname;
    c.Label.Position = [1 0.5 0];
    if opts.effectBinarize
        c.Ticks = [0.25, 0.75];
        c.TickLabels = ["-", "+"];
    else
        c.Ticks = [0, 1];
        c.TickLabels = round([min(edges.effect), max(edges.effect)], 2, "decimals");
    end
   
end

% send edges to the back and behind nodes
pax = findobj(ax2, 'Tag', "Node");
pax_idx = ismember(ax2.Children, pax);
ax2.Children = [ax2.Children(pax_idx); ax2.Children(~pax_idx)];

% add node labels
text(ax2, X, Y, ax.NodeLabel, FontName=opts.fontname, ...
    VerticalAlignment="bottom", FontSize=opts.fontsize,...
    HorizontalAlignment="center");
ax2.Visible = "off";
close(ax.Parent.Parent)

if opts.save
    exportgraphics(ax2.Parent, opts.output + "." + opts.format, 'Resolution', opts.resolution)
end

if opts.savefig % doesn't work with "scatter"
    savefig(ax2.Parent, opts.output + ".fig")
end

gax = ax2.Parent;
end % END