function gNetworkPlot(tab, opts)
% visualizes LDSC cell-type specific/group enrichment results. Inspired by
% https://www.nature.com/articles/s41588-018-0047-6.
% Oveis Jamialahmadi, University of Gothenburg, 27 Februrary 2023.

arguments
    tab {mustBeA(tab, "table")} % results of cell-type specific analysis from LDSC. Example table can be found within the same directory: ldsc_cell_type_specific.mat
    opts.pCol {mustBeTextScalar} = "Coefficient P" % p-values column
    opts.cellCol {mustBeTextScalar} = "Cell-type" % cell types column
    opts.traitCol {mustBeTextScalar} = "Trait" % traits column
    opts.traitCatCol {mustBeTextScalar} = "Trait category" % trait categories column
    opts.tissueCol {mustBeTextScalar} = "Tissue" % tissue type column
    opts.fdr = 0.05 % FDR-CUTOFF
    opts.edgeFactor = 3 % max width of edges
    opts.nodeFactor = 9 % adds this value to node size
    opts.maxNodeSize = 17 % max node size: if any node size (upon the addition of 'nodeFactor') > this value, node sizes will be capped by this value.
    opts.mapCell (1,1) {mustBeMember(opts.mapCell, ["jet", "turbo", "parula", "hot", "winter", "copper", "bone"])} = "turbo"
    opts.mapTrait(1,1) {mustBeMember(opts.mapTrait, ["jet", "turbo", "parula", "hot", "winter", "copper", "bone"])} = "jet"
    opts.fontname {mustBeTextScalar} = "Garamond" % node font name
    opts.fontsize (1,1) double = 11 % node font size
    opts.fontsizeLegend (1,1) double = 9 % legend's font size
    opts.layout {mustBeMember(opts.layout, ["force", "circle", "layered"])} = "force"
    opts.BezierPoints (1,1) double = 60 % number of points used for Bezier curve
    opts.d (1,1) double = 0.2 % Bezier control point distance factor

    opts.output {mustBeTextScalar} = "CellNetworkPlot"
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "jpg" % format of output plots
    opts.resolution (1, 1) double = 600 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well

end

% tab = load("ldsc_cell_type_specific.mat").res;

% filter based on 'fdr' cutoff
tab.q = mafdr(tab.(opts.pCol), "BHFDR", true);
tab(tab.q >= opts.fdr, :) = [];

if isempty(tab)
    disp('no significant results remained after adjusting for FDR!')
end

% tab.q = tab.("Coefficient P").^10; 

% correct cell-tyepe names
tab.(opts.cellCol) = replace(tab.(opts.cellCol), "_", " ");

% trait and cell-type counts 
cnt.trait = groupsummary(tab, opts.traitCol);
cnt.ct = groupsummary(tab, opts.cellCol);

% add colors
utissues = unique(tab.(opts.tissueCol), "stable");
col_tissue = feval(opts.mapCell, numel(utissues));
[~, idx] = ismember(tab.(opts.tissueCol), utissues);
tab.col_tissue = col_tissue(idx, :);

ucats = unique(tab.(opts.traitCatCol), "stable");
col_traits = feval(opts.mapTrait, numel(ucats));
[~, idx] = ismember(tab.(opts.traitCatCol), ucats);
tab.col_trait = col_traits(idx, :);

% create a digraph
G = digraph(tab.(opts.cellCol), tab.(opts.traitCol), -log10(tab.q));

% for each trait-cell pair, keep only the strongest association per each marker 
G = simplify(G, 'max');

% edge withs correpond to -log10 FDR
LWidths = opts.edgeFactor*G.Edges.Weight/max(G.Edges.Weight);

% plot the graph
close all force
ax = plot(G, LineWidth=LWidths, Layout=opts.layout);

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

[f1, f2] = ismember(ax.NodeLabel, cnt.ct.(opts.cellCol)); 
ax.MarkerSize(f1) = cnt.ct.GroupCount(f2(f1));

ax.MarkerSize = ax.MarkerSize + opts.nodeFactor;

% cap node sizes
if any(ax.MarkerSize) > opts.maxNodeSize
    ax.MarkerSize = rescale(ax.MarkerSize, min(ax.MarkerSize), opts.maxNodeSize);
end

% node colors
[f1, f2] = ismember(ax.NodeLabel, tab.(opts.traitCol));
ax.NodeColor = repmat(tab.col_tissue(1, :), numel(f1), 1);
ax.NodeColor(f1, :) = tab.col_trait(f2(f1), :);

[f1, f2] = ismember(ax.NodeLabel, tab.(opts.cellCol));
ax.NodeColor(f1, :) = tab.col_tissue(f2(f1), :);

%% generate a new graph plot for better control over the figure behaviors--
fig = figure(WindowState="maximized");
h = axes(fig);
hold(h, "on")
X = ax.XData;
Y = ax.YData;

h = gscatter(h, X, Y, 1:numel(X));
ax2 = h(1).Parent;
ax2.Position(1) = 0.075; % move the axis to the left, to make extra space for legends to the right
ax2.Legend.Visible = "off";

for k = 1:numel(h)
    h(k).Marker = ax.Marker(k);
    h(k).MarkerSize = ax.MarkerSize(k);
    h(k).MarkerFaceColor = ax.NodeColor(k, :);
    h(k).MarkerEdgeColor = 'k';
    h(k).Tag = "Node";
end
hold(ax2, "off")

% add legends for nodes: traits and cell-type groups ----------------------
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
leg.Title.String = "Trait category";
hold(ax2, "off")

tmp = unique(tab(:, [opts.tissueCol, "col_tissue"]), "rows");
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
    lg2(k) = plot(ax3, NaN, NaN, MarkerFaceColor=tmp.col_tissue(k,:),...
            LineStyle='none', Marker='s', MarkerEdgeColor='k',...
            Color=tmp.col_tissue(k, :), MarkerSize=mean(ax.MarkerSize));
end

leg2 = legend(lg2, tmp.(1), Color=[0.93 0.93 0.93], AutoUpdate='off',...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    Location='southeastoutside', Orientation='Vertical');
leg2.Title.String = "Cell-type group";
hold(ax3, "off")

drawnow
ax3.Position = ax2.Position;

% shift the second legend upwards
drawnow;
leg2.Position(2) = leg.Position(2) - leg2.Position(4) - 0.01;

% legend for marker size
nodeinfo.marker = string(ax.Marker);
nodeinfo.size = ax.MarkerSize;

% create a new axis for trait/cell-type marker sizes ----------------------
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

% add two annotation text boxes for traits and cell-types
% nodew = getMarkerBounds(h); % non-overlapping (horizontally) markers for legend
drawnow;
axtest = plot(ax4, 1, 0.5, Marker="O", MarkerSize=max(nodeinfo.size));
te =  text(ax4, 1, 0.5, "Enriched cell-types (per trait)", FontSize=opts.fontsizeLegend, FontName=opts.fontname); % test text title
teExt = te.Extent; % how much space do we need for text title?
teExt = 1.05*teExt(4); % the space we need for the title string
nodew = getMarkerBounds(ax4.Children);
delete(axtest)
delete(te)
nodew.w = nodew.w.*.5;
nodew.h = nodew.h.*.5;

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
text(ax4, midX, 2.2*maxMarginY+1.5*teExt, "Enriched cell-types (per trait)", ...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    VerticalAlignment="middle", HorizontalAlignment="center", FontWeight="bold");

rectangle(ax4, Position=pos2);
yline(ax4, maxMarginY);
text(ax4, midX, maxMarginY+0.5*teExt, "Traits (per cell-type)", ...
    FontSize=opts.fontsizeLegend, FontName=opts.fontname, ...
    VerticalAlignment="middle", HorizontalAlignment="center", FontWeight="bold");

idx = nodeinfo.marker == "o";
msize1 = [max(nodeinfo.size(idx)), median(nodeinfo.size(idx)), min(nodeinfo.size(idx))];
msize2 = [max(nodeinfo.size(~idx)), median(nodeinfo.size(~idx)), min(nodeinfo.size(~idx))];
mlabel1 = round([max(cnt.trait.GroupCount), median(cnt.trait.GroupCount), min(cnt.trait.GroupCount)]);
mlabel2 = round([max(cnt.ct.GroupCount), median(cnt.ct.GroupCount), min(cnt.ct.GroupCount)]);

hold(ax4, "on")
for k = 1:numel(msize1)
    plot(ax4, 1 + maxMarginX*(k-1), teExt+1.7*maxMarginY, Marker="O", MarkerSize=msize1(k), Color="k");
    text(ax4, 1 + maxMarginX*(k-1), teExt+1.7*maxMarginY, string(mlabel1(k)), ...
        VerticalAlignment="middle", HorizontalAlignment="center", ...
        FontName=opts.fontname, FontSize=opts.fontsizeLegend-2);

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
ax5.Position(4) = ax4.Position(2) - ax5.Position(2) - 0.015;
df = ax5.Position(4) - ax5.Position(3);
if df > 0
    ax5.Position(4) = ax5.Position(3);
    ax5.Position(2) = ax5.Position(2) + df;
end
log10fdr = -log10(tab.q);
log10fdr = ceil([min(log10fdr), median(log10fdr), max(log10fdr)]);
[log10fdr, uidx] = unique(log10fdr, "stable");
log10fdrw = [min(ax.LineWidth), median(ax.LineWidth), max(ax.LineWidth)];
log10fdrw = log10fdrw(uidx);

for k = 1:numel(log10fdr)
    yline(ax5, log10fdr(k), LineWidth=log10fdrw(k), Color="#4B0082")
end

ax5.FontName = opts.fontname;
drawnow;
axis(ax5, "padded")
ax5.YTick = log10fdr;
ax5.YAxis.Label.String = "-log_{10} FDR";

% fit the Y label
ax5.PositionConstraint = "outerposition";
df = ax5.Position(1) - ax5.OuterPosition(1);
ax5.OuterPosition(1) = ax5.Position(1);
ax5.Position(3) = ax5.Position(3) - df;
ax5.Position(4) = ax5.Position(4)/2;
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

    bx1 = bx(1:idx-1); bx2 = bx(idx-2:end);
    by1 = by(1:idx-1); by2 = by(idx-2:end);

    lax = plot(h(1).Parent, bx1, by1, LineWidth=ax.LineWidth(i), Color="#4B0082");
    lax = line2arrow(lax);
    lax.HeadStyle = "cback1";
    plot(h(1).Parent, bx2, by2, LineWidth=ax.LineWidth(i), Color="#4B0082");
    
end

hold(ax2, "off")

% send edges to the back and behind nodes
pax = findobj(ax2, 'Tag', "Node");
pax_idx = ismember(ax2.Children, pax);
ax2.Children = [ax2.Children(pax_idx); ax2.Children(~pax_idx)];

% add node labels
text(ax2, X, Y, ax.NodeLabel, FontName=opts.fontname, ...
    VerticalAlignment="bottom", FontSize=opts.fontsize);
ax2.Visible = "off";
close(ax.Parent.Parent)

if opts.save
    exportgraphics(ax2.Parent, opts.output + "." + opts.format, 'Resolution', opts.resolution)
end

if opts.savefig % doesn't work with "scatter"
    savefig(ax2.Parent, opts.output + ".fig", 'compact')
end

end % END