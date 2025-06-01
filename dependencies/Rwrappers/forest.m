function forest(p, meta, weights, hetero)
% draws forest plot for data coming from meta_analysis
% Oveis Jamialahmadi, GU, Nov 2020.

close all force

% construct default option struct
optsDef = struct('study', '', 'pool', 'both', 'poolRef', 'random', ...
    'studyFont', 'Arial', 'studyFontSize', 13, 'studyBold', false, ...
    'label', '', 'labelFont', 'Arial', 'labelFontSize', 15, ...
    'labelBold', true, 'lineColor', [0 0 0], 'studyMarker', 's', ...
    'studyMarkerColor', [0 0 0], 'poolMarkerColor', 'b', ...
    'showhet', false, 'markerSize', 300, 'savefig', false, ...
    'save', true, 'out', "forestplot", 'res', 600, 'format', 'tiff',...
    'aratio', 1, 'gap', 0.02, 'winstate', 'normal', 'box', 'off', ...
    'axWidth', 1.3, 'lineWidth', 1.5);

if nargin < 1
    help forest
    disp('default fields')
    disp(optsDef)
    return
end
    
% forest plot options 
if ~all(ismember(fieldnames(optsDef), fieldnames(p.forestopt)))
    diffFields = setdiff(fieldnames(optsDef), fieldnames(p.forestopt));
    for i = 1:numel(diffFields)
        p.forestopt.(diffFields{i}) = optsDef.(diffFields{i});
    end
end

if nargin < 4 || isempty(hetero) || ~p.forestopt.showhet
    hetero = [];
end
% if isempty(hetero) || ~p.forestopt.showhet
%     hetero = [];
% end


if p.forestopt.studyBold
    p.forestopt.studyBold = 'bold';
else
    p.forestopt.studyBold = 'normal';
end

if p.forestopt.labelBold
    p.forestopt.labelBold = 'bold';
else
    p.forestopt.labelBold = 'normal';
end

if strcmpi(p.sm, "or")
    p.ES = exp(p.ES);
end

% calculated X-axis limits
XLim = [min([p.ci(:); meta.("CI Lower")]), ...
    max([p.ci(:); meta.("CI Upper")])];
XLim = [0.95*XLim(1), 1.05*XLim(2)];

% calculate Y-axis limist
if strcmpi(p.forestopt.pool, 'both')
    Nticks = size(weights, 1) + 3; % 2 for overall (fixed-/random-) and 1 for gap between.
    labels = ["Random effects", "Fixed effect"];
    p.ci = [p.ci; [meta.("CI Lower"), meta.("CI Upper")]];
    p.ES = [p.ES; meta.Estimate];
elseif any(strcmpi(p.forestopt.pool, {'fixed', 'random'}))
    Nticks = size(weights, 1) + 2; % 1 for overall and 1 for gap between.
    if strcmpi(p.forestopt.pool, 'fixed')
        labels = "Fixed effect";
        idx = contains(lower(meta.Properties.RowNames), 'fixed');
        idxW = contains(lower(weights.Properties.VariableNames), 'fixed');
    else
        labels = "Random effects";
        idx = contains(lower(meta.Properties.RowNames), 'random');
        idxW = contains(lower(weights.Properties.VariableNames), 'random');
    end
    weights = weights(:, idxW);
    p.ci = [p.ci; [meta.("CI Lower")(idx), meta.("CI Upper")(idx)]];
    p.ES = [p.ES; meta.Estimate(idx)];
else
    error('unknown pool option!')
end

% start drawing ===========================================================
if strcmp(p.forestopt.winstate, 'normal')
    fig = figure('visible', 'off', 'WindowState', p.forestopt.winstate);
else
    fig = figure('visible', 'on', 'WindowState', p.forestopt.winstate);
end
ax1 = gca;
ax1.Box = p.forestopt.box;
ax1.YAxis.TickLength = [0, 0];
ax1.XAxis.TickLength = [0.01, 0.01];
if strcmp(ax1.Box, 'off')
    ax1.YAxis.Visible = 'off';
end
ax1.FontSize = p.forestopt.studyFontSize;
ax1.LineWidth = p.forestopt.axWidth;
axis(ax1, 'square')
axis(ax1, 'tight')
ax1.YLim = [.5, Nticks + .5];
ax1.XLim = XLim;    
ax1.YTick = 1:Nticks;
ax1.YTickLabel = [];

%% add study labels -------------------------------------------------------
if isempty(p.forestopt.study)
    p.forestopt.study = [labels, "", flip("Study " + (1:size(weights, 1)))];
else
    if numel(p.forestopt.study) ~= size(weights, 1)
        error('wrong number of study labels!')
    end
    studies = flip(string(p.forestopt.study));
    if iscolumn(studies), studies = studies'; end
    p.forestopt.study = [labels, "", studies];
end

p.forestopt.study = arrayfun(@sprintf, p.forestopt.study); % in case of newline 

newpos = ax1.Position;
newpos(1) = 0.02;
newpos(3) = 0.001;
ax2 = axes('Position', newpos); % for study labels 
ax2.YTick = ax1.YTick;
ax2.YLim = ax1.YLim;

t = text(zeros(numel(ax2.YTick), 1), ...
    (ax2.YTick-ax2.YLim(1))./diff(ax2.YLim), p.forestopt.study, ...
    'HorizontalAlignment', 'left', 'Units', 'normalized', 'FontName',...
    p.forestopt.studyFont, 'FontSize', p.forestopt.studyFontSize,...
    'FontWeight', p.forestopt.studyBold);
poolIdx = find(contains(p.forestopt.study,...
    {'Fixed effect', 'Random effects'}));
for i = 1:numel(poolIdx)
    t(poolIdx(i)).FontWeight = 'bold';
end

if p.forestopt.showhet % show heterogeneity info
    hetero = convertvars(hetero, 1:width(hetero), @(x)double(string(x)));
    if hetero.tau2 == 0
        heterotau2 = '%.2g, ';
    else
        heterotau2 = '%.2f, ';
    end
    p.forestopt.het = "Heterogenity: {\it I}^2 = " + ...
        compose('%.2g%%, ', hetero.I2.*100) + "{\tau}^{2} = " + ...
        compose(heterotau2, hetero.tau2) + "{\it p} = " + ...
        compose('%.2f', hetero.("pval.Q"));
end

[~, maxLabelLen] = max(arrayfun(@(x)x.Extent(3), t));

while(abs(t(maxLabelLen).Extent(3) - 1) > 1e-6) % fit ax2 to contain study labels
    ax2.Position(3) = ax2.Position(3)*t(maxLabelLen).Extent(3);
    pause(0.2);
end
ax2.Visible = 'off';

% move ax1 to left
newpos(1) = ax2.Position(1) + ax2.Position(3); % new start position for ax1
newpos(2) = ax1.Position(1) - newpos(1) + ax1.Position(3);
ax1.Position(1) = newpos(1); %ax1.Position(3) = newpos(2);
pos = plotboxpos(ax1);% get real position of ax1 (since we squared it)
ax1.Position(1) = 2.*ax1.Position(1) - pos(1);

%% plot effect sizes ------------------------------------------------------
ax1.FontWeight = p.forestopt.labelBold;
ax1.LineWidth = 1.3;

if isempty(p.forestopt.label) % on X-axis
    if strcmpi(p.sm, "or")
        p.forestopt.label = "Odds ratio";
        ax1.XScale = 'log';
        ax1.XMinorTick = 'off';
    else
        p.forestopt.label = "Effect size";
    end
end

if p.forestopt.showhet
    title(ax1, p.forestopt.label, 'FontWeight', p.forestopt.labelBold, ...
        'FontSize', p.forestopt.labelFontSize, ...
        'FontName', p.forestopt.labelFont);
    xlabel(ax1, p.forestopt.het, ...
        'Units', 'normalized', 'FontName',...
        p.forestopt.studyFont, 'FontSize', p.forestopt.studyFontSize,...
        'FontWeight', p.forestopt.studyBold);
else
    xlabel(ax1, p.forestopt.label, 'FontWeight', p.forestopt.labelBold, ...
        'FontSize', p.forestopt.labelFontSize, ...
        'FontName', p.forestopt.labelFont);
end

if strcmpi(p.sm, "or") % OR
    line(ax1, [1, 1], ax1.YLim, 'LineStyle', '-', 'LineWidth', ...
        p.forestopt.lineWidth, 'Color', [128 128 128]./256)
else % beta
    line(ax1, [0, 0], ax1.YLim, 'LineStyle', '-', 'LineWidth', ...
        p.forestopt.lineWidth, 'Color', [128 128 128]./256)
end

if strcmpi(p.forestopt.pool, 'both') % construct a vecotr with CI coordinates
    ciCoor = flip([1, 2, 4:ax1.YTick(end)]);
    [poolFlag, poolLine] = deal(flip([true true, false(1, numel(ciCoor) - 2)]));
%     if strcmpi(p.forestopt.poolRef, 'fixed')
%         poolLine(end) = false;
%     else
%         poolLine(end - 1) = false;
%     end
else
    ciCoor = flip([1, 3:ax1.YTick(end)]);
    [poolFlag, poolLine] = deal(flip([true, false(1, numel(ciCoor) - 1)]));
end

% set study marker sizes based on pooled weights
if size(weights, 2) > 1 
    idx = ismember(lower(weights.Properties.VariableNames), ...
        p.forestopt.poolRef);
    W = weights.(find(idx));
else
    W = weights{:, :};
end
% markerSize = ( (W - min(W)) .* abs(diff(p.forestopt.markerSize)) )./ ...
%     (max(W) - min(W)) + p.forestopt.markerSize(1);
markerSize = (W) .* p.forestopt.markerSize;

hold(ax1, 'on')
for ii = 1:numel(ciCoor)
    
    if poolFlag(ii) 
        patch(ax1, 'Vertices', [p.ES(ii) ciCoor(ii)+.12; ...
            p.ci(ii, 1) ciCoor(ii); p.ES(ii) ciCoor(ii)-.12; ...
            p.ci(ii, 2) ciCoor(ii)], 'Faces', 1:4,...
            'FaceColor', p.forestopt.poolMarkerColor, ...
            'EdgeColor', p.forestopt.lineColor, 'FaceAlpha', 0.9);
        
        if poolLine(ii) % vertical line passing pooled effect
            line(ax1, [p.ES(ii), p.ES(ii)], ax1.YLim, ...
                'LineStyle', '--', ...
                'LineWidth', 0.8*p.forestopt.lineWidth,...
                'Color', [128 128 128]./256)
        end
        
    else
        line(ax1, p.ci(ii, :),[ciCoor(ii), ciCoor(ii)], 'LineStyle', '-',...
            'LineWidth', p.forestopt.lineWidth, 'Color', p.forestopt.lineColor)
    
        line(ax1, [p.ci(ii, 1), p.ci(ii, 1)], ...
            [ciCoor(ii)-.05, ciCoor(ii)+.05], 'LineStyle', '-', ...
            'LineWidth', p.forestopt.lineWidth, 'Color', p.forestopt.lineColor)

        line(ax1, [p.ci(ii, 2), p.ci(ii, 2)], ...
            [ciCoor(ii)-.05, ciCoor(ii)+.05], 'LineStyle', '-', ...
            'LineWidth', p.forestopt.lineWidth, 'Color', p.forestopt.lineColor)
        
        scatter(ax1, p.ES(ii) ,ciCoor(ii), ....
            'Marker', p.forestopt.studyMarker, ...
            'MarkerEdgeColor', p.forestopt.lineColor, ...
            'MarkerFaceColor', p.forestopt.studyMarkerColor, ...
            'MarkerFaceAlpha', 0.9, ....
            'SizeData', markerSize(ii))
    end
    
end
hold(ax1, 'off')

%% add CI labels -----------------------------------------------------------
Ytitle = 1.05; % Title Y-axis position
% add ax3 to right or ax1
newpos = plotboxpos(ax1);
newpos(1) = newpos(1) + newpos(3); % new start position for ax1
newpos(3) = 1 - ax2.Position(3) - ax1.Position(3);
if newpos(3) < 0
    newpos(3) = 0.1;
end
ax3 = axes('Position', newpos);
ax3.YTick = ax1.YTick;
ax3.YLim = ax1.YLim;
ax3.YTickLabel = [];

% create effect-[CI] labels
p.ES_ci = strings(numel(p.forestopt.study), 1);
p.ES_ci(flip(p.forestopt.study) ~= "") = compose("%.2f [%.2f, %.2f]", p.ES, p.ci);

if strcmpi(p.sm, "or") % OR
    p.ES_ci_title = compose('Odds ratio\n[95% CI]');
else % beta
    p.ES_ci_title = compose('Effect\n[95% CI]');
end

% print effect-[CI] labels
projectLabels(ax3, zeros(numel(ax3.YTick), 1), ...
    flip((ax3.YTick-ax3.YLim(1))./diff(ax3.YLim)), p, p.ES_ci, Ytitle,...
    p.ES_ci_title);

%% add wieghts labels -----------------------------------------------------

% add ax4 to right or ax3
newpos = ax3.Position;
newpos(1) = ax3.Position(1) + ax3.Position(3); % new start position for ax1
newpos(3) = ax3.Position(3);
ax4 = axes('Position', newpos); % for study labels 
ax4.YTick = ax3.YTick;
ax4.YLim = ax3.YLim;
ax4.YTickLabel = [];

if strcmpi(p.forestopt.pool, 'both') % add both fixed-/random-weights
    p.weights_title = {sprintf('Weight\n(fixed)'), sprintf('Weight\n(random)')};
    p.weights.fixed = compose('%.1f%%', weights.fixed.*100);
    p.weights.random = compose('%.1f%%', weights.random.*100);
else
    p.weights_title = sprintf('Weight\n(%%)');
    if strcmpi(p.forestopt.pool, 'fixed')
        p.weights.overall = compose('%.1f%%', weights.fixed.*100);
    else
        p.weights.overall = compose('%.1f%%', weights.random.*100);
    end
end

% print weight labels
weightLabelX = zeros(size(weights, 1), 1);
weightLabelY = (ax4.YTick(end - size(weights, 1) + 1:end)...
    - ax4.YLim(1))./diff(ax4.YLim);

if strcmpi(p.forestopt.pool, 'both') % print both fixed-/random-weights   
    
     projectLabels(ax4, weightLabelX, weightLabelY, p, ...
         flip(p.weights.fixed), Ytitle, p.weights_title{1});
    
    % add ax5 to right or ax4
    newpos = ax4.Position;
    newpos(1) = ax4.Position(1) + ax4.Position(3); % new start position for ax5
    newpos(3) = ax4.Position(3);
    ax5 = axes('Position', newpos); % for study labels 
    ax5.YTick = ax4.YTick;
    ax5.YLim = ax4.YLim;
    ax5.YTickLabel = [];
    
    projectLabels(ax5, weightLabelX, weightLabelY, p, ...
         flip(p.weights.random), Ytitle, p.weights_title{2});

else % either fixed- or random-effects weights
    
    projectLabels(ax4, weightLabelX, weightLabelY, p, ...
         flip(p.weights.overall), Ytitle, p.weights_title);
    
end

%% adjust axes width ------------------------------------------------------

% this is the figure width containing text labels/info etc. So, this width
% is fixed. keepWidth may go above 1 (e.g. with large font size), but we
% assumed it will remain under 1 most of the time (with default settings).

keepWidth = ax2.Position(3) + ax3.Position(3) + ax4.Position(3) + ...
    2.*p.forestopt.gap;
if exist('ax5', 'var')
    keepWidth = keepWidth + ax5.Position(3) + p.forestopt.gap;
end

ax1Pos = plotboxpos(ax1);
ax5Modify = false;
if (ax1Pos(3) + keepWidth) > 1 % if does not fit in current fig
    ax5Modify = true;
    
    % modify ax1 width
    ax1.Position(3) = 1 - keepWidth;
    
    % check again if ax1 overlaps ax2
    ax1.Position(1) = ax2.Position(1) + ax2.Position(3);
    pos = plotboxpos(ax1);% get real position of ax1 (since we squared it)
    ax1.Position(1) = 2.*ax1.Position(1) - pos(1);
    
    % slide ax3, ax4 or ax5 to left
    newX = ax1.Position(3) + ax2.Position(3) + p.forestopt.gap;
    ax3.Position(1) = newX;
    newX = newX + ax3.Position(3) + p.forestopt.gap;
    ax4.Position(1) = newX;
    
    % modify Y-pos and height of ax2, 3, 4 and 5 (because we squared ax1)
    ax1Pos = plotboxpos(ax1);
    ax2.Position(2) = ax1Pos(2); ax2.Position(4) = ax1Pos(4);
    ax3.Position(2) = ax1Pos(2); ax3.Position(4) = ax1Pos(4);
    ax4.Position(2) = ax1Pos(2); ax4.Position(4) = ax1Pos(4);
end
ax3.Visible = 'off'; ax4.Visible = 'off';

if exist('ax5', 'var')
    
    if ax5Modify
        ax5.Position(1) = newX + ax4.Position(3) + p.forestopt.gap;
        ax5.Position(2) = ax1Pos(2); ax5.Position(4) = ax1Pos(4);
    end
    ax5.Visible = 'off';
    ax5.PlotBoxAspectRatio(2) = ax5.PlotBoxAspectRatio(2)/p.forestopt.aratio;
end

% not recommended when ax1 sclae is set to square
ax1.PlotBoxAspectRatio(2) = ax1.PlotBoxAspectRatio(2)/p.forestopt.aratio;
ax2.PlotBoxAspectRatio(2) = ax2.PlotBoxAspectRatio(2)/p.forestopt.aratio;
ax3.PlotBoxAspectRatio(2) = ax3.PlotBoxAspectRatio(2)/p.forestopt.aratio;
ax4.PlotBoxAspectRatio(2) = ax4.PlotBoxAspectRatio(2)/p.forestopt.aratio;

fig.Color = [1 1 1];
fig.Visible = 'on';
%% save forest plot -------------------------------------------------------
if p.forestopt.savefig
    savefig(p.forestopt.out + ".fig")
end

if p.forestopt.save
    print(fig, p.forestopt.out, "-d" + p.forestopt.format,...
        "-r" + p.forestopt.res);
end


end % END

%% subfunctions ===========================================================
function textHandle = projectLabels(ax, X, Y, p, vals, Ytitle, Title)
textHandle = text(ax, X, Y, vals, ...
    'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontName',...
    p.forestopt.studyFont, 'FontSize', p.forestopt.studyFontSize,...
    'FontWeight', p.forestopt.studyBold);

% make overall (fixed-/random-effects) values bold
idx = find(vals == ""):numel(vals);
if ~isempty(idx)
    idx(1) = [];
    for i = idx
        textHandle(i).FontWeight = 'bold';
    end
end

% print title
textHandle(end + 1) = text(ax, ...
    max(arrayfun(@(x)x.Extent(3), textHandle))/2,...
    Ytitle, Title, ...
    'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontName',...
    p.forestopt.studyFont, 'FontSize', p.forestopt.studyFontSize + 0.5,...
    'FontWeight', 'bold');

[~, maxLabelLen] = max(arrayfun(@(x)x.Extent(3), textHandle));
while(abs(textHandle(maxLabelLen).Extent(3) - 1) > 1e-6) % fit ax3 to contain effect-[CI] labels
    ax.Position(3) = ax.Position(3)*textHandle(maxLabelLen).Extent(3);
end

for i = 1:numel(textHandle)
    textHandle(i).Position(1) = 0.5;
end

end % END

%%
function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h
% Copyright 2010 Kelly Kearney
% Check input
if nargin < 1
    h = gca;
end
if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end
% Get position of axis in pixels
currunit = get(h, 'units');
axisPos  = getpixelposition(h);
% Calculate box position based axis limits and aspect ratios
darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');
if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else
    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis limits auto-set via Inf/-Inf use
    
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end
    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');
    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);
    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end
% Convert plot box position to the units used by the axis
hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');
temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);
set(hfig, 'currentaxes', currax);
end