function barPlotter(res, opts)
% draws bar plots for stratification analyses in 'gwasrunner' function.
% INPUTs
%   res: a struct with descriptive statistics for variant and
%        stratification trait under test.
%   opts: option structs of gwasrunner (only few fields are used here).
% 
% Note: this function is still under development, so all aesthetics are
%       hard coded for the time being.
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, Sep 2021.
% 
% @10/12/2021 some bugs were fixed.
% @11/01/2022 'adjLegendPos' and 'adjNfont' subfunctions have been added.

arguments
    res {mustBeA(res, 'struct')}
    opts.barfunc {mustBeMember(opts.barfunc, ["median", "mean"])} = "mean"
    opts.stratify {mustBeText, mustBeVector} = "strata"
    opts.output {mustBeTextScalar} = "barplotter"
    opts.saveplot (1,1) logical = true
    opts.savefig (1,1) logical = false
    opts.format  {mustBeTextScalar} = "png"
    opts.resolution (1,1) double {mustBeGreaterThan(opts.resolution, 96)} = 300;
    opts.fontsize (1,1) double {mustBeGreaterThan(opts.fontsize, 0)} = 9;
end

n = size(res.Pheno, 1);

for i = 1:n % draw bar plots for every variant
    
    in = struct; % input struct for the bar plot
    % prepare labels (genotypes) 
    in.geno = squeeze(res.usnp(i, 1, :));
    if ~isempty(res.A1(i)) && res.A1(i) ~= "" % if maj/min alleles are present, use them
        in.geno = [res.A1(i)+res.A1(i), res.A1(i)+res.A2(i), res.A2(i)+res.A2(i)];
    end
    in.snp = res.SNP(i);
    
    % phenotype values on y-axis 
    if opts.barfunc == "mean"
        try
            in.y = squeeze(res.mean_s(i, :, :));
            in.ysd = squeeze(res.sd_s(i, :, :)); % shown as error bars
        catch % for binary traits
            in.y = squeeze(res.n_perc_s(i, :, :));
            in.ysd = nan(size(in.y));
        end
    else
        try
            in.y = squeeze(res.median_s(i, :, :));
            in.ysd = squeeze(res.iqr_s(i, :, :)); % shown on each bar
            in.ysd = round(in.ysd, 1, 'decimals');
        catch % for binary traits
            in.y = squeeze(res.n_perc_s(i, :, :));
            in.ysd = nan(size(in.y));
        end
    end
    in.func = opts.barfunc;
    
    try % numbers per each stratum
        in.n = squeeze(res.n_s(i, :, :));
    catch % numbers of cases per each stratum
        in.n = squeeze(res.n_case_s(i, :, :));
    end
    
    % x-/y-axis labels
    in.ylab = res.Pheno(i); % + " (" + opts.barfunc + ")";
    in.ylab = replace(in.ylab, ", n(%)", " (%)"); % for binary traits
    in.xlab = opts.stratify;
    in.xtick = replace(res.Strata(i, :), opts.stratify, '');
    try % compatible with gwasrunner input
        tmp = cellfun(@(x)round(double(string(x)), 1), regexp(in.xtick, '\d+\.?\d*', 'match'), 'uni', false);
        if any(contains(in.xtick, ["<", "≥", ">"]))
            in.xtick(1) = "<" + string(tmp(1));
            in.xtick(end) = "≥" + string(tmp(end));
            for j = 2:numel(in.xtick)-1
                in.xtick(j) = join(string(tmp{j}), '-');
            end
        else % a trait like sex with 0/1 but no relational operators
            for j = 1:numel(in.xtick)
                in.xtick(j) = join(string(tmp{j}), '-');
            end

            % case(0) control(1)
            if numel(in.xtick) < 3 && all(ismember(in.xtick, ["0", "1"]))
                in.xtick = replace(in.xtick, ["0", "1"], ["control", "case"]);
            end
        end
        % in.xtick(2:end-1) = strip(replace(in.xtick(2:end-1), ["<", ">", "≤", "≥"], ' '));
        % in.xtick(2:end-1) = replace(in.xtick(2:end-1), '  ', '-');
    catch % may not have the same format as defined in gwasrunner (custom inputs)
    end
    
    % add p-values (if present)
    try
        in.p = res.padj;
    catch
        in.p = nan;
    end
    
    in.sfontsize = opts.fontsize;
    
    h = drawBarPlot(in);
    
    % apply aesthetics
    % h(1).Parent.Parent.WindowState = 'maximized';
    h(1).Parent.LineWidth = 1.5;
    h(1).Parent.YAxis.FontSize = opts.fontsize + 6;
    h(1).Parent.YAxis.Label.FontSize = opts.fontsize + 12;
    drawnow;
    h(1).Parent.OuterPosition(1) = 0; % fix to zero (may vary with font size change)
    adjFontSize(h(1).Parent.YAxis.Label); % check if y-axis label is too large
    h = adjLegendPos(h); % adjust legend position
    h = adjNfont(h); % adjust N font size to avoid overlaps
    
    savename = matlab.lang.makeValidName(opts.output) + "." + ...
        matlab.lang.makeValidName(res.Pheno(i)) + "." + ...
        matlab.lang.makeValidName(res.SNP(i)) + ".bar.";
    if opts.saveplot
        try
            exportgraphics(h(1).Parent, savename + opts.format, ...
                'Resolution', opts.resolution)
        catch % invalid name
            savename = matlab.lang.makeValidName(savename) + ".";
            savename = regexprep(replace(savename, '_', '.'), '\.{2,}', '.');
            exportgraphics(h(1).Parent, savename + opts.format, ...
                'Resolution', opts.resolution)
        end
    end
    
    if opts.savefig
        try
            savefig(h(1).Parent.Parent, savename + "fig")
        catch % invalid name
            savename = matlab.lang.makeValidName(savename) + ".";
            savename = regexprep(replace(savename, '_', '.'), '\.{2,}', '.');
            savefig(h(1).Parent.Parent, savename + "fig")
        end
    end
    close(h(1).Parent.Parent)
end

end % END

%% subfunctions ===========================================================
function hBar = drawBarPlot(in)

close all force 
f = figure;%('WindowState', 'maximized');
if isvector(in.y) && isrow(in.y) % overall bar graph (non-stratified per snp)
    in.y = in.y.';
    try in.n = in.n.'; catch; end
end
ax = axes(f);
hBar = bar(ax, in.y);

% add error bars
if in.func == "mean"
    hold on
    errorbar([hBar.XEndPoints], in.y(:), [], in.ysd(:),...
        'k', 'linestyle', 'none', 'LineWidth', 1.5)
    hold off
    ysd = in.ysd;
else
    ysd = zeros(size(in.ysd));
end

hBar(1).Parent.TickLength = [0 0];

% increase YLim
if in.func == "mean" && ~any(isnan(ysd), 'all')
    hBar(1).Parent.YLim(2) = 1.15*max(in.y(:) + ysd(:)); 
else
    hBar(1).Parent.YLim(2) = 1.065*max(in.y(:));
end

% Labeling ---------------------------------------------------------------- 
% Insert p-values
if isvector(ysd) && isrow(ysd) % overall bar graph (non-stratified per snp)
    ysd = ysd.';
end
[ctr, ydt] = deal(zeros(1, size(ysd, 1)));
if in.func == "mean"
    if all(~isnan(in.p))
        text(1:size(in.y, 1), 1.04*max(in.y + ysd, [], 2), ...
            strcat('P= ', num2str(in.p', '%0.2g')), ...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom', 'Tag', "p")
    end
    
else
    if all(~isnan(in.p))
        text(1:size(N_perc,1), max(N_perc,[],2), ...
            strcat('P= ',num2str(in.p','%0.2g')), ...
            'HorizontalAlignment', 'center', ....
            'VerticalAlignment', 'bottom', 'Tag', "p")
    end
    
    % show IQR values on each bar
    for k = 1:size(in.y, 2)
        ctr(k,:) = bsxfun(@plus, hBar(k).XData, hBar(k).XOffset'); 
        ydt(k,:) = hBar(k).YData + ysd(:, k).';
        text(ctr(k,:), ydt(k,:), string(in.ysd(:, k)),...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', in.sfontsize, ...
            'Tag', "iqr")
    end
    clear ctr ydt
end

% shift current axes upwards to make some space for lables
hBar(1).Parent.Position(2) = 0.95 - hBar(1).Parent.Position(4);

% Label X-axis ------------------------------------------------------------
% add numbers per each stratum
nX = hBar(1).Parent.XLim(1);
nY = -hBar(1).Parent.Position(2)/6.*hBar(1).Parent.YLim(2);
text(nX, nY, 'N =','HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontSize', in.sfontsize+3, ...
    'Color','k', 'Tag', "n");

for k = 1:size(in.y, 2)
    ctr(k,:) = bsxfun(@plus, hBar(1).XData, hBar(k).XOffset'); 
    ydt(k,:) = repmat(nY, 1, numel(hBar(k).YData));
     text(ctr(k,:), ydt(k,:), sprintfc('%.0f', in.n(:, k)),...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', in.sfontsize, ...
        'Color', 'k', 'Tag', "n");
end

% add stratum labels
% hBar(1).Parent.XTickLabel = in.xtick;
% text(hBar(1).Parent.XLim(1), 0, ...
%     in.xlab, 'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'top', 'FontSize', in.sfontsize+2, 'Color','k')

hBar(1).Parent.XTickLabel = '';
nY = nY*3;
text(hBar(1).Parent.XLim(1), nY, ...
    in.xlab, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'cap', 'FontSize', in.sfontsize+5, 'Color','k',...
    'Interpreter', 'none', 'Tag', "strata")

text(1:size(in.y ,1), repmat(nY, 1, numel(hBar(1).YData)), in.xtick,...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'cap', 'FontSize', in.sfontsize+5, 'Tag', "strata")

% % shift also x tick labels upwards
% % thanks to https://undocumentedmatlab.com/articles/customizing-histogram-plots
% hBar(1).Parent.XRuler.TickLabelGapOffset = -2;

% Color
CfaceCodes = [0 0 0;128 128 128;256 256 256];
if isscalar(hBar), CfaceCodes = 128.*ones(1, 3); end
for ii = 1:numel(hBar)
    hBar(ii).FaceColor = CfaceCodes(ii,:)./256;
    hBar(ii).LineWidth = 1;
end

%Legend
if ~any(ismissing(in.geno) | in.geno == "") % don't show the legend if nothing's to show!
    h = legend(hBar, in.geno, 'Orientation', 'Horizontal', 'Box', 'on', ...
        'FontSize', in.sfontsize+4);
    h.Title.String = in.snp;
    h.Title.Interpreter = 'none';
    set(h,'Location','northwest')
    % h.Position = [1.037*h.Position(1) 1.024*h.Position(2) h.Position(3:4)];
end

hBar(1).Parent.YLabel.String = in.ylab;

end

%% ========================================================================
function h = adjFontSize(h)
% checks the font size of axis h, if too large, reduced the font size
drawnow;
h.Units = 'normalized';
while h.Extent(4) > 1
    h.FontSize = h.FontSize - 0.5;
    drawnow;
end
drawnow;
end 

%% ========================================================================
function h = adjLegendPos(h)
% checks if legend overalps data and moves it upwards. This means the
% position should be north, northeast, or northwest
drawnow;
leg = h(1).Parent.Legend;
if isempty(leg); return; end
pos = leg.Position;
loc = leg.Location;
ylim = h(1).Parent.YLim;
xlim = h(1).Parent.XLim;

% check which bars are covered by the legend's extension
xoverlap = ([h.XEndPoints] - xlim(1))./diff(xlim) <= (pos(1) + pos(3));

% check max length of bars
errEnds = ([h.YEndPoints] - ylim(1))./diff(ylim); % normalized
errEnds = max(errEnds(xoverlap)); % tallest bar within legend area

% check iqr text objects (for median)
iqr = findobj(h(1).Parent, 'Tag', 'iqr');
if ~isempty(iqr)
    xiqr = [iqr.Position]; xiqr = xiqr(1:3:end);
    xoverlapiqr = (xiqr - xlim(1))./diff(xlim) <= (pos(1) + pos(3));
    errEnds = reshape([iqr.Extent], 4, []).';
    errEnds = errEnds(:, 2) + errEnds(:, 4); % y-point + height
    errEnds = (errEnds - ylim(1))./diff(ylim); % normalized
    errEnds = max(errEnds(xoverlapiqr));
end

% similarly, check the error bars (for mean)
err = findobj(h(1).Parent, 'type', 'ErrorBar');
if ~isempty(err)
    xoverlaperr = (err.XData - xlim(1))./diff(xlim) <= (pos(1) + pos(3));
    errEnds = (err.YData + err.YPositiveDelta - ylim(1))./diff(ylim); % normalized
    errEnds = max(errEnds(xoverlaperr));
end

drawnow;
if errEnds <= pos(2) % no overalps
    return 
end

% increase YLim till it encloses the legend
ax = h(1).Parent;
step = diff(ylim)/100;
while errEnds > pos(2)
    drawnow;
    ax.YLim(2) = ax.YLim(2) + step;
    leg.Location = loc;
    pos = leg.Position;
    if ~isempty(err)
        errEnds = (err.YData + err.YPositiveDelta - ax.YLim(1))./diff(ax.YLim); % normalized
        errEnds = max(errEnds(xoverlaperr));
    elseif ~isempty(iqr) % err and iqr are mutually exclusive 
        errEnds = reshape([iqr.Extent], 4, []).';
        errEnds = errEnds(:, 2) + errEnds(:, 4); % y-point + height
        errEnds = (errEnds - ax.YLim(1))./diff(ax.YLim); % normalized
        errEnds = max(errEnds(xoverlapiqr));
    else
        errEnds = ([h.YEndPoints] - ax.YLim(1))./diff(ax.YLim); % normalized
        errEnds = max(errEnds(xoverlap)); % tallest bar within legend area
    end
end
end

%% ========================================================================
function h = adjNfont(h)
drawnow;
t = findobj(h(1).Parent.Children, 'Tag', 'n');
nIdx = ismember({t.String}, 'N =');
npos = t(nIdx).Position;
idx = false(numel(t), 1);
for i = 1:numel(t)
    if t(i).Position(2) == npos(2)
        idx(i) = true;
    end
end
[t(idx).Units] = deal('normalized');
idx(nIdx) = false;
idx = find(idx);

xpos = [t(idx).Position]; xpos = xpos(1:3:end);
[~, xpos] = sort(xpos);
idx = idx(xpos);

fontsz = t(idx(1)).FontSize - 1;

% check if text margins overlap
ext = reshape([t(idx).Extent], 4, []).';
while any(ext(1:end-1, 1) + ext(1:end-1, 3) >= ext(2:end, 1))
    drawnow;
    [t(idx).FontSize] = deal(fontsz);
    ext = reshape([t(idx).Extent], 4, []).';
    fontsz = fontsz - 1;
end

% use abbreviation for strata
t = findobj(h(1).Parent.Children, 'Tag', 'strata');
if ~isempty(t)
    nIdx = numel(t);
    if contains(t(nIdx).String, "(") && contains(t(nIdx).String, ")") && strlength(t(nIdx).String) > 25
        t(nIdx).String = extractBetween(t(nIdx).String, "(", ")");
    end

    if strlength(t(nIdx).String) > 8
        t(nIdx).FontSize = .6*t(nIdx).FontSize;
    end
end
end