function GeneTestPlotter(out, opts)
% plots summary statistics from gene-basted tests, where 'out' is a table
% from REGENIE gene-based rare variant analysis, and 'test' can be tests
% supported in REGENIE (and found in 'out' table), such as "ADD",
% "ADD-SKAT", ...
% Oveis Jamialahmadi, 06 Feb 2023. University of Gothenburg.
% 
% @08MAY2023: some bugs were fixed.

arguments
    out {mustBeA(out, 'table')}
    opts.test {mustBeVector, mustBeText} = "ADD" % test(s) to be kept
    opts.markerSize (1,1) double = 220 % scatter marker size
    opts.graph (1,1) {mustBeA(opts.graph, ["digraph", "graph"])}
    opts.resolution (1,1) double = 400
    opts.fontname {mustBeTextScalar} = "Garamond"
    opts.fontsize (1,1) double = 12
    opts.legendfontsize (1,1) double % legend font size. if left empty round(fontsize/1.5) will be used.
    opts.tickRotation (1,1) double = 65 % Xtick rotation
    opts.output {mustBeText}
    opts.savefig (1,1) logical = false
    opts.save (1,1) logical = true
    opts.arrow {mustBeMember(opts.arrow, ["normal", "bold"])} = "normal" % arrow line with
    opts.numcols (1,1) double = 3 % total columns for tiledlayout, where last column is kept for legend. Increase this to allow wider range for effect plots

    % for  genome-wide, 5e-8 will be used; for fdr, BH adjusted p-value
    % will be used. If set to none, no line will be plotted.
    opts.significance {mustBeTextScalar, mustBeMember(opts.significance, ["genome-wide", "bonferroni", "fdr", "none", "nominal"])} = "fdr" 
    opts.hideNonSignificant (1,1) logical = false % hide non-significant hits

    opts.pruneP (1,1) = 1 % remove rows with a P above this, useful when the plot is crowded. This option will be applied AFTER multiple testing, so it doesn't affect the adjusted P-values
end

if isfield(opts, 'output'), opts.output = string(opts.output); end
if isfield(opts, 'pruneP'), opts.pruneP = -log10(opts.pruneP); end

pheno = unique(out.Pheno);
for i = 1:numel(pheno)
    tmp = out(out.Pheno == pheno(i), :);
    if height(tmp) < 2, continue; end
    opts.pheno = pheno(i);
    fig = drawPerPheno(tmp, opts);
    if ~isfield(opts, 'output')
        outFileName = matlab.lang.makeValidName(pheno(i)) + ".genePlot";
    else
        pth = fileparts(opts.output);
        if pth ~= "" && ~isfolder(pth)
            mkdir(pth)
        end
        outFileName = opts.output + "." + matlab.lang.makeValidName(pheno(i));
    end
    
    if opts.save
        exportgraphics(fig.Parent, outFileName + ".jpeg", "Resolution", opts.resolution)
    end

    if opts.savefig
        savefig(fig.Parent, outFileName + ".fig")
    end
end

end % END

%% subfunctions ===========================================================
function fig = drawPerPheno(out, opts)

idx = out.ID.contains("(");
out.Gene = out.ID;
out.Gene(idx) = extractBefore(out.ID(idx), "(");
out.Gene(ismissing(out.Gene)) = out.ID(ismissing(out.Gene));
out.ID = out.Gene;
% out.Gene = out.ID;

convcols = ["Gene", "Mask", "GO"];
hasGO = true;
numcols = opts.numcols;
if ~any(colnames(out) == "GO")
    convcols(end) = [];
    hasGO = false;
end
out = convertvars(out, convcols, @categorical);

% convert P to -log10P
tests = strings(numel(opts.test), 1);
for k = 1:numel(opts.test)
    if contains(opts.test(k), "-") % remove "ADD-" prefix 
        test = "-log_{10} P_{" + opts.test(k).extractAfter("-") + "}";
    else
        test = "-log_{10} P_{" + opts.test(k) + "}";
    end
    
    out.(test) = -log10(out.(opts.test(k)));
    tests(k) = test;
end

close gcf force
% fig = figure("WindowState", "maximized");
% ax = axes(fig);
fig = tiledlayout(numel(opts.test), numcols, TileSpacing="tight", Padding="compact");
fig.Parent.WindowState = "maximized";

for i = 1:numel(opts.test)
    ax1 = nexttile(fig, tilenum(fig, i, 1), [1, numcols-1]);
    tmp = out;
    
    test = tests(i);
    
    idx = find(startsWith(lower(tmp.Properties.VariableNames), ...
        ["beta|or", "or|beta", "β|or", "or|β"]));
    if isempty(idx)
        idx = find(startsWith(lower(tmp.Properties.VariableNames), ["or", "beta", "β"]));
    end
    
    if isempty(idx) % beta column is missing?
        tmp.Beta = zeros(height(tmp), 1);
    else
        tmp.Beta = tmp.(idx); % beta or log OR
    end

    % OR --> log OR
    if any(colnames(tmp) == "N.case") && all(~isnan(tmp.("N.case")))
        tmp.Beta = log(tmp.Beta);
    end
    tmp.Beta(ismissing(tmp.Beta) | isinf(tmp.Beta)) = 0;

    % Bubblechart
%     tmp.Beta = categorical(sign(tmp.Beta));
%     tmp.Beta(tmp.Beta == "-1") = "<0";
%     tmp.Beta(tmp.Beta == "1") = ">0";
    
%     tmp.sz = ones(height(tmp), 1);
%     b = bubblechart(tmp, "Gene", "-log_{10} P_{ADD}", "sz", "Beta");
%     ylabel(b.YVariable, "Interpreter", "tex")
%     bubblesize([8 14]);
%     blgd = bubblelegend('Beta', 'NumBubbles', 2);
%     blgd.LimitLabels = {'<0', '>0'};

    % scatter
   
%     % get RGB for beta
%     cmap = parula;
%     cmin = min(tmp.Beta);
%     cmax = max(tmp.Beta);
%     m = length(cmap);
%     index = fix((tmp.Beta-cmin)/(cmax-cmin)*m)+1;
%     Cdata = squeeze(ind2rgb(index,cmap));
    
    tmp.Gene = categorical(string(tmp.Gene));

    if hasGO
        % color GO
        gou = string(unique(tmp.GO));
        cmap = parula(numel(gou));
    
        % group by GO: some gene/mask combinations may map to the same GO. So,
        % in order to have color coded GO, we keep all GOs per each gene/mask
        % combination an add colors on top with different marker size for
        % visualization
        tmpg = groupsummary(tmp, ["ID", "Mask", test], @(x) join(string(x),";"), "GO");
        [f1, f2] = ismember(tmp.ID+string(tmp.Mask), tmpg.ID+string(tmpg.Mask));
        tmp.GO = strings(height(tmp), 1);
        tmp.GO(f1) = tmpg.fun1_GO(f2(f1));
    else
        cmap = [1, 0, 0];
    end

    tmp.Mask = string(tmp.Mask);
    mask = unique(tmp.Mask);
    mrkr = ["o", "square", "diamond", "^", "v", ">", "<", "pentagram", "hexagram", "x", "*", "+", "."];
    mrkr = mrkr(1:numel(mask));
    hold on

    % multiple testing adjustment for each test
    opts.sig_all = applyMultiTest(tmp.(opts.test(i)), opts);
    
    pruned_idx_outer = any(tmp{:, tests} >= opts.pruneP, 2); % to keep union of genes for all tests
    tmp_clean = tmp(pruned_idx_outer, :);
    tmp_clean.Gene = removecats(tmp_clean.Gene);
    
    % plot gene vs -log10 P-val with GO color codes and different markers
    % for masks
    s = cell(numel(mask), 1);

    sig_idx_join = false(height(tmp), 1); % for non-significant hits after multiple testing adjustment
    for j = 1:numel(mask)
        mask_idx = tmp.Mask == mask(j);
        tmp2 = tmp(mask_idx, :);
    
        % multiple testing adjustment for each mask (suggestive threshold)
        [opts.sig, sig_idx] = applyMultiTest(tmp2.(opts.test(i)), opts);
        pruned_idx = any(tmp2{:, tests} >= opts.pruneP, 2);
        tmp2_clean = tmp2(pruned_idx, :);
        sig_idx(~pruned_idx) = [];
        mask_idx_tmp = mask_idx;
        mask_idx_tmp = find(mask_idx_tmp);
        mask_idx_tmp(~pruned_idx) = [];
        mask_idx(:) = false;
        mask_idx(mask_idx_tmp) = true;
        tmp2_clean.Gene = setcats(tmp2_clean.Gene, string(tmp_clean.Gene));

        if hasGO
            [f1, f2] = ismember(string(tmp2.GO), gou);
            h = scatter(ax1, tmp2_clean(f1, :),'Gene', test, 'filled', ...
               'Marker', mrkr(j), 'CData', cmap(f2(f1), :), ...
               'MarkerEdgeColor', 'k', 'SizeData', opts.markerSize, ...
               'MarkerFaceAlpha', 0.8, "MarkerEdgeAlpha", 0.4);
    
            if any(~f1) % mult GO per each mask/ID pair
                tmp3 = tmp2_clean(~f1, :);
                [~, uidx] = unique(tmp3.ID+string(tmp3.Mask)+string(tmp3.GO));
                tmp3 = tmp3(uidx, :);
                for k = 1:height(tmp3)
                    gGO = split(tmp3.GO(k), ";");
                    [f1, f2] = ismember(gGO, gou);
                    cmapGO = cmap(f2(f1), :);
                    for k2 = 1:size(cmapGO, 1)
                        if k2 > 1, ecol = "none"; mscale = opts.markerSize*k2/70; else, ecol = "k"; mscale = 1; end
                        scatter(ax1, tmp3(k, :),'Gene', test, 'filled', ...
                           'Marker', mrkr(j), 'CData', cmapGO(k2, :), ...
                           'MarkerEdgeColor', ecol, 'SizeData', opts.markerSize/mscale, ...
                           'MarkerFaceAlpha', 0.9, "MarkerEdgeAlpha", 0.4);
                    end
                end
    
            end

        else % no GO
            h = scatter(ax1, tmp2_clean,'Gene', test, 'filled', ...
               'Marker', mrkr(j), 'CData', cmap, ...
               'MarkerEdgeColor', 'k', 'SizeData', opts.markerSize, ...
               'MarkerFaceAlpha', 0.8, "MarkerEdgeAlpha", 0.4);
        end
        
        % modify hits based on multiple testing method
        if any(~sig_idx)
            h.SizeData = ones(height(tmp2_clean), 1).*opts.markerSize;
            h.AlphaData = ones(height(tmp2_clean), 1).*h.MarkerFaceAlpha; 
            h.MarkerFaceAlpha = 'flat';
            if opts.hideNonSignificant
                h.MarkerEdgeAlpha = 'flat'; % to hide non-significant hits
            end

            if ~any(sig_idx) % nothing's significant
                h.MarkerFaceAlpha = 0;
            else
                h.AlphaData(~sig_idx) = 0.4;
            end
        end
        sig_idx_join(mask_idx) = sig_idx;

        ylabel(ax1, test, "Interpreter", "tex")
        s{j, 1} = h;

    end % end of enumerating masks
    
    % Significance threshold for all masks jointly
    if ~isnan(opts.sig_all)
        yline(ax1, [opts.sig_all, opts.sig_all],...
            'LineWidth', 1, ...
            'color', 'r', 'LineStyle', '--') 
    end
    
    % if numel(opts.test) == 1, axis(ax1, "square"); end
    if i == 1, title(ax1, opts.pheno, "FontSize", opts.fontsize + 2); end
    ax1.FontName = opts.fontname;
    ax1.XAxis.FontSmoothing = "on";
    ax1.XAxis.FontSize = opts.fontsize;
    ax1.YAxis.FontSmoothing = "on";
    ax1.YAxis.FontSize = opts.fontsize;
    ax1.XLabel.String = "";
    if i < numel(opts.test) % && (~isfield(opts, "pruneP") || opts.pruneP <= 0)
        ax1.XTickLabel = {};
    else
        ax1.XTickLabelRotation = opts.tickRotation;
        if any(sig_idx_join) % bold significant hits after adjusting for multiple testing
            sigIDs = tmp.ID(sig_idx_join);
            ax1.XTickLabel = replace(ax1.XTickLabel, sigIDs, "\bf " + sigIDs);
            % [~, min_p_idx] = max(tmp.(test));
            min_p_idx = tmp.(test) >= opts.sig_all;
            sigIDTop = tmp.ID(min_p_idx);
            ax1.XTickLabel = replace(ax1.XTickLabel, sigIDTop, "\bf \color{red} " + sigIDTop);
        end
    end

    % show beta directions
    idxup = tmp_clean.Beta > 0;
    idxdown = tmp_clean.Beta < 0;
    if opts.hideNonSignificant
        idxup = idxup & sig_idx_join(pruned_idx_outer);
        idxdown = idxdown & sig_idx_join(pruned_idx_outer);
    end
    text(ax1, tmp_clean.Gene(idxup), tmp_clean.(test)(idxup), ...
        "\uparrow", "HorizontalAlignment", "center", ...
        "VerticalAlignment", "bottom", "FontWeight", opts.arrow);

    text(ax1, tmp_clean.Gene(idxdown), tmp_clean.(test)(idxdown), ...
        "\downarrow", "HorizontalAlignment", "center", ...
        "VerticalAlignment", "cap", "FontWeight", opts.arrow);   

    shade(h.Parent, 'dir', 'x');
end

% add legends and GO graph 
span = floor(numel(opts.test)/2);
ax1 = nexttile(fig, tilenum(fig, 1, numcols), [max(1, span), 1]);
ax1.Visible = 'off';
hold(ax1, "on")
lg = (0);
for j = 1:numel(mask)
    % for legend
    lg(j) = plot(ax1, nan, nan, 'MarkerFaceColor', [0.5 0.5 0.5],...
        'MarkerEdgeColor', 'k', 'LineStyle', 'none',...
        'Marker', mrkr(j));
end
hold(ax1, "off")

% mask legend
mtch1 = ["^mask_", "lof_", "lof", "(\W{1}|_|&|\|AND|OR)(c|CAD|CADD)(?![a-zA-Z])", "(&|\|)sr", "(&|\|)pr", "(\W{1}|_|&|\|AND|OR)(r|REVEL)(?=(0|5))", "(&|\|)cvar"];
mtch2 = ["", "LoF & ", "LoF", " $1 CADD≥", " $1 SIFT", " $1 PolyPhen"," $1 REVEL≥", " $1 ClinVar"];
mskstr = regexprep(mask, mtch1, mtch2);
mskstr = regexprep(mskstr, "\s+", " ");
mafstr = extractAfter(mskstr, ".");
mafstr_cats = unique(rmmissing(mafstr));
idx = ismissing(mafstr); % Joint test?
mskstr(~idx) = erase(mskstr(~idx), "." + mafstr(~idx) + textBoundary("end"));
if ~isscalar(mafstr_cats) % keep MAFs in the mask
    mskstr(~idx) = mskstr(~idx) + " (MAF:" + mafstr(~idx) + ")";
end

mskstr = regexprep(mskstr, "(REVEL(≥|>))(\d+)", "$10.$2"); % convert REVEL cutoff to floating point

if ~isfield(opts, 'legendfontsize')
    opts.legendfontsize = round(opts.fontsize/1.5);
end
leg = legend(lg, 'AutoUpdate', 'off',...
    'String', mskstr, 'Location', 'northwest',...
    'Color', '#BAD8E0', "FontSize", opts.legendfontsize, ...
    "FontName", opts.fontname);

if isscalar(mafstr_cats)
    leg.Title.String = "Mask (MAF " + mafstr_cats + ")";
else
    leg.Title.String = "Mask";
end
leg.TextColor = '#3F3F3F';

% GO legend
if hasGO
    % ax2 = copyobj(ax1, gcf);
    ax2 = axes('position', ax1.Position,'visible','off');
    drawnow;
    % ax2.Visible = 'off';
    for j = 1:numel(ax2.Children)
        ax2.Children(j).Visible = "off";
        try ax2.Children(j).Position = ax1.Children(j).Position; catch, end
    end
    
    lg2 = (0);
    hold(ax2, "on")
    for j = 1:size(cmap, 1)
        lg2(j) = plot(ax2, NaN, NaN, 'MarkerFaceColor', cmap(j,:),...
                'LineStyle', 'none', 'Marker', 'O',...
                'MarkerEdgeColor',  'k',...
                'Color',  cmap(j,:), 'MarkerSize', 5);
    end
    hold(ax2, "off")
    leg2 = legend(lg2, gou, 'Color', [0.93 0.93 0.93], 'AutoUpdate', 'off',...
        'Location', 'southwest', 'Orientation', 'Vertical', ...
        'Color', '#BAD8E0', "FontSize", opts.legendfontsize, ...
        "FontName", opts.fontname);
    leg2.Title.String = "GO:BP";
    leg2.TextColor = '#3F3F3F';
    ax2.Position = ax1.Position;
    drawnow;
    leg2.Position(1) = leg.Position(1);
    leg2.Position(2) = leg.Position(2) - leg2.Position(4) - 0.005;
end

% add graph to the axis
if ~isfield(opts, "graph"), return; end

ax3 = nexttile(fig, tilenum(fig, span+1, 3), [numel(opts.test)-span, 1]);
h = plot(ax3, opts.graph, "Layout","layered", "EdgeColor", "#4682B4", ...
    "NodeColor", "#DAA520", "LineWidth", 1, "MarkerSize", 5);
l = string(h.NodeLabel);
h.NodeLabel = {};
ax3.Visible = "off";

% tiledlayout doesn't allow to modify these properties
% drawnow;
% try    
%     ax3.Position(1) = leg2.Position(1);
% catch
%     ax3.Position(1) = leg.Position(1);
% end
% ax3.OuterPosition(2) = 0;
% ax3.Position(3) = 0.9 - ax3.Position(1);
% try
%     ax3.OuterPosition(4) = leg2.Position(2) - ax3.OuterPosition(2); % not sure?
% catch
%     ax3.OuterPosition(4) = leg.Position(2) - ax3.OuterPosition(2);
% end
% l_go = extract(l, "GO:"+digitsPattern+":");
l = replace(l, ["negative", "positive", "regulation"], ["neg", "pos", "reg"]);
l = regexp(l', "(?<=GO:\d+):", "split");
l = vertcat(l{:});
text(ax3, h.XData, h.YData, l(:,1) + newline + l(:, 2), ...
    "FontSize", 5, "VerticalAlignment", "cap", "Rotation", -15);

end % END

%% subfunctions ===========================================================
function [sig, sig_idx] = applyMultiTest(pvec, opts)

opts.significance = string(opts.significance);
if opts.significance == "bonferroni"
    sig = 0.05/height(tmp2);
elseif opts.significance == "genome-wide"
    sig = 5e-8;
elseif opts.significance == "fdr"
    sig = mafdr(pvec, 'BHFDR', true);
    sig = max(pvec(sig <= 0.05)); % max p-value passing BH FDR cutoff
    if isempty(sig); sig = 0; end
elseif opts.significance == "nominal"
    sig = 0.05;
else
    sig = nan;
end

sig_idx = pvec <= sig;
if isnan(sig)
    sig_idx = true(numel(pvec), 1);
end
sig = -log10(sig);

end