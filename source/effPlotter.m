function effPlotter(tab, opts)
% draws forest plot(s) for a set of input traits and variants (on y-axis).
% effPlotter can infer type of the input trait (binary/continuous) from
% the summary statistics, and offers multiple plots (e.g. displaying
% p-values in addition to effect size or merging all into one single plot).
% Oveis Jamialahmadi, University of Gothenburg, March 2022.
% 
% @25NOV2022: tab can ba single table for all traits. Note that in this
%             case 'binary' must have the same size as height ot table.
% @26NOV2022: 'groupCol' option was added to identify grouping variable in
%              the input table. The grouping varialbe will be shown in
%              legend (default it "Pheno").
%              'pheno' is also deprecated and 'groupLabels' should be used
%              instead for labeling grouping variable. 
%              'yCol' option was added to identify the variable to be shown
%              on y-axis (default is "SNP"). 'yLabels' was added which is
%              similar to 'groupLabels' but for 'yCol' variable. 'yLabels'
%              is displayed on y-axis.
% @15DEC2022: a bug was fixed.
% @16DEC2022: 'ignoreCI' flag was added (default: false). If true, does not
%             use, estimate and plot 95% CI. Some bugs were fixed with
%             'binary' option. Now should be in the same size/structure as
%             'tab'.
% @2JUNE2023: 'betaColumn' was added for cases where effect sizes are
%             stored in a non-canonical variable name. 
% @18JULY2023: 'groupOrder' was added to change order of 'groupCol' values.
%             This is useful when ther order of strings in 'groupCol'
%             matters for the visualization.
%             'marker' option was added for customized markers.
%             'markersize' can be a vector of different sizes.
% @26JULY2023: 'squeeze' flag was added to squeeze the plot whenever the
%             beta or CIs are too wide. Use isoutlier built-in function
%             with and removes those 2 SD away from the mean. For CIs,
%             arrows show if a certain CI extends further from the current
%             axes.
% @29AUG2023: 'yOrder' option was added. Note that 'sortEffect' overrides
%             this option. 'yOrder' acts similarly to 'groupOrder'.
% @11MAR2024: 'hidealpha' option was added to control the transparency when
%             'hide' flag is true.
% @08OCT2024: 'xline' option was added. 'xline' decides the position of
%             vertical x-line for effect size. 'hideCutoffLower'
%             (default=0.05) was added to hide values > this cut-off.
%             'hideCutoffUpper' was added (default=0) to hide values <=
%             this cut-off.

arguments
    tab {mustBeAllTab(tab)} % table for each phenotype (snp column must be unique)
    opts.groupCol = "Pheno" % grouping variable. "Pheno" shows variants on y-axis grouping all pheno for each variant together.
    opts.groupLabels {mustBeText, mustBeVector} % group variable label. Must have the same size as tab
    opts.yCol = "SNP" % the variable to be shown on y-axis. Usually this is the main term in the statistical model for which effect sizes are provided.
    opts.yLabels  {mustBeText, mustBeVector} % similar to groupLabels but for 'yCol' variable.
    opts.pheno {mustBeText, mustBeVector} % phenotype label must have the same size as tab (DEPRECATED. 'groupLabels' should be used now)
    opts.binary logical {mustBeVector} % a logical vector/cell (same structure of 'tab'). If tab is table, it should be a vector of the same height, if cell, should be a cell array of logical vectors. If absent, P column must be present.
    opts.markersize {mustBeNumeric, mustBeVector} % marker size (either a vector corresponding to lables in 'groupLabels' or a single value). default is 55 for all methods, except "scatter" which is 5
    opts.markeralpha (1,1) double = 1
    opts.showp (1,1) logical = false 
    opts.showpmethod {mustBeTextScalar, mustBeMember(opts.showpmethod, ["bubble", "color", "scatter"])} = "color"
    opts.bubbleSize (1, 2) double = [5 10] % bubble size if showpmethod is "bubble"
    opts.merge (1,1) logical = true % to merge all traits into a single plot. overrides showp flag
    opts.overlap (1,1) logical = false % to allow (true) overlapping effect sizes or non-overlapping markers (false) when using 'merge' flag with multiple traits. 
    opts.fontsize (1,1) double = 11
    opts.fontname {mustBeTextScalar, mustBeMember(opts.fontname, ["Bauhaus", "Futura", "Bodoni", "Garamond", "Corbel", "Corbel Light", "Arial", "Tahoma"])} = "Arial"
    opts.xbold (1,1) logical = false % x-axis font weight
    opts.ybold (1,1) logical = false % y-axis font weight
    opts.output {mustBeTextScalar} = ""
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "png" % format of output plots
    opts.resolution (1, 1) double = 300 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    opts.sortEffect (1,1) logical = true % sort effect sizes
    opts.betaColumn {mustBeTextScalar} % effect sizes are stored in this variable name. By default look for beta column.
    opts.xlabel {mustBeTextScalar} = "" % xlabel. By default it's beta or OR
    
    % developer options
    % changes fontsize or position of colorbar's or legend's title so that
    % it doesn't cross the upper x-axis (only if 'showp' is true)
    opts.adjmethod {mustBeMember(opts.adjmethod, ["fontsize", "pos"])} = "pos"
    opts.color (1,:) {mustBeMember(opts.color, ["Coral", "OrangeRed",...
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
        "BlueLight", "BrutalBlue", "CrownYellow", "YassQueen", "SisterSister"])} = "DeepSkyBlue" % can be > 1 if 'merge' is true. If 'merge' is true and numel is 1, 'colormap' will be used.
    opts.WindowState {mustBeTextScalar, mustBeMember(opts.WindowState, ["maximized", "normal"])} = "maximized" 
    opts.WindowPosition (1,4) double = [1, 41, 1280, 720] % windows position in pixels, overrides WindowState
    opts.titleBuffer (1,1) double = 4 % buffer between title and axis (DEPRECATED). Now 'TitleHorizontalAlignment' is used.
    opts.TileSpacing {mustBeTextScalar, mustBeMember(opts.TileSpacing, ["tight", "loose", "compact"])} = "tight"
    opts.visible {mustBeTextScalar, mustBeMember(opts.visible, ["on", "off"])} = "on" % use off only for debugging/rendering issues
    opts.square (1,1) logical = true % make axis square
    opts.verbose (1,1) logical = true
    opts.shade (1,1) logical = false 
    opts.cbarSqueeze (1,1) double = 2 % factor to squeeze p-value colorbar horizontally.
    opts.addTitle (1,1) logical = true % add title per each plot (merge should be false)
    opts.Title {mustBeTextScalar} % optional title instead of group labels
    opts.titleFontSize (1,1) double = 12
    opts.grid (1,1) logical = true
    opts.box (1,1) logical = true
    opts.ignoreCI (1,1) logical = false % If true, does not use, estimate and plot 95% CI.
    opts.legendOrientation {mustBeTextScalar, mustBeMember(opts.legendOrientation, ["horizontal", "vertical"])} = "vertical" % legend orientation
    opts.legendCols (1,1) double = 1 % legend columns
    opts.legLocation {mustBeMember(opts.legLocation, ["north", "south", ...
        "east", "west", "northeast", "northwest", "southeast", "southwest", ...
        "northoutside", "southoutside", "eastoutside", "westoutside", ...
        "northeastoutside", "northwestoutside", "southeastoutside", ...
        "southwestoutside", "best", "bestoutside", "none", "layout"])} = "northeastoutside"
    opts.legTitle {mustBeTextScalar} % legend title
    opts.legFontOffset (1,1) double = 2 % legend font size offset (e.g. if fontsize is 16, and offset is 2, legend fontsize is 14)
    opts.ciLineWidth (1,1) double = 1.1 % CI line width
    
    opts.hide (1,1) logical = false % to empty the color of non-significant markers
    opts.hidealpha (1,1) double = 0.3
    opts.colormap {mustBeMember(opts.colormap, ["jet", "turbo", "abyss", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", "summer", ...
        "colorcube", "lines", "magma", "plasma", "inferno", "cividis", ...
        "viridis", "twilight", "twilight_shifted", "turbo", "Blues", ...
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
        "flare", "crest", "colorblind", "deep", "muted"])} = "turbo" % colormap for when 'color' is not used.
    opts.groupOrder {mustBeText, mustBeVector} % order of group labels in 'groupCol' column. This option is still under development and for now works only when tab is a table and "groupLabels" is not used.
    opts.yOrder {mustBeText, mustBeVector}
    opts.marker % a table or string vector of marker types (default is 'O'). In case of a table, first column should be group lables and second column is markers.
    opts.squeeze (1,1) logical = false % for better visualization removes outliers
    opts.zeroEffect (1,1) logical = true % to show a vertical line (1 for OR/HR and 0 for beta) denoting the null effect.

    %@17JULY2024
    opts.ycol2 {mustBeTextScalar} % the variable to be shown on the secondary y-axis to the right (only works with 'merge' as true)
    opts.fontsizeYcol2 (1,1) double = 8 % fontsize of ycol2

    %@22JULY2024
    opts.df {mustBeTextScalar} % degree of freedom to calculate 95%CI (stutent t-dist)

    %@08OCT2024
    opts.xline (1,1) double
    opts.hideCutoffLower (1,1) double = 0.05
    opts.hideCutoffUpper (1,1) double = 0
    opts.ycol2Label {mustBeTextScalar} = "" % ycol2 label

    %@27MARCH2025
    opts.xscale {mustBeMember(opts.xscale, ["linear", "log"])} = "linear"
    opts.squeezeOffset (1,1) double % squeeze offset around effect -/+ CI
end

% check 'groupCol' and 'yCol' options
if all(isfield(opts, ["groupCol", "yCol"]))
    if strcmp(opts.groupCol, opts.yCol)
        if opts.verbose
            disp('groupCol and yCol cannot point to the same column!')
        end
        opts.yCol = setdiff(["SNP", "Pheno"], opts.groupCol);
        opts.yCol = opts.yCol(1);
    end
end

% check and munge if input is a table
if istable(tab)
    if any(ismissing(tab.(opts.groupCol)))
        tab = fillmissing(tab, "previous", "DataVariables", opts.groupCol);
    end

    if isfield(opts, 'yCol') && any(ismissing(tab.(opts.yCol)))
        tab = fillmissing(tab, "previous", "DataVariables", opts.yCol);
    end

    if isfield(opts, 'binary') && numel(opts.binary) ~= height(tab)
        error('effPlotter:binary', 'binary and tab must have the same size!')
    end
    
    % group stats based on 'groupCol'
    ugrp = string(unique(tab.(opts.groupCol), "stable"));
    
    % check ordering of groups 
    if isfield(opts, "groupOrder")
        [~, idx] = ismember(flip(opts.groupOrder), ugrp);
        ugrp = ugrp(idx(idx > 0));
    end

    ttab = cell(numel(ugrp), 1);
    if isfield(opts, 'binary'), bbinary = ttab; end
    for i = 1:numel(ugrp)
        idx = ismember(tab.(opts.groupCol), ugrp(i));
        ttab{i, 1} = tab(idx, :);

        if isfield(opts, 'binary')
            bbinary{i, 1} = opts.binary(idx);
        end
    end
    
    if isfield(opts, 'binary'), opts.binary = bbinary; end
    tab = ttab;
end

% opts.pheno is deprecated
if isfield(opts, "pheno"), opts.groupLabels = opts.pheno; end

% check groupLabels
if isfield(opts, 'groupLabels') && numel(opts.groupLabels) ~= numel(tab)
    error('effPlotter:groupLabels', 'groupLabels and tab must have the same size!')
elseif ~isfield(opts, 'groupLabels')
    opts.groupLabels(1:numel(tab)) = missing; 
else
    opts.groupLabels = string(opts.groupLabels);
    missIdx = opts.groupLabels == "" | ismissing(opts.groupLabels);
    opts.groupLabels(missIdx) = missing;
end
opts.groupLabels = string(opts.groupLabels);

% check yLables (compares to height of first tab)
if isfield(opts, 'yLabels') && numel(opts.yLabels) ~= height(tab{1})
    error('effPlotter:yLabels', 'yLabels should be of the same size as the underlying tables!')
elseif ~isfield(opts, 'yLabels')
    opts.yLabels(1:height(tab{1})) = missing;
else
    opts.yLabels = string(opts.yLabels);
    missIdx = opts.yLabels == "" | ismissing(opts.yLabels);
    opts.yLabels(missIdx) = missing;
end
opts.yLabels = string(opts.yLabels);

if isfield(opts, 'binary') && numel(opts.binary) ~= numel(tab)
    error('effPlotter:binary', 'binary and tab must have the same size!')
elseif ~isfield(opts, 'binary')
    for i = 1:numel(tab)
        opts.binary{i, 1} = nan(height(tab{i}), 1);
    end
end

opts.output = string(opts.output);
if ismissing(opts.output) || opts.output == ""
    opts.output = getRandomName("effPlot", 3);
end

plt = palette;
opts.color = string(opts.color);
sz = numel(tab);
if opts.merge
    if numel(opts.color) > sz, opts.color = opts.color(1:sz); end
    if numel(opts.color) ~= sz % pick random colors
%         idx = randsample(height(plt), sz); % throws error if numel(tab) > height(plt) --> quite rare
%         opts.color = plt.color(idx);
        try
            opts.color = feval(opts.colormap, sz);
        catch
            opts.color = get_colormap(opts.colormap, sz);
        end
    else
        [~, idx] = ismember(opts.color, plt.color);
        opts.color = plt.code(idx);
    end
else
    if numel(opts.color) > 1; opts.color = opts.color(1); end
    opts.color = plt.code(ismember(plt.color, opts.color));
end

% unify table headers and check if CI is present or not
for i = 1:numel(tab)
    cols = lower(tab{i}.Properties.VariableNames);
    oldcol.ci = find(contains(cols, "ci"), 1);

    if opts.ignoreCI || isempty(oldcol.ci) % if 95% is absent and should be estimated (P must be present in this case)
        opts.ci(i) = false;
    else
        opts.ci(i) = true;
%         newcols(5) = "ci";
    end
    
    if isfield(opts, 'yCol')
        checkSNP = lower(opts.yCol);
    else
        checkSNP = {'snp', 'annotate', 'variant', 'variantid', 'id'};
    end
    oldcol.snp = find(ismember(cols, checkSNP), 1);
    if isfield(opts, "betaColumn")
        oldcol.b = find(cols == string(lower(opts.betaColumn)));
    else
        oldcol.b = find(ismember(cols, {'b', 'β'}) | contains(cols, ["beta", "or"]), 1);
    end
    oldcol.se = find(ismember(cols, {'se'}), 1);

    if ismissing(opts.groupLabels(i)) % check if there is a pheno column for title
        if isfield(opts, 'groupCol')
            checkPheno = lower(opts.groupCol);
        else
            checkPheno = ["pheno", "trait", "tag", "label", "title"];
        end
        phenocol = find(contains(cols, checkPheno), 1);
        if isempty(phenocol)
            opts.groupLabels(i) = "";
        else
            phenocol = string(tab{i}.(phenocol));
            phenocol(ismissing(phenocol) | phenocol == "") = [];
            if isempty(phenocol)
                phenocol = "";
            else
                phenocol = phenocol(1); % sould be unique
            end
            opts.groupLabels(i) = phenocol;
        end
    end
    
    if opts.ignoreCI && isempty(oldcol.b)
        error('effPlotter:emptyCols', 'Beta column is needed!')
    elseif ~opts.ignoreCI && numel([oldcol.b, oldcol.se]) < 2
        error('effPlotter:emptyCols', 'Beta and SE columns are needed!')
    end

    oldcol.p = find(ismember(cols, {'padj', 'p-adj', 'p_adj'}), 1);
    if isempty(oldcol.p)
        oldcol.p = find(ismember(cols, {'p', 'p_value', 'pval', 'pvalue'}), 1);
    end

    if ~opts.ignoreCI
        if isempty(oldcol.p)
            if ~opts.ci(i)
                error('effPlotter:emptyP', 'P column must be present when CI doesn''t exist!')
            elseif any(isnan(vertcat(opts.binary{:})))
                error('effPlotter:emptyP', 'P column must be present when binary flag is not provided!')
                
            end
        else
            opts.p(i) = true;
%             newcols(4) = "p";
        end
    else
        opts.p(i) = true;
    end

    emptfis = structfun(@isempty, oldcol);
    fis = string(fieldnames(oldcol));
    oldcol = rmfield(oldcol, fis(emptfis));
    
    %@17JULY2024: allow for ycol2 and df
    selCols = struct2array(oldcol);
    if isfield(opts, "ycol2")
        keepCols = [selCols; find(colnames(tab{i}) == opts.ycol2)];
    else
        keepCols = selCols;
    end

    if isfield(opts, "df")
        keepCols = union(keepCols, find(colnames(tab{i}) == opts.df), "stable");
    end

    tab{i} = tab{i}(:, keepCols);

    tab{i} = renamevars(tab{i}, colnames(tab{i}, index=1:numel(selCols)), fieldnames(oldcol));
    if isfield(opts, "df")
        tab{i} = renamevars(tab{i}, opts.df, "df");
    end
    % check if binary trait (in this case OR must be provided instead of
    % Beta)
    if ~opts.ignoreCI
        if any(isnan(opts.binary{i}))
            tab{i}.p(tab{i}.p <= realmin) = realmin;
            check_p = -log10(chi2pval((tab{i}.b./tab{i}.se).^2, 1));
            check_p(isinf(check_p)) = -log10(realmin);
            check_or = max([check_p ,-log10(tab{i}.p)], [], 2)./...
                    min([check_p ,-log10(tab{i}.p)], [], 2) < 1.05;
            opts.binary{i} = ~check_or;
        end
    
        % add CI
        tab{i} = getCI(tab{i}, opts.ci(i), opts.binary{i});

        % if a mixture of binary/quantitative traits are present, show beta
        % instead (yCol is phenotype)
        if any(opts.binary{i}) && any(~opts.binary{i})
            tab{i}(opts.binary{i}, :) = convertvars(tab{i}(opts.binary{i}, :), ["b", "ci_l", "ci_h"], @log);
            opts.binary{i} = false(height(tab{i}), 1);
        end
    end
    
    if i == 1
        snps = tab{i}.snp;
    else
        snps = intersect(snps, tab{i}.snp, "stable");
    end
end

if all(~ismissing(opts.yLabels))
    if numel(unique(cellfun(@height, tab))) > 1
        if opts.verbose
            disp('WARNING: underlying tables are of different size, cannot set yLabels!')
        end
    else
        snps = opts.yLabels;
        if isrow(snps), snps = snps'; end
        for i = 1:numel(tab), tab{i}.snp = snps; end
    end
end

% check order of yLabels
if isfield(opts, "yOrder")
    opts.yOrder = flip(string(opts.yOrder));
    for k = 1:numel(tab)
        [~, idx] = ismember(opts.yOrder, tab{k}.snp);
        tab{k} = tab{k}(idx(idx>0), :);

        if isfield(opts, "binary")
            opts.binary{k} = opts.binary{k}(idx(idx>0), :);
        end

        if k == 1
            snps = tab{k}.snp;
        else
            snps = intersect(snps, tab{k}.snp, "stable");
        end
    end

end

% if a mixture of binary/quantiative traits are present, show beta instead
% (yCol is snp)
for i = 1:numel(tab)
    if ~opts.ignoreCI && any(opts.binary{i}) && any(~opts.binary{i})
        if all(opts.binary{i})
            opts.binary{i} = false(height(tab{i}), 1);
            tab{i} = convertvars(tab{i}, ["b", "ci_l", "ci_h"], @log);
        end
    end
end

if ~any(opts.p) || opts.merge; opts.showp = false; end

if isfield(opts, 'markersize')
    
    if numel(opts.markersize) < numel(opts.groupLabels)
        opts.markersize = ones(numel(opts.groupLabels), 1).*opts.markersize(1);
    else
        opts.markersize = opts.markersize(1:numel(opts.groupLabels));
    end
    opts.markersize = flip(opts.markersize);

else
    opts.markersize = ones(numel(opts.groupLabels), 1).*55;
    if opts.showp && opts.showpmethod == "scatter"
        opts.markersize = ones(numel(opts.groupLabels), 1).*3;
    end
end

%% draw now ///////////////////////////////////////////////////////////////
% warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
close all force
if all(~isempty(opts.WindowPosition) & ~ismissing(opts.WindowPosition))
    fig = figure(Position=opts.WindowPosition, Units="pixels", Visible=opts.visible);
elseif isfield(opts, "WindowState")
    fig = figure("WindowState", opts.WindowState, Visible=opts.visible);
else
    fig = figure(Visible=opts.visible);
end

if ~(opts.showp && opts.showpmethod == "bubble") && ~opts.merge
    ti = tiledlayout(fig, 1, numel(tab), 'TileSpacing', opts.TileSpacing);
end

if opts.merge
    % % under development: for meta-analysis
    % ti = tiledlayout(fig, 1, 4, 'TileSpacing', opts.TileSpacing);
    % h = nexttile(ti, [1, 3]);

    h = gca(fig);
    sax = (0);
    if isscalar(tab); opts.overlap = true; end
else
    opts.overlap = true;
end

if ~opts.overlap % generate y vector (different for each trait)
    yvec = linspace(0.65, 1.35, numel(tab));
end

% @26JULY2023: for better visualization, axes won't cover too wide CIs.
if opts.squeeze
    biflag = vertcat(opts.binary{:});
    bvec_raw = cellfun(@(x)x.b, tab, uni=false);
    bvec_raw = vertcat(bvec_raw{:});
    % % nanIdx = isnan(bvec_raw);
    % bvec = bvec_raw;
    % bvec(biflag) = log(bvec_raw(biflag)); % convert back to log OR
    % boutlier = isoutlier(bvec, "mean", ThresholdFactor=2);
    % bvec_raw(boutlier) = nan;
    % boutlier = reshape(boutlier, [], numel(tab));
    % for k = 1:numel(tab)
    %     tab{k}.b(boutlier(:, k)) = nan;
    %     if ~opts.ignoreCI
    %         tab{k}.ci_l(boutlier(:, k)) = nan;
    %         tab{k}.ci_h(boutlier(:, k)) = nan;
    %     end
    % end

    if ~opts.ignoreCI
        ci_l = cellfun(@(x)x.ci_l, tab, uni=false);
        ci_h = cellfun(@(x)x.ci_h, tab, uni=false);
        ci_l = vertcat(ci_l{:});
        ci_h = vertcat(ci_h{:});
        ci_l_b = ci_l;
        ci_h_b = ci_h;
        ci_l_b(biflag) = log(ci_l(biflag));
        ci_h_b(biflag) = log(ci_h(biflag));
        opts.xlim1 = isoutlier(ci_l_b, "mean", ThresholdFactor=2);
        opts.xlim2 = isoutlier(ci_h_b, "mean", ThresholdFactor=2);

        % make sure that new xlim doesn't cover betas
        xlim_b = [min(bvec_raw), max(bvec_raw)];

        % some offset
        if isfield(opts, "squeezeOffset")
            oset = opts.squeezeOffset;
        else
            oset = abs(diff(xlim_b))*.1;
        end

        if all(biflag)
            if (xlim_b(1)-oset) < 0
                xlim_b1 = xlim_b(1);
            else
                xlim_b1 = xlim_b(1)-oset;
            end
            xlim_b = [xlim_b1, xlim_b(2)+oset]; % don't cover negative values
        else
            xlim_b = [xlim_b(1), xlim_b(2)+oset];
        end

        xlim = [min(ci_l(~opts.xlim1)), max(ci_h(~opts.xlim2))];
        opts.xlim = [min(xlim(1), xlim_b(1)), ...
            max(xlim(2), xlim_b(2))];

        % opts.xlim(1) = opts.xlim(1) - 0.1*abs(diff(opts.xlim));
        % opts.xlim(2) = opts.xlim(2) + 0.1*abs(diff(opts.xlim));

        opts.xlim_cl = opts.xlim(1);
        opts.xlim_ch = opts.xlim(2);
        opts.xlim1(ci_l >= opts.xlim_cl) = false;
        opts.xlim2(ci_h <= opts.xlim_ch) = false;
        
        opts.xlim1(ci_l < opts.xlim_cl) = true;
        opts.xlim2(ci_h > opts.xlim_ch) = true;

        % opts.xlim(1) = opts.xlim(1) - 0.01;%*abs(diff(opts.xlim));
        % opts.xlim(2) = opts.xlim(2) + 0.01;%*abs(diff(opts.xlim));
    
        opts.xlim1 = reshape(opts.xlim1, [], numel(tab));
        opts.xlim2 = reshape(opts.xlim2, [], numel(tab));
    
    end
end

% check group markers
if isfield(opts, "marker")

    if istable(opts.marker)
        % table of markers with first column being the group labels and
        % second column being the marker type     
        [~, idx] = ismember(opts.groupLabels, opts.marker.(1)); 
        opts.marker = opts.marker.(2)(idx(idx > 0));

        opts.marker = flip(opts.marker); %@09APR2025: a bug was fixed
    end

    if numel(opts.marker) < numel(opts.groupLabels)
        % use only the first marker
        opts.marker = string(opts.marker);
        opts.marker = repmat(opts.marker(1), numel(opts.groupLabels), 1);
    elseif numel(opts.marker) > numel(opts.groupLabels)
        opts.marker = string(opts.marker);
        opts.marker = opts.marker(1:numel(opts.groupLabels));
    end
    opts.marker = flip(opts.marker);
else
    opts.marker = repmat("O", numel(opts.groupLabels), 1);
end

%@27MARCH2025: xscale option
h.XAxis.Scale = opts.xscale;

% loop over tables 
for i = 1:numel(tab)

    if i == 1 % sort effect sizes only for first phenotype
        if opts.sortEffect
            tab{i} = sortrows(tab{i}, 'b');
        end
        [~, idx] = ismember(tab{i}.snp, snps);
        snps = snps(idx(idx > 0));

        if opts.squeeze && ~opts.ignoreCI
            opts.xlim1(:, i) = opts.xlim1(idx(idx > 0), i);
            opts.xlim2(:, i) = opts.xlim2(idx(idx > 0), i);
        end
    end
    [~, idx2] = ismember(snps, tab{i}.snp); idx2(idx2 < 1) = [];
    tab{i} = tab{i}(idx2, :);

    if opts.squeeze && ~opts.ignoreCI
        opts.xlim1(:, i) = opts.xlim1(idx2, i);
        opts.xlim2(:, i) = opts.xlim2(idx2, i);
    end

    if i > 1 % double check tables' consistency
        if ~all(tab{i}.snp == tab{i-1}.snp)
            error("effPlotter:mismatch between tables rownames!")
        end
    end
    
    if ~opts.merge
        if opts.showp && opts.showpmethod == "bubble"
            h = subplot(1, numel(tab), i);
        else
            h = nexttile(ti);
        end
    end

    hold(h, 'on')
    h.FontName = opts.fontname;
    h.XLabel.Interpreter = 'latex';
    if opts.xbold; h.XAxis.FontWeight = 'bold'; end
    h.XAxis.FontSize = opts.fontsize;

    % draw null effect line 
    if ~opts.ignoreCI && all(opts.binary{i})
        nullLine = 1;
        if opts.xlabel == ""
            h.XLabel.String = "\textbf{OR}";
        else
            h.XLabel.String = opts.xlabel;
            h.XLabel.Interpreter = "tex";
        end
        h.XAxis.Scale = 'log';
    else
        nullLine = 0;
        if opts.xlabel == ""
            h.XLabel.String = "\mbox{\boldmath {$\beta$}}";
        else
            h.XLabel.String = opts.xlabel;
            h.XLabel.Interpreter = "tex";
        end
    end
    h.XLabel.FontSize = opts.fontsize + 2;

    if isfield(opts, "xline")
        nullLine = opts.xline;
    end

    if opts.zeroEffect
        xline(h, nullLine, '--', 'LineWidth', 1.3, 'Color', [105 105 105]./256);
    end
    
    y = 1:height(tab{i});
    h.YTick = y;

    if ~opts.overlap
        if i == numel(tab)
            h.YLim = [0, y(end) + 1];
        end
        y = y + yvec(i) - 1;
    else
        h.YLim = [0, y(end) + 1]; % should it be 0.5 and y(end) + 0.5 instead?  
    end

    %@17JULY2024: for ycol2 we need y coordinates
    if isfield(opts, "ycol2")
        opts.y2{i, 1} = y';
        opts.y2lab{i, 1} = tab{i}.(opts.ycol2);
    end

    if ~opts.ignoreCI
        % draw line for 95% CI
        drawnow
        lax = line([tab{i}.ci_l, tab{i}.ci_h]',[y; y], 'color', 'k', ...
            'LineWidth', opts.ciLineWidth, "LineStyle", "-");
        for m = 1:numel(lax)
            lax(m).Tag = "ci" + i;
        end
    
        % @15DEC2022: check XLim so that it has a margin on both ends
        drawnow
        if opts.squeeze
            XLim = opts.xlim;
            % if nullLine == 1 && XLim(1) < 0
            %     XLim(1) = min(tab{i}.ci_l);
            % end

            if any(opts.xlim1(:, i))
                [x_l, y_l] = deal(nan(numel(y), 1));
                cl_idx = opts.xlim1(:, i);
                y_l(cl_idx) = y(cl_idx);
                x_l(cl_idx) = opts.xlim_cl; % max(opts.xlim_cl, tab{i}.ci_l(cl_idx));
                cax1 = line([x_l, x_l]', [y_l, y_l]', "Marker", "<", "MarkerFaceColor", "k", ...
                    'LineWidth', opts.ciLineWidth, "LineStyle", "none",...
                    "MarkerEdgeColor", "k", "MarkerSize", 4);
                for m = 1:numel(cax1), cax1(m).Tag = "cil" + i; end
            end

            if any(opts.xlim2(:, i))
                [x_h, y_h] = deal(nan(numel(y), 1));
                ch_idx = opts.xlim2(:, i);
                y_h(ch_idx) = y(ch_idx);
                x_h(ch_idx) = opts.xlim_ch; % min(opts.xlim_ch, tab{i}.ci_h(ch_idx));
                cax2 = line([x_h, x_h]', [y_h, y_h]', "Marker", ">", "MarkerFaceColor", "k", ...
                    'LineWidth', opts.ciLineWidth, "LineStyle", "none",...
                    "MarkerEdgeColor", "k", "MarkerSize", 4);
                for m = 1:numel(cax2), cax2(m).Tag = "cih" + i; end
            end

        else
            XLim = h.XLim;
            if XLim(1) >= min(tab{i}.ci_l)
                XLim(1) = min(tab{i}.ci_l) - 0.01*abs(diff(XLim));
                if nullLine == 1 && XLim(1) < 0
                    XLim(1) = min(tab{i}.ci_l);
                end
            end
        
            if XLim(2) <= max(tab{i}.ci_h)
                XLim(2) = max(tab{i}.ci_h) + 0.01*abs(diff(XLim));
            end
        end
    else
        drawnow
        XLim = h.XLim;
        beta_range = [min(tab{i}.b), max(tab{i}.b)];
        if XLim(1) >= beta_range(1)
            XLim(1) = beta_range(1) - 0.01*abs(diff(XLim));
        end

        % check if xscale is log, and make the XLim(1) strictly positive
        if opts.xscale == "log" && XLim(1) < 0
            XLim(1) = max(eps, beta_range(1));
        end
    
        if XLim(2) <= beta_range(2)
            XLim(2) = beta_range(2) + 0.01*abs(diff(XLim));
        end
    end

    h.XLim = XLim;
    
    if i > 1
        if ~opts.merge
            h.YTickLabel = [];
        end
    else
        h.YTickLabel = tab{i}.snp;
        h.YAxis.TickLabelInterpreter = 'tex';
        if opts.ybold; h.YAxis.FontWeight = 'bold'; end
        h.YAxis.FontSize = opts.fontsize + 1;
    end
    if opts.square
        drawnow
        xlim = h.XLim;
        axis(h, 'square');
        drawnow
        h.XLim = xlim;
    end

    % add title 
    if ~opts.merge && opts.addTitle
        h.Title.String = opts.groupLabels(i);
        h.Title.FontSize = opts.titleFontSize;
    elseif isfield(opts, 'Title') && i == numel(tab)
        h.Title.String = opts.Title;
        h.Title.FontSize = opts.titleFontSize;
    end
    
    % draw effect sizes
    if opts.p(i) && opts.showp
        sax = scatter(h, tab{i}.b, y, opts.markersize(i), 'filled',...
            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', opts.markeralpha);
    elseif opts.merge
        sax(i) = scatter(h, tab{i}.b, y, opts.markersize(i), 'filled',...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', opts.color(i, :), ...
            'MarkerFaceAlpha', opts.markeralpha, 'Marker', opts.marker(i));
    else
        sax = scatter(h, tab{i}.b, y, opts.markersize(i), 'filled',...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', opts.color(1, :), ...
            'MarkerFaceAlpha', opts.markeralpha);
    end

    % @30MAY2023: make non-significant hits transparent
    if opts.hide
        nsidx = tab{i}.p > opts.hideCutoffLower | tab{i}.p <= opts.hideCutoffUpper;
        hs = findobj(h, "Type", "Scatter");
        hs(1).MarkerFaceAlpha = "flat";
        hs(1).AlphaDataMapping = "none";
        hs(1).SizeData = opts.markersize(i).*ones(numel(nsidx), 1);
        hs(1).AlphaData = ones(numel(nsidx), 1).*opts.markeralpha;
        hs(1).AlphaData(nsidx) = opts.hidealpha;
        hs(1).MarkerEdgeAlpha = opts.hidealpha + 0.1;
        
        % hs(1).SizeData(nsidx) = hs(1).SizeData(nsidx).*.7;
        
        tags = ["ci" "cil", "cih"] + i;
        for j = 1:numel(tags)
            l = findobj(h, "Tag", tags(j));
            if ~isempty(l)
                l(~flip(nsidx)) = [];
                for k = 1:numel(l)
                    l(k).Color = (1-opts.hidealpha).*[1 1 1];
                    l(k).MarkerFaceColor = (1-opts.hidealpha).*[1 1 1];
                    l(k).MarkerEdgeColor = (1-opts.hidealpha).*[1 1 1];
                end
            end
        end

    end

    if opts.p(i) && opts.showp
        h.TitleHorizontalAlignment = "right";
        if opts.showpmethod == "color" % ----------------------------------
            sax.CData = -log10(tab{i}.p);
            if max(-log10(tab{i}.p)) > 100
                set(h, 'ColorScale', 'log');
            end
            leg = addPcbar(h, opts, i);
        elseif opts.showpmethod == "scatter" % ----------------------------
            sax.Visible = 'off';
%             pscaled = rescale(-log10(tab{i}.p), opts.markersize*0.75, opts.markersize*2);
%             scatter(h, tab{i}.b, y, pscaled, 'O', ...
%                 'filled', 'MarkerEdgeColor', 'k', ...
%                 'MarkerFaceColor', opts.color);
            
            % pscaled = rescale(-log10(tab{i}.p), opts.markersize, opts.markersize*2);
            pscaled = opts.markersize(i).*(-log10(tab{i}.p)).^0.2 + 3; % 3-> if p is 1
            pscaled = modifyMarkerSize(pscaled, tab{i}, h, 'minus'); % marker shouldn't cover all 95% CI line
            for j = 1:numel(pscaled)
                plot(h, tab{i}.b(j), y(j), 'MarkerFaceColor', opts.color,...
                    'MarkerEdgeColor', 'k', 'LineStyle', 'none',...
                    'Marker', 'O', 'MarkerSize', pscaled(j));
            end

            % for legend
            opts.pLeg = round(quantile(pscaled, [0 .5 1]));
            opts.pLegTag = string(round(quantile(-log10(tab{i}.p), [0 .5 1])));
            leg = modifyLeg(h, opts, i);
        else % bubbles ----------------------------------------------------
            sax.Visible = 'off';
            bubblechart(h, tab{i}.b, y, -log10(tab{i}.p), ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', opts.color, ...
                'MarkerFaceAlpha', 1);
            bubblesize(h, opts.bubbleSize)
            pLegTag = round(quantile(-log10(tab{i}.p), [0 0.5 1]));
            leg = bubblelegend(h, 'Style', 'horizontal', 'Color', '#BAD8E0', ...
                'LimitLabels', string(pLegTag), 'FontSize', opts.fontsize - opts.legFontOffset);
            leg.Title.String = "-log_{10} \it{p}";
            leg.TextColor = '#3F3F3F';
            drawnow
            pos = plotboxpos(h);
            leg.Position(1) = pos(1);
            leg.Position(2) = pos(2) + pos(4);
        end
        
        % shift title upwards to avoid colorbar/legend overlapping
        h.Tag = "ax" + i;
        h = checkTitle(h, leg); % make sure that title is within current axis limits
        
    end
    
    if opts.grid
        h.YGrid = 'on';
        h.XGrid = 'on';
    end

    if opts.box, h.Box = 'on'; end
    if ~opts.merge; hold(h, 'off'); end
end


if opts.merge
    drawnow    
    h.YLim = [h.YLim(1)+0.5, h.YLim(2)-0.5];

    %@17JULY2024
    if isfield(opts, "ycol2")
        yyaxis(h, "right");
        h.YAxis(2).Color = "k";
        h.YAxis(2).FontSize = opts.fontsizeYcol2;
        h.YAxis(2).FontName = opts.fontname;
        h.YAxis(2).TickDirection = "out";
        h.YAxis(2).Limits = h.YAxis(1).Limits;
        y2 = vertcat(opts.y2{:});
        y2lab = vertcat(opts.y2lab{:});
        [y2, idx] = sort(y2);
        y2lab = y2lab(idx);
        h.YAxis(2).TickValues = y2;
        h.YAxis(2).TickLabels = y2lab;
        
        %@08OCT2024: add label
        h.YAxis(2).Label.String = opts.ycol2Label;
        h.YAxis(2).Label.Rotation = 0;
        drawnow;
        un = h.YAxis(2).Label.Units;
        h.YAxis(2).Label.Units = "normalized";
        h.YAxis(2).Label.Position(2) = 1;
        h.YAxis(2).Label.Position(1) = 1.1; % h.YAxis(2).Label.Position(1) - 2e-3;
        h.YAxis(2).Label.VerticalAlignment = "bottom";
        h.YAxis(2).Label.HorizontalAlignment = "right";
        
        drawnow;
        h.YAxis(2).Label.Units = un;
        yyaxis(h, "left");
    end
    
    h = shade(h, 'length', 1, "FaceAlpha", 0.1);

    if opts.addTitle % && numel(tab) > 1
        leg = legend(sax);
        leg.String = opts.groupLabels;
        if isfield(opts, 'legTitle')
            leg.Title.String = opts.legTitle;
        end
        leg.Orientation = opts.legendOrientation;
        leg.NumColumns = opts.legendCols;
        leg.FontSize = opts.fontsize - opts.legFontOffset;
        leg.Location = opts.legLocation;
        leg.AutoUpdate = false;
        drawnow;

        if any(contains(lower(leg.Location), "eastoutside"))            
            pos = h.tightPosition; % pos after axis square
            while (pos(1) + pos(3)) > leg.Position(1)
                leg.Position(1) = pos(1) + pos(3) + 0.0025;
                drawnow;
                pos = h.tightPosition;
            end
        end

        if leg.Position(1) < 0
            leg.Orientation = 'vertical';
            leg.Location = 'bestoutside';
        end
    end
    hold(h, 'off');
elseif numel(tab) < 2 || opts.shade
    drawnow
    h.YLim = [h.YLim(1)+0.5, h.YLim(2)-0.5];
    h = shade(h, 'length', 1);
end

% update colorbar positions (tend to change with 'tight' in tiledlayout)
if ~(opts.showp && opts.showpmethod == "bubble") && ~opts.merge
    drawnow
    ti = updatePCbar(ti);
else

    % under development: for meta-analysis
    % h = nexttile(ti);
    % hold(h, "on")
    % for k = 1:numel(opts.groupLabels)
    %     xx = k.*ones(numel(tab{k}), 1);
    %     yy = (1:numel(tab{k}))';
    %     idx = tab{k}.b < 0;
    %     scatter(h, xx(idx), yy(idx), Marker="_", MarkerFaceColor="k", MarkerEdgeColor="k");
    %     scatter(h, xx(~idx), yy(~idx), Marker="+", MarkerFaceColor="k", MarkerEdgeColor="k");
    % end
    % hold(h, "off")
    % h.XTick = 1:numel(opts.groupLabels);
    % h.XTickLabel = opts.groupLabels;
    % h.YAxis.Visible = "off";

    ti = h;
end

fig.Visible = 'on';
drawnow

if opts.save

    for k = 1:numel(opts.format)
        if any(ismember(opts.format(k).lower, ["pdf", "eps"]))
            exportgraphics(ti.Parent, opts.output + "." + opts.format(k), ...
                "ContentType", "vector")
        else
            exportgraphics(ti.Parent, opts.output + "." + opts.format(k), ...
                'Resolution', opts.resolution)
        end
    end
end

if opts.savefig % doesn't work with "scatter"
    savefig(ti.Parent, opts.output + ".fig")
end

% warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode')
end % END

%% subfunctions ===========================================================
function c = addPcbar(h, opts, cnt)
c = colorbar(h, 'Location', 'north');
try
    cc = get_colormap(opts.colormap, 256);
    colormap(h, cc);
catch
    colormap(h, opts.colormap)
end
drawnow
ticks = c.Ticks;
% ticklabs = c.TickLabels;
c.Position(1) = h.Position(1);
c.Position(2) = h.Position(2) + h.Position(4);% - c.Position(4)*.8;
c.Position(4) = c.Position(4)*.8;
c.Position(3) = c.Position(3)/opts.cbarSqueeze;

c.AxisLocation  = 'out';
c.Label.String = "-log_{10} \it{p}";
c.Label.Units = 'normalized';
c.Label.HorizontalAlignment = 'left';
c.Label.VerticalAlignment = 'middle';
c.FontSize = opts.fontsize - opts.legFontOffset - 4;
c.Label.FontSize = opts.fontsize - opts.legFontOffset;
c.Label.Position = [1 0.5 0];
if numel(ticks) > 3
    if all(ticks < 1)
        dig = 2;
    else
        dig = 1;
    end
    ticks = [ticks(1), round(median(ticks), dig, "significant"), ticks(end)];
%     ticklabs = string(ticks);
end
ticks = unique(ticks);
c.Ticks = ticks;
c.TickLabels = string(ticks);
c.Tag = "ax" + cnt;

while c.Label.Extent(2) < 0
    if opts.adjmethod == "pos"
        c.Label.Position(2) = c.Label.Position(2) + 3e-2;
    else
        c.Label.FontSize = c.Label.FontSize - opts.legFontOffset;
    end
    drawnow
end
end

%% ------------------------------------------------------------------------
function leg = modifyLeg(h, opts, cnt)
lg = (0);
for j = 1:numel(opts.pLeg)
%     lg(j) = scatter(h, nan, nan, opts.pLeg(j), 'O', ...
%                 'filled', 'MarkerEdgeColor', 'k', ...
%                 'MarkerFaceColor', opts.color);
    lg(j) = plot(h, nan, nan, 'MarkerFaceColor', opts.color,...
        'MarkerEdgeColor', 'k', 'LineStyle', 'none',...
        'Marker', 'O', 'MarkerSize', opts.pLeg(j));
end

% drawnow
leg = legend(lg);
leg.String = opts.pLegTag;
leg.Location = 'north';
leg.Color = "#BAD8E0";
leg.Orientation = 'horizontal';
leg.FontSize = opts.fontsize - 4;
leg.Title.FontSize = leg.FontSize - opts.legFontOffset;
title(leg, "-log_{10} \it{p}")

legEntry = leg.EntryContainer.Children;
for j = 1:numel(legEntry)
    drawnow
    legEntry(j).Label.NodeChildren.HorizontalAlignment = 'center';
end

% % Doesn't need since scatter was replaced by plot
% % pLeg = sort(rescale(opts.pLeg, 5, 10), 'descend');
% for j = 1:numel(legEntry)
%     drawnow
%     legEntry(j).Icon.Transform.Children.Children.Size = (opts.pLeg(j));
% end

% DEPRECATED ##############################################################
% textmark = findobj(legmark, 'Type', 'Text');
% legmark = findobj(legmark, 'Type', 'Line');
% 
% for j = 2:2:numel(legmark)
%     legmark(j).MarkerSize = pLeg(j/2);
% end
% for i = 1:numel(textmark)
%     textmark(i).HorizontalAlignment = 'center';
% end
% % title doesn't work well when legend nargout > 1
% te = text(1, 1, "-log_{10} \it{p}", 'Units', 'normalized',...
%     'Parent', leg.DecorationContainer, 'VerticalAlignment', 'cap', ...
%     'Interpreter', 'tex', 'FontSize', opts.fontsize - 2);
% #########################################################################

leg.TextColor = '#3F3F3F';
leg.Units = 'normalized';
drawnow
leg.Position(1) = h.Position(1);
leg.Position(2) = h.Position(2) + h.Position(4);
leg.Position(3) = leg.Position(3)*.9;
if leg.Position(3) > h.Position(3)*0.45 % if too wide
    drawnow
    leg.Position(3) = h.Position(3)*0.45;
end
leg.Tag = "ax" + cnt;

% DEPRECATED ##############################################################
% while te.Extent(2) < 0
%     if opts.adjmethod == "pos"
%         te.Position(2) = te.Position(2) + 3e-2;
%     else
%         te.FontSize = te.FontSize - 1;
%     end
%     drawnow
% end
% #########################################################################

end

%% ------------------------------------------------------------------------
function ti = updatePCbar(ti)

c = findobj(ti.Children, 'Type', 'Colorbar');
if isempty(c)
    c = findobj(ti.Children, 'Type', 'Legend');
    if isempty(c); return; end
end
h = findobj(ti.Children, 'Type', 'Axes');
ctag = string({c.Tag});

for i = 1:numel(h)
    ax = h(i);
    idx = ax.Tag == ctag;
    cbar = c(idx);
    drawnow

    cbar.Position(1) = ax.Position(1);
    cbar.Position(2) = ax.Position(2) + ax.Position(4);
end
end

%% ------------------------------------------------------------------------
function h = checkTitle(h, leg)
% makes sure that title is within current axis limits (when legend/colorbar
% pushes title away)

h.Title.Units = 'normalized';
if string(h.Title.String) == "" || isempty(h.Title.String)
    return % no title exists
end

drawnow
if strcmp(leg.Type, 'bubblelegend')
    pos = plotboxpos(h);
else
    pos = h.Position;
end

% if ~strcmp(leg.Type, {'bubblelegend', 'legend'}) % bubblelegend is compact enough
%     h.Title.Position(2) = leg.Position(4)*opts.titleBuffer + 1;
% end 

% ratio = h.PlotBoxAspectRatio;

% check vertical overlap
% while (pos(2) +pos(4) + h.Title.Position(2) - 1) > 0.98
%     drawnow
%     ratio(2) = ratio(2).*0.99;
%     h.PlotBoxAspectRatio = ratio;
%     if strcmp(leg.Type, 'bubblelegend')
%         pos = plotboxpos(h);
%     else
%         pos = h.Position;
%     end
% end

while (pos(2) +pos(4) + h.Title.Extent(2) + h.Title.Extent(4) - 1) > 0.98
    drawnow
    h.Title.Position(2) = h.Title.Position(2).*.9999;
    if strcmp(leg.Type, 'bubblelegend')
        pos = plotboxpos(h);
    else
        pos = h.Position;
    end
end

% check horizontal overlap
while (pos(1) + pos(3)*h.Title.Extent(1)) < leg.Position(1) + leg.Position(3)
    drawnow
    h.Title.Position(1) = h.Title.Position(1) + 1e-2;
end

end

%% ------------------------------------------------------------------------
function pscaled = modifyMarkerSize(pscaled, tab, h, method)
% plot MarkerSize 72 fits into an axis of 1 inch wide
if strcmp(h.Parent.Parent.Visible, 'off')
    h.Parent.Parent.Visible = 'on';
    wasoff = true;
else
    wasoff = false;
end
drawnow
h.Units = 'inches';
pos = h.Position;
drawnow
h.Units = 'normalized';
if wasoff; h.Parent.Parent.Visible = 'off'; end
ci = [tab.ci_l, tab.ci_h];
ci = diff(ci, [], 2)/diff(h.XLim); % CI width of each marker

p = (pscaled + 2)./72; % +2 for MarkerEdge
p = p./min(pos(3:4)); % diameter of each marker

while any(p >= ci)
    if strcmp(method, 'minus')
        pscaled = pscaled - 0.1;
    elseif strcmp(method, 'rescale')
        pscaled = rescale(tab.p, min(pscaled), max(pscaled)*0.99);
    else
        error("effPlotter:modifyMarkerSize", "unknown method!")
    end

    p = (pscaled + 2)./72; % +2 for MarkerEdge
    p = p./min(pos(3:4)); % diameter of each marker
end

if any(pscaled < 0)
    pscaled = rescale(pscaled, 1, max(pscaled));
end

end

%% ------------------------------------------------------------------------
function mustBeAllTab(tab)
if istable(tab), tab = {tab}; end
if ~all(cellfun(@(x) istable(x), tab))
    eid = 'input:notTable';
    msg = 'first argument must either be a cell with all elements being a table or a single table';
    throwAsCaller(MException(eid,msg))
end
end