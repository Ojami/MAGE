function clumped = manplotter(infile, opts)
% description goes here %
% Oveis Jamialahmadi, July 2020. Sahlgrenska Akademy, GU.
% 
% @30/12/2021: some bugs were fixed with Ybreak.
% @11/01/2022: a bug was fixed.
% @13/06/2022: when using a non-clumped GWAS summary stat file, and 'label'
%              is true, top loci are detected via findpeaks built-in
%              function, and then annotated to find nearest genes. 
% @11FEB2023: two clumping methods are now supported: PLINK or internal
%             rough physical clumping using signal toolbox by identifying
%             local maxima. Default is PLINK clumping. Note that, for
%             'peak', use a small window for 'kb' argument, as windows for
%             'peak' method is different from that of the 'plink' method.
% @21MAR2023: 'magma' flag was added. If true, input file will be treated
%              as a MAGMA gene-based result table. Gene positions will be
%              calculated as mean of START and STOP columns.

arguments
    infile {mustBeFile}  % the input summary stat file
    opts.parallel (1,1) logical = false % use tall/gather for reading 'infile' summary stat file
    opts.col (1,1) {mustBeMember(opts.col, ["b", "m"])}= "b" % color pattern: 'b' for binary and 'm' for multiple
    opts.break (1,1) double {mustBeNonNan, mustBeNonnegative} =  20 % Y-axis break
    
    % Y-axis damp, the number can be greated than 1, when max -log10 p-value is
    % less than twice of breaking point (e.g. strongest association is 27, but
    % breaking point is 20, so if Ydamp is 4, it shrinks -log10 p-values in
    % range 20-27 four times more compared to -log10 p-values below 20).
    opts.Ydamp (1,1) double {mustBeGreaterThanOrEqual(opts.Ydamp, 1)} = 1 
    opts.maf (1,1) double {mustBeInRange(opts.maf, 0, 0.5)} = 0
    opts.title {mustBeTextScalar} = "Manhattan-manplotter" 
    opts.light (1,1) logical = true % light version of Manhattan plot by excluding variants with P < 0.05
    opts.sig (1,1) double {mustBeInRange(opts.sig, 0, 1)} = 5e-8 % significance threshold
    opts.label (1,1) logical = false % label hits. If 'clump' or 'gene' column are present, this can be set to true.
    opts.maxlabel (1,1) double = 20 % max number of labels to show (if opts.label is true)
    opts.chrSep (1,1) logical = false % draw dashed vertical lines to separate chromosomes 
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "tif" % format of output plots
    opts.resolution (1, 1) double = 400 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    
    % aesthetics
    opts.markerSize (1,1) double = NaN 
    opts.edgeAlpha (1,1) double = NaN
    opts.markerAlpha (1,1) double = NaN
    opts.markerAlpha_sig (1,1) double = NaN
    opts.fontsize (1,1) double = 14
    opts.fontname {mustBeTextScalar} = "Garamond"
    opts.xfontsize (1,1) double = 16
    opts.yfontsize (1,1) double = 17
    opts.box (1,1) {mustBeMember(opts.box, ["on", "off"])} = "off" % enclose the plot within a box
    opts.magma (1,1) logical = false % treat the input as MAGMA gene-based summ stat file

    % clumping options
    opts.bgenhome {mustBeFolder} 
    opts.clumpingMethod {mustBeMember(opts.clumpingMethod, ["plink", "peak"])} = "plink" % clumping method: only if 'label' is true
    opts.clumpingFolder {mustBeFolder} % bed/bim/fam files for clumping should be present here
    opts.kb (1,1) double = 1000 % for peak use 1e2 or 1e3 depending on the number of independent loci.
    opts.r2 (1,1) double = 0.01

    %@01NOV2024
    opts.refGenome {mustBeTextScalar, mustBeMember(opts.refGenome, ["GRCh37", "GRCh38"])} =  "GRCh37" % mapVariant2gene argument
    opts.checkld (1,1) logical = true % effective if clumpingMethod is peak, if false (doesn't prune based on LD friends).
end

fclose('all');
close all

% read GWAS summary stat file ---------------------------------------------
if opts.label && opts.clumpingMethod == "plink"
    opts.fullgwas = true; % PLINK needs this
else
    opts.fullgwas = false;
end

infile = readGWASfile(infile, 'parallel', opts.parallel, ...
    'full', opts.fullgwas, ...
    'light', opts.light, 'p', 0.05);
infile(ismissing(infile.p), :) = [];

% reconcile zero p-values passing machine precision
f_zero = infile.p == 0;
if any(f_zero)
    infile.p(f_zero) = realmin*eps;
    infFlag = true;
else
    infFlag = false;
end

if any(infile.p <= 0) || any(infile.p > 1)
    error('P values must be in [0, 1)!')
end

infile.p = -log10(infile.p); % log10 transofrmation of p-values
opts.sig = -log10(opts.sig);

% check if clump column exists: if so, only annotated labels for strongest
% association per each clump is shown. If clump column does not exist (e.g.
% raw input from PheWeb or other large-scale statistical tools such as
% PLINK, SAIGE,...), manplotter checks for duplicates in annotate column
% (note that this only works for gene names and not rs ids) and keeps the
% label with strongest association (e.g. one label for each gene). Only
% labels passing defaultSig are shown.
% @06/12/2022: now finds local maxima in a margin of 1 mbp around each
% peak and annotates variants to find the nearest gene (locus). 

colNames = infile.Properties.VariableNames;
opts.snp = false;

if opts.magma
    % create 'pos' column
    infile.pos = (infile.start + infile.stop)./2;
    infile(:, ["start", "stop"]) = [];
    infile.annotate = infile.symbol;
end

% sort chromosomes
infile = sortrows(infile, "pos", "ascend"); % findpeaks needs strictly increasing points
if isstring(infile.chr)
    [~, idx] = natsort(infile.chr);
    infile = infile(idx, :);
else
    infile = sortrows(infile, "chr", "ascend");
end

if any(ismember(colNames, 'clump'))
    infile.annotate(~logical(infile.clump)) = "";
elseif opts.label && ~opts.magma % no clump column
    
    if ~any(infile.p >= opts.sig)
        opts.label = false;
    else
        opts.snp = true; % no clump/gene
        
        if opts.clumpingMethod == "peak"
            infile = ldclumpPeaks(infile, sig=opts.sig, win=opts.kb, r2=opts.r2, ...
                parallel=opts.parallel, log10p=true, bgenhome=opts.bgenhome, ...
                checkld=opts.checkld);
            infile = renamevars(infile, ["annotate", "locus"], ["snp", "annotate"]);

        elseif opts.clumpingMethod == "plink"
            bedfiles = getfilenames(opts.clumpingFolder, "bed", "fullpath", true).bed;
            if isempty(bedfiles) % bgen?
                bedfiles = getfilenames(opts.clumpingFolder, "bgen", "fullpath", true).bgen;
            end
            % bedfiles = fullfile(opts.clumpingFolder, natsort(bedfiles));
            clumped = ldclump(infile, bedfiles, parallel=opts.parallel, ...
                r2=opts.r2, kb=opts.kb, save=false, log10p=true, p1=10.^-opts.sig);
            vep = mapVariant2gene(clumped, verbose=true, ...
                refGenome=opts.refGenome, bgenhome=opts.bgenhome);
            [f1, f2] = ismember(clumped.SNP, vep.id); f2(f2<1) = [];
            clumped.locus(:) = "";
            clumped.locus(f1) = vep.nearestGene(f2);
            
            infile.snp = infile.annotate;
            infile.annotate(:) = "";
            [f1, f2] = ismember(infile.snp, clumped.SNP); f2(f2<1) = [];
            infile.annotate(f1) = clumped.locus(f2);
        end

    end
end

% return clumped variants 
if opts.label && ~opts.magma
    infileTopLoci = infile(infile.annotate ~= "", :);
    infile(infile.annotate ~= "", :) = [];

    if opts.clumpingMethod == "peak"
        clumped = infileTopLoci;  
        clumped = renamevars(clumped, "annotate", "locus");
    end

    % only keep 1 gene name for each peak
    infileTopLoci = groupfilter(infileTopLoci, "annotate", @(x)x == max(x), "p");
    infile = [infile; infileTopLoci];
else
    clumped = [];
end

infile = sortrows(infile, 'p', 'descend');
if opts.label
    f_annotate_keep = find(~cellfun(@isempty,infile.annotate));
    try
        infile.annotate(f_annotate_keep(opts.maxlabel + 1):end) = "";
    catch
    end
end

if isstring(infile.chr)
    chrU = natsort(unique(infile.chr));
else
    chrU = unique(infile.chr);
end

if strcmp(opts.col,'m')
    [colMap, colMap_sig] = deal(lines(numel(chrU))); % color codes for each chr
    if isnan(opts.markerAlpha)
        opts.markerAlpha = 0.25;
    end
    if isnan(opts.markerAlpha_sig)
        opts.markerAlpha_sig = 0.6;
    end
elseif strcmp(opts.col,'b')
    colCode = [0 0 0;128 128 128]./256;
    colColde_sig = [0,0,139;0,0,205]./256;
    colMap = repmat(colCode, round(numel(chrU)/2), 1);
    colMap_sig = repmat(colColde_sig, round(numel(chrU)/2), 1);
    if isnan(opts.markerAlpha)
        opts.markerAlpha = 0.5;
    end
    if isnan(opts.markerAlpha_sig)
        opts.markerAlpha_sig = 0.5;
    end
end

% Define settings
if size(infile, 1) > 1e6
    if isnan(opts.markerSize); opts.markerSize = 15; end
    if isnan(opts.edgeAlpha); opts.edgeAlpha = 0.025; end
else
    if isnan(opts.markerSize); opts.markerSize = 20; end
    if isnan(opts.edgeAlpha); opts.edgeAlpha = 0.2; end
end

offsetCoeff = (max(infile.pos) - min(infile.pos)) / (numel(chrU));
maxY = 1.2*max(infile.p);
breakMe = opts.break; 

infile.p_raw = infile.p; % for data tips

% zoom-out for very low p-values (-log10 > breakMe)
f_break = infile.p > breakMe; % zoom-out for very low p-values
if any(f_break) %%&& (max(infile.p) > 2*breakMe)
    maxOldP = max(infile.p); % to maintain Y-axis upper limit close to maximum hit -log_p 
    breakFlag = true;
    
    % <---Ybefore--->  <----Yafter--->
    % Ylim(1)-----|breakMe|-----Ylim(2) 
    % main idea: length of Ybefore should be equal to Yafter
    Y.after = [min(infile.p(f_break)), max(infile.p(f_break))];
    Y.after(1) = breakMe; % strictly set the first point after break to breaking point
    Y.before = [min(infile.p(~f_break)), max(infile.p(~f_break))];
    if sum(f_break) < 2 % only one hit passes breaking point
        Y.after(1) = (breakMe + Y.after(2))./2;
    end

    Y.beta = (diff(Y.after)/diff(Y.before)).*opts.Ydamp; % scaling coeff.
    if Y.after(2) > 300
        step = 100;
    elseif Y.after(2) > 200
        step = 40;
    elseif Y.after(2) > 100
        step = 20;
    elseif Y.after(2) > 40
        step = 10;
    else
       step = 5;
    end
    
    if infFlag % shift downwards the last point in Y.after to avoid overlap with 'inf' symbol on the plot
        Y.after(2) = Y.after(2) - 10;
    end

    breakYLim = ceil(unique([Y.after(1):step:Y.after(2), Y.after(2)]));

    if breakYLim(1) < Y.after(1); breakYLim(1) = ceil(Y.after(1)) + eps; end
    breakYLim = [breakYLim; (breakYLim-Y.after(1))./Y.beta + breakMe];
    infile.p(f_break) = (infile.p(f_break)-Y.after(1))./Y.beta + breakMe;

    % cap Y-axis upper limit
    idx = find(breakYLim(1, :) > maxOldP & breakYLim(2, :) > max(infile.p));
    if numel(idx) > 1
        idx(1) = [];
    end
    breakYLim(:, idx) = [];
    
else
    breakFlag = false;
end

% Create a figure
g = figure('Visible', 'on');
h = gca;

hold on
xoffset = 0; getpos = (0); maxX = (0); 
maxlabelpos = -1; % max label on Y-axis (label of the strongest hit)

infile.sig = infile.p > opts.sig;
infile.pos_plot = infile.pos; % for drawing 

for ii = 1:numel(chrU) % loop over chromosomes and plot -log10 P vs pos
    fidx = infile.chr == chrU(ii);
    infile.pos_plot(fidx) = infile.pos(fidx) + xoffset; % Update pos
    
    if ii == 1 % for first chromosome
        Xoffset = .5*median(infile.pos_plot(fidx));
        infile.pos_plot(fidx) = infile.pos_plot(fidx) + Xoffset; % shift towards right
    end

    maxX(ii) = max(infile.pos_plot((fidx)));

    scatter(h, infile.pos_plot(~infile.sig & fidx), infile.p(~infile.sig & fidx), 'MarkerFaceColor', colMap(ii,:),...
        'MarkerEdgeColor', 'k', 'SizeData', opts.markerSize,...
        'MarkeredgeAlpha', opts.edgeAlpha, 'MarkerFaceAlpha', opts.markerAlpha)

    ax = scatter(h, infile.pos_plot(infile.sig & fidx), infile.p(infile.sig & fidx), 'MarkerFaceColor', colMap_sig(ii, :),...
        'MarkerEdgeColor', 'k', 'SizeData', opts.markerSize,...
        'MarkeredgeAlpha', opts.edgeAlpha, 'MarkerFaceAlpha', opts.markerAlpha_sig);

      top_hits = infile.annotate(infile.sig & fidx);
      f_empty = cellfun(@(x) isempty(x) | strcmp(x, ''), top_hits);
    
    if ~isempty(top_hits)
        sig_pos = infile.pos_plot(infile.sig & fidx);
        sig_p = infile.p(infile.sig & fidx);

        % add data-tip
        if opts.snp
            dtrow = dataTipTextRow('snp', infile.snp(infile.sig & fidx));            
        else
            dtrow = dataTipTextRow('annotation', top_hits);
        end
        if numel(dtrow.Value) < 2, dtrow.Value = [dtrow.Value, dtrow.Value]; end
        ax.DataTipTemplate.DataTipRows = matlab.graphics.datatip.DataTipTextRow;
        ax.DataTipTemplate.DataTipRows(1) = dtrow;
        dtrow = dataTipTextRow('p', 10.^-infile.p_raw(infile.sig & fidx));
        ax.DataTipTemplate.DataTipRows(2) = dtrow;

        if opts.label
            te = text(h, sig_pos(~f_empty), sig_p(~f_empty),...
                top_hits(~f_empty), 'FontAngle', 'italic',...
                'HorizontalAlignment', 'center', 'FontSize', opts.fontsize-2,...
                'VerticalAlignment', 'bottom', 'Margin', 10, ...
                'FontName', opts.fontname);
            tepos = [te.Position];
            [tepos, temaxIdx] = max(tepos(2:3:end));
            if tepos > maxlabelpos
                maxlabelpos = te(temaxIdx).Extent(2) + te(temaxIdx).Extent(4);
            end
        end

    end
    
    getpos(ii) = (max(infile.pos_plot(fidx)) + min(infile.pos_plot(fidx)))/2;
    xoffset = max(infile.pos_plot(fidx)) + offsetCoeff; % Set new starting point
    
    if opts.chrSep
        if ii < numel(chrU)
            line([xoffset, xoffset], [0, maxY], 'LineStyle','--', 'color', 'c')
        end
    end
end


line(h.XLim, [opts.sig, opts.sig], 'LineWidth', 1, 'color', 'r',...
    'LineStyle', '--') % Significance threshold
hold off
h.XTick = getpos;

h.XTickLabel = strtrim(cellstr(num2str(chrU)));
set(h,'box',opts.box,'LineWidth', 1.2, 'TickLength', [0.001 0],...
    'FontSize', opts.fontsize, 'FontWeight', 'bold', 'TickDir', 'out')
xlabel('Chromosome', 'FontSize', opts.xfontsize, 'FontWeight', 'bold')
ylabel('-log_{10} (\itp)', 'FontSize', opts.yfontsize, 'FontWeight', 'bold')
h.XAxis.TickLength = [0 0];
h.YLim(1) = min(infile.p) - 0.5;
h.XLim(2) = max(max(maxX)) + Xoffset;
if breakFlag
    h.YLim(2) = breakMe*2 + 2; %maxY
else
    h.YLim(2) = max(infile.p) + 2;
end

set(g, 'Visible', 'on')
set(gcf, 'WindowState', 'maximized');

% modify YTickLabels and add annotations
if exist('breakYLim', 'var')
    drawnow;
    oldYTick = h.YTick; oldYTick(oldYTick > breakMe) = [];
    oldYTick(end) = [];
    h.YTick = [oldYTick, breakYLim(2,:)];
    h.YLim(2) = h.YTick(end) + 2;
    h.YTickLabel = [string(oldYTick), string(breakYLim(1, :))];
end

if opts.label && opts.box == "on"
    if h.YLim(2) < maxlabelpos % axis doesn't enclose the strongest hit's label
        h.YLim(2) = maxlabelpos;
    end
end

if infFlag % breakFlag must be true, right?
    drawnow;
    h.YTick(end + 1) = max(infile.p(f_break));
    h.YTickLabel(end + 1) = {'Inf'};
    drawnow;
    H = h.Position(4)*(max(infile.p(f_break))-h.YLim(1))/diff(h.YLim)...
        + h.Position(2);
    drawnow;
    annotation(gcf, 'line', [h.Position(1)-4.5e-3, h.Position(1)+4.5e-3],...
        [H-4.5e-3, H+4.5e-3], 'LineWidth', 2, 'color', [0,102,102]./256)
end

if breakFlag
    drawnow;
    H = h.Position(4)*(breakMe-h.YLim(1))/diff(h.YLim) + h.Position(2);
    drawnow;
    annotation(gcf, 'line', [h.Position(1)-4.5e-3, h.Position(1)+4.5e-3],...
        [H-4.5e-3, H+4.5e-3], 'LineWidth', 2, 'color', [0,102,102]./256)
end
h.FontName = opts.fontname;

% check if text objects pass Ylim
te = findobj(h, "Type", "Text");
for k = 1:numel(te)
    if te(k).Position(2) > h.YLim(2)
        te(k).Position(2) = h.YLim(2);
    end
end

if opts.save
    exportgraphics(h, opts.title + "." + opts.format, 'Resolution', opts.resolution)
end

if opts.savefig
    savefig(h.Parent, opts.title + ".fig")
end

end

