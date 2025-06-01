function ax2 = rarePlotter(restab, opts)
% geneEffPlotter is similar to genePlotter but instead of p-values,
% visualizes effect sizes and their 95% CI on a region of interest (it must
% contain data for one chromosome).
% INPUTS:
%   - restab: a table with similar columns from gwasrunner, and necessary
%             columns are "95% CI", "POS", "CHR", "OR" or
%             "beta"/"β"/"BETA". if "consequence" is present, markers will
%             be classified in different groups accordingly.
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, Aug 2021.
% 
% @08/27/2021 a bug was fixed in reading "OR" column, and inf
%                      values in CI were handled. 

arguments
    restab {mustBeA(restab, 'table')} % a table with same columns as the output of gwasrunner function
    opts.refGenome {mustBeMember(opts.refGenome, ["GRCh37", "GRCh38"])} = "GRCh38"
    opts.raremaf (1,1) double {mustBeInRange(opts.raremaf, 0, 0.5)} = 0.5 % variants >= this MAF will not be shown.
    opts.sz (1,1) double = 55 % marker size
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "png" % format of output plots
    opts.resolution (1, 1) double = 400 % resolution of output plot
    opts.save (1,1) logical = true % only if plot is true
    opts.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    opts.output {mustBeTextScalar} = ""
    opts.verbose (1,1) logical = true
    opts.shade (1,1) logical = true % shade exonic regions on effect plot

    %@12OCT2024: column name can be set
    opts.legacy (1,1) logical = true % old implementation which works with gwasrunner (should be changed in future)

    % Some of the following flags work only with 'legacy' set to false
    opts.mafCol {mustBeTextScalar} = "A2Freq";
    opts.pcol {mustBeTextScalar} = "P-adj";
    opts.CIcol {mustBeTextScalar} % CI column name, if left empty, tries to find, if doesn't exist, will generate CI internally.
    opts.consequence {mustBeTextScalar} % column name of variant consequences (e.g. LoF, missense, ...)
    opts.binary (1,1) logical = false % are estimates in 'betacol' OR or beta? default is false (beta or log OR).
    opts.betacol {mustBeTextScalar} = "BETA" % column name of estimates, if value under 'binary' column is tru, log of this value will be used.
    opts.secol {mustBeTextScalar} = "SE"
    opts.poscol {mustBeTextScalar} = "POS"
    opts.chrcol {mustBeTextScalar} = "CHR"
    opts.title {mustBeTextScalar} % name of phenotype
    opts.ylabel = "Beta"
    opts.padjust {mustBeMember(opts.padjust, ["bonferroni", "none", "BH"])} = "none" % adjustment for number of tests
    opts.pcutoff (1,1) double = 0.05 % variants with an adjusted p using 'padjust' < 'pcutoff' will be hidden.
    opts.colormap {mustBeMember(opts.colormap, ["jet", "turbo", "abyss", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", "summer", ...
        "colorcube", "lines", "magma", "plasma", "inferno", "cividis", ...
        "viridis"])} = "plasma" % colormap for when 'color' is not used.
    opts.method {mustBeMember(opts.method, ["scatter", "bubble"])} = "bubble" % if bubble (default) also shows allele freqnecies of variants
    opts.bsize (1, 2) double = [5, 8] % bubble size range for AAF if plot method is "bubble"
    opts.hidealpha (1,1) double {mustBeInRange(opts.hidealpha, 1e-7, 1)} = 0.1
end

close all force

if opts.refGenome == "GRCh38"
    opts.dataset = "gnomad_r3";
elseif opts.refGenome == "GRCh37"
    opts.dataset = "gnomad_r2_1";
end

restab(restab.(opts.mafCol) >= opts.raremaf, :) = [];

% find beta column and remove missing values
if opts.legacy
    restab.Properties.VariableNames = replace(colnames(restab), "(Firth)", "");
    if any(strcmpi('or', restab.Properties.VariableNames))
        restab(isnan(restab.OR), :) = [];
    else
        betaCol = find(ismember(lower(restab.Properties.VariableNames), {'beta', 'β'}));
        restab(isnan(restab.(betaCol)), :) = [];
    end

else
    if isfield(opts, "betacol")
        betaCol = opts.betacol;
    else
        betaCol = find(ismember(lower(restab.Properties.VariableNames), {'beta', 'β'}));
    end
    restab(isnan(restab.(betaCol)), :) = [];

end

pcol = opts.pcol;
restab = renamevars(restab, [opts.chrcol, opts.poscol], ["CHR", "POS"]);

% extract CI
if ~isfield(opts, "CIcol")
    CIcol = find(ismember(lower(restab.Properties.VariableNames), {'95% ci', 'ci'}));
else
    CIcol = opts.CIcol;
end

if isempty(CIcol) % failed to find: calculate
    ciflag = false;
else
    ciflag = true;
end

if opts.legacy % CI must be present
    ci = restab.(CIcol);
    ci = cellfun(@(x)sscanf(x, '[%f, %f]').', ci, 'uni', false);
    ci = vertcat(ci{:});
else
    
    if opts.binary
        check_or = true(height(restab), 1);
    else % assumes beta and not exp(beta)
        check_or = false(height(restab), 1);
    end
    
    ci_tmp = restab;
    ci_tmp = renamevars(ci_tmp, [pcol, betaCol, opts.secol], ["p", "b", "se"]);
    if ciflag
        ci_tmp = renamevars(ci_tmp, CIcol, "ci");
    end

    ci_tmp = getCI(ci_tmp, ciflag, check_or);
    ci = ci_tmp{:, ["ci_l", "ci_h"]};
end

rmIdx = any(isinf(ci), 2) | any(isnan(ci), 2);
ci(rmIdx, :) = [];
restab(rmIdx, :) = [];

if isempty(restab)
    disp("no variant left in the table!")
    ax2 = [];
    return
end

% convert estimates if necessary (OR --> logOR) only if 'legacy' is true.
% Otherwise keep the original underlying format
if opts.legacy
    if any(strcmpi('or', restab.Properties.VariableNames)) || ...
        ~isnan(restab.AF_control(1)) % for calls from gwasrunner
        try
            restab.beta = log(restab.OR);
        catch
            restab.beta = log(restab.(betaCol));
        end
        ci = log(ci);
    else
        restab.beta = restab.(betaCol);
    end
else
    restab = renamevars(restab, betaCol, "beta");
end

% adjust P-values
padjusted = padjust(restab.(pcol), method=opts.padjust);
restab.hidep(:) = false; 
restab.hidep(padjusted > opts.pcutoff) = true;

if all(restab.hidep)
    disp("no variant remained significant after the multiple testing!")
    ax2 = [];
    return
end

% fetch region boundaries from gnomAD
if opts.verbose 
    fprintf('fetching required data from gnomAD...')
end
opts.neighbour = 1e3;
tryME = 1;
while tryME
    try
        r = gnomad('start', min(restab.POS)-opts.neighbour, ...
            'stop', max(restab.POS)+opts.neighbour, 'chrom', string(restab.CHR(1)),...
            'compact', true, 'dataset',char(opts.dataset),'reference_genome',char(opts.refGenome));
        geneList = {r.data.region.genes.gene_id};
        tryME = 0;
    catch
        if tryME >= 5
            error('rarePlotter:gnomAD server error:\n %s', r)
        else
            tryME = tryME + 1;
        end
    end
end
r = ({});
for i = 1:numel(geneList) % loop over genes in the region and fetch their transcripts info
    maxcount = 10; cnt = 1;
    while cnt <= maxcount
        rtmp = gnomad('gene_id', geneList{i}, 'compact', true,...
            'dataset',char(opts.dataset),'reference_genome',char(opts.refGenome));
        if isstruct(rtmp)
            r{i} = rtmp.data.gene; cnt = 11;
        else
            if cnt == maxcount % still server doesn't allow :(
                error('%s', rtmp)
            end
            pause(5) % too many genes to catch?
        end
        cnt = cnt + 1;
    end
    clear rtmp
end
r = cell2mat(r);

if opts.verbose; fprintf('\b\b Done.\n'); end

% check if genes have CDS/utr features
remove_idx = false(numel(r), 1);
for i = 1:numel(r)
    canon_idx = strcmp({r(i).transcripts.transcript_id},...
        r(i).canonical_transcript_id);
    checkR = r(i).transcripts(canon_idx).exons;
    if ~any(ismember({checkR.feature_type}, {'CDS', 'UTR'}))
        remove_idx(i) = true;
    end
end
r(remove_idx) = [];
[~, f_sorted] = sort([r.start]);
r = r(f_sorted);

% check gene overlap
boundary = 2e3;
[r.row] = deal(1);

for i = 1:numel(r)
    for j = i+1:numel(r)
        if (r(j).start >= r(i).start && r(j).start <= r(i).stop) || ...
                (r(j).stop >= r(i).start && r(j).stop <= r(i).stop)
            if r(i).row == r(j).row
                r(j).row = r(i).row - 1; % overlapping genes
            end
            
        elseif r(j).stop+boundary <= r(i).start || ...
                r(j).start-boundary <= r(i).stop
            if r(i).row == r(j).row
                r(j).row = r(i).row - 1; % adjacent genes
            end
        
        end
    end
end

% draw now
figure(Position=[1, 41, 1280, 720], Units="pixels")
hold on
textOffset = 0.45;
boxOffset = 0.25;

for j = 1:numel(r)
    plot([r(j).start r(j).stop], [r(j).row, r(j).row], ...
        'Color', '#252a34', 'linewidth', 4);

    canonicalIdx = ismember({r(j).transcripts.transcript_id}, ...
        r(j).canonical_transcript_id);
    rIN = r(j).transcripts(canonicalIdx).exons;

    if strcmp(r(j).strand, '+')
        txtarrow = r(j).symbol + ">"; % "\rightarrow";
    else
        txtarrow = "<" + r(j).symbol; % "\leftarrow"
    end

    text(r(j).start+(r(j).stop-r(j).start)/2, r(j).row+textOffset,...
        txtarrow, 'FontAngle', 'italic',...
            'HorizontalAlignment', 'center',...
            'FontSize', 8, 'Color', '#F04520', 'FontWeight', 'bold');

    for i = 1:numel(rIN)
        X = [rIN(i).start, rIN(i).stop, rIN(i).stop, rIN(i).start];
        Y = [r(j).row-boxOffset, r(j).row-boxOffset,...
            r(j).row+boxOffset, r(j).row+boxOffset];

        if strcmp(rIN(i).feature_type, 'CDS')
            patch('XData', X, 'YData', Y, 'FaceColor', '#252a34', ...
                'FaceAlpha', 1,...
                'EdgeColor', '#252a34');
        elseif strcmp(rIN(i).feature_type, 'UTR')
            patch('XData', X, 'YData', Y, 'FaceColor', '#ff2e63', ...
                'FaceAlpha', 1,...
                'EdgeColor', '#ff2e63');
        end
    end
end

ax1 = gca;
ax1.YLim(2) = max([r.row]) + textOffset + 0.1;
ax1.YLim(1) = ax1.YLim(1) - 0.05;
ax1.Position(1) = 0.05;
ax1.Position(2) = 0.05; ax1.Position(3) = 0.9;
ax1.Position(4) = ax1.Position(4)/8;
ax1.XAxis.Visible = 'off';
ax1.YAxis.Visible = 'off';
newpos = ax1.Position;
newpos(2) = newpos(2) + newpos(4) + 0.075;
newpos(4) = ax1.Position(4)*7;
ax2 = axes('Position', newpos);
ax2.Box = 'off';

hold on

grid(ax2, 'on')

if opts.binary
    restab.isprotected = restab.beta < 1;
else
    restab.isprotected = restab.beta < 0;
end

if opts.hidealpha > 5e-2 % otherwise almost white
    opts.hide = false;
else
    opts.hide = true;
end

% hide non-significant associations
if ~opts.hide % otherwise almost white
    hide_color = (1-opts.hidealpha).*ones(1, 3);
    line(ax2, [restab.POS(restab.hidep, :), restab.POS(restab.hidep, :)].', ...
        ci(restab.hidep, :).', Color=hide_color, MarkerFaceColor=hide_color, ...
        MarkerEdgeColor=hide_color);
end

idx_increase = ~restab.isprotected & ~restab.hidep;
line(ax2, [restab.POS(idx_increase, :), restab.POS(idx_increase, :)].', ci(idx_increase, :).', 'color', 'b');
idx_protect = restab.isprotected & ~restab.hidep;
line(ax2, [restab.POS(idx_protect, :), restab.POS(idx_protect, :)].', ci(idx_protect, :).', 'color', 'r');

% % limix yaxis to min/max ci
% mmYlim = [min(ci(~restab.hidep, :), [], "all"), max(ci(~restab.hidep, :), [], "all")];
% mmYlim_offset = abs(0.05.*mmYlim);
% ax2.YLim = [mmYlim(1) - mmYlim_offset(1), mmYlim(2) + mmYlim_offset(2)];

if isfield(opts, "consequence") && any(colnames(restab) == opts.consequence)% set markers based on consequences
    restab = renamevars(restab, opts.consequence, "consequence");

    % definition of these categories are based on gnomAD
    pLoF = ["transcript_ablation", "splice_acceptor_variant", ...
        "splice_donor_variant", "stop_gained", "frameshift_variant"];
    Missense_inframeINDEL = ["stop_lost", "start_lost", ...
        "inframe_insertion", "inframe_deletion", "missense_variant"];
    Synonymous = "synonymous_variant";
    restab.consCat = ones(height(restab), 1);
    restab.consCat(contains(restab.consequence, Synonymous)) = 2;
    restab.consCat(contains(restab.consequence, Missense_inframeINDEL)) = 3;
    restab.consCat(contains(restab.consequence, pLoF)) = 4;
    
    consCat = unique(restab.consCat);
    markerMeaning = ["Other", "Synonymous", "Missense / Inframe indel", "pLoF"];
    markerList = {'v', 'square', 'O', 'diamond'};
    markerListb = {'downtriangle', 'square', 'circle', 'diamond'};

    lg = (0);
    for i = 1:numel(consCat)
        subTab = restab(restab.consCat == consCat(i), :); 

        idx_pro_hide = subTab.isprotected & subTab.hidep;
        idx_pro = subTab.isprotected & ~subTab.hidep;
        idx_inc_hide = ~subTab.isprotected & subTab.hidep;
        idx_inc = ~subTab.isprotected & ~subTab.hidep;

        % hide non-significant effects
        if opts.method == "scatter"

            if ~opts.hide % otherwise almost white
                scatter(ax2, subTab.POS(idx_pro_hide), ...
                    subTab.beta(idx_pro_hide),...
                    opts.sz, -log10(subTab.(pcol)(idx_pro_hide)), 'filled',...
                    'MarkerEdgeColor', 'r', 'Marker', markerList{consCat(i)}, ...
                    'MarkerEdgeAlpha', opts.hidealpha, 'MarkerFaceAlpha', opts.hidealpha);
        
                scatter(ax2, subTab.POS(idx_inc_hide), ...
                    subTab.beta(idx_inc_hide),...
                    opts.sz, -log10(subTab.(pcol)(idx_inc_hide)), 'filled',...
                    'MarkerEdgeColor', 'b', 'Marker', markerList{consCat(i)}, ...
                    'MarkerEdgeAlpha', opts.hidealpha, 'MarkerFaceAlpha', opts.hidealpha);
            end
    
            scatter(ax2, subTab.POS(idx_pro), ...
                subTab.beta(idx_pro),...
                opts.sz, -log10(subTab.(pcol)(idx_pro)), 'filled',...
                'MarkerEdgeColor', 'r', 'Marker', markerList{consCat(i)});
    
            scatter(ax2, subTab.POS(idx_inc), ...
                subTab.beta(idx_inc),...
                opts.sz, -log10(subTab.(pcol)(idx_inc)), 'filled',...
                'MarkerEdgeColor', 'b', 'Marker', markerList{consCat(i)});
        else % bubblechart
           try
               if ~opts.hide % otherwise almost white
                   bax = bubblechart(ax2, subTab.POS(idx_pro_hide), ...
                        subTab.beta(idx_pro_hide), subTab.(opts.mafCol)(idx_pro_hide), ...
                        -log10(subTab.(pcol)(idx_pro_hide)), ...
                        "MarkerEdgeColor", 'r', "MarkerFaceAlpha", opts.hidealpha, ...
                        'MarkerEdgeAlpha', opts.hidealpha);
                    bax.MarkerHandle.Style = markerListb{consCat(i)};
               end
           catch % 1 empty vector
           end
           
           try
               if ~opts.hide % otherwise almost white
                   bax = bubblechart(ax2, subTab.POS(idx_inc_hide), ...
                        subTab.beta(idx_inc_hide), subTab.(opts.mafCol)(idx_inc_hide), ...
                        -log10(subTab.(pcol)(idx_inc_hide)), ...
                        "MarkerEdgeColor", 'b', "MarkerFaceAlpha", opts.hidealpha, ...
                        'MarkerEdgeAlpha', opts.hidealpha);
                    bax.MarkerHandle.Style = markerListb{consCat(i)};
               end
           catch
           end
            
           try
                bax = bubblechart(ax2, subTab.POS(idx_pro), ...
                    subTab.beta(idx_pro), subTab.(opts.mafCol)(idx_pro), ...
                    -log10(subTab.(pcol)(idx_pro)), ...
                    "MarkerEdgeColor", 'r', "MarkerFaceAlpha", 1);
                bax.MarkerHandle.Style = markerListb{consCat(i)};
           catch
           end
            
           try
                bax = bubblechart(ax2, subTab.POS(idx_inc), ...
                    subTab.beta(idx_inc), subTab.(opts.mafCol)(idx_inc), ...
                    -log10(subTab.(pcol)(idx_inc)), ...
                    "MarkerEdgeColor", 'b', "MarkerFaceAlpha", 1);
                bax.MarkerHandle.Style = markerListb{consCat(i)};
           catch
           end

        end

        
        % for legend
        lg(i) = plot(ax2, nan, nan, 'MarkerFaceColor', '#fce38a',...
            'MarkerEdgeColor', 'k', 'LineStyle', 'none',...
            'Marker', markerList{consCat(i)});
    end
    
    leg = legend(lg, 'AutoUpdate', 'off',...
        'String', markerMeaning(consCat), 'Location', 'northeastoutside',...
        'Color', '#BAD8E0');
    leg.Title.String = 'VEP annotation';
    leg.TextColor = '#3F3F3F';
 
else
    
    idx_pro_hide = restab.isprotected & restab.hidep;
    idx_pro = restab.isprotected & ~restab.hidep;
    idx_inc_hide = ~restab.isprotected & restab.hidep;
    idx_inc = ~restab.isprotected & ~restab.hidep;

    
    if opts.method == "scatter"
        % hide non-significant effects
        if ~opts.hide % otherwise almost white
            scatter(ax2, restab.POS(idx_pro_hide), ...
                restab.beta(idx_pro_hide), opts.sz,...
                -log10(restab.(pcol)(idx_pro_hide)), ...
                'filled', 'MarkerEdgeColor', 'r', ...
                'MarkerEdgeAlpha', opts.hidealpha, 'MarkerFaceAlpha', opts.hidealpha);
        
            scatter(ax2, restab.POS(idx_inc_hide), ...
                restab.beta(idx_inc_hide), opts.sz,...
                -log10(restab.(pcol)(idx_inc_hide)), ...
                'filled', 'MarkerEdgeColor', 'b', ...
                'MarkerEdgeAlpha', opts.hidealpha, 'MarkerFaceAlpha', opts.hidealpha);
        end
    
        scatter(ax2, restab.POS(idx_pro), ...
            restab.beta(idx_pro), opts.sz,...
            -log10(restab.(pcol)(idx_pro)), ...
            'filled', 'MarkerEdgeColor', 'r');
    
        scatter(ax2, restab.POS(idx_inc), ...
            restab.beta(idx_inc), opts.sz,...
            -log10(restab.(pcol)(idx_inc)), ...
            'filled', 'MarkerEdgeColor', 'b');
    else % bubblechart
        
        if ~opts.hide % otherwise almost white
            bubblechart(ax2, restab.POS(idx_pro_hide), ...
                restab.beta(idx_pro_hide), restab.(opts.mafCol)(idx_pro_hide), ...
                -log10(restab.(pcol)(idx_pro_hide)), ...
                "MarkerEdgeColor", 'r', "MarkerFaceAlpha", opts.hidealpha, ...
                'MarkerEdgeAlpha', opts.hidealpha);
    
            bubblechart(ax2, restab.POS(idx_inc_hide), ...
                restab.beta(idx_inc_hide), restab.(opts.mafCol)(idx_inc_hide), ...
                -log10(restab.(pcol)(idx_inc_hide)), ...
                "MarkerEdgeColor", 'b', "MarkerFaceAlpha", opts.hidealpha, ...
                'MarkerEdgeAlpha', opts.hidealpha);
        end

        bubblechart(ax2, restab.POS(idx_pro), ...
            restab.beta(idx_pro), restab.(opts.mafCol)(idx_pro), ...
            -log10(restab.(pcol)(idx_pro)), ...
            "MarkerEdgeColor", 'r', "MarkerFaceAlpha", 1);

        bubblechart(ax2, restab.POS(idx_inc), ...
            restab.beta(idx_inc), restab.(opts.mafCol)(idx_inc), ...
            -log10(restab.(pcol)(idx_inc)), ...
            "MarkerEdgeColor", 'b', "MarkerFaceAlpha", 1);

    end
    
end

c = colorbar(ax2, 'Location', 'north');
c.Position(1) = ax2.Position(1);
c.Position(2) = ax2.Position(2) + ax2.Position(4);
c.Position(4) = c.Position(4)/2;
c.Position(3) = ax2.Position(3)/4;
c.AxisLocation  = 'out';
c.Label.String = "-log_{10} (\itp)";
c.Label.Units = 'normalized';
c.Label.HorizontalAlignment = 'left';
c.Label.VerticalAlignment = 'middle';
c.Label.FontSize = 11;
c.Label.Position = [1 0.5 0];

set(ax2,'box','off','LineWidth',1.1, 'TickLength', [0.001 0], ...
    'FontSize', 11,...
    'FontWeight', 'bold', 'TickDir', 'out')

try
    opts.color = feval(opts.colormap);
catch
    opts.color = get_colormap(opts.colormap, 256);
end
colormap(ax2, opts.color);

ylabel(opts.ylabel, 'FontSize', 14, 'FontWeight', 'bold');
drawnow;
ax2.XLim = ax1.XLim;
ax2.XTickLabel = string(ax2.XTick./1e6);
ax2.XAxis.FontSize = 10;
xlabel('chromosome ' + string(restab.CHR(1)) + ' (Mb)', 'FontSize', 11, ...
    'FontWeight', 'bold');
ax2.GridAlpha = 0.025;

% check if new ax1 limits covers gene names
tobj = findobj(ax1,'Type','text');
for i = 1:numel(tobj)
    if tobj(i).Position(1) <= ax1.XLim(1)
        if isscalar(tobj)
            tobj(i).Position(1) = mean(ax1.XLim);
            tobj(i).HorizontalAlignment = 'center';
        else
            tobj(i).Position(1) = ax1.XLim(1);
            tobj(i).HorizontalAlignment = 'left';
        end
    elseif tobj(i).Position(1) >= ax1.XLim(2)
        if isscalar(tobj)
            tobj(i).Position(1) = mean(ax1.XLim);
            tobj(i).HorizontalAlignment = 'center';
        else
            tobj(i).Position(1) = ax1.XLim(2);
            tobj(i).HorizontalAlignment = 'right';
        end
    end
end

if isfield(opts, "title")
    title(ax2, opts.title, 'HorizontalAlignment', 'left');
elseif any(colnames(restab) == "Pheno")
    title(ax2, restab.Pheno(1), 'HorizontalAlignment', 'left');
end

% shade exonic regions (CDS/UTR)
if opts.shade
    ax2 = shade(ax2, 'step', sort([rIN.start, rIN.stop]), 'dir', 'x', ...
        'FaceAlpha', 0.1);
end

% set bubblesize and add bubble legend
if opts.method == "bubble"
    bubblesize(ax2, opts.bsize);
    bax = bubblelegend(ax2, NumBubbles=2, Location="northeastoutside");
    bax.Color = '#BAD8E0';
    bax.Title.String = "MAF";
    bax.TextColor = '#3F3F3F';
    bax.LineWidth = 1.1;
    drawnow;
    limlabs = bax.LimitLabels;
    bax.LimitLabels = compose("%.2g", double(string(limlabs)));

    % align below current legend
    leg = findobj(ax2.Parent, "Type", "Legend");
    if ~isempty(leg)
        drawnow;
        bax.Position(1) = leg.Position(1);
        if bax.Position(1) ~= leg.Position(1)
            drawnow; bax.Position(1) = leg.Position(1);
        end
        bax.Position(2) = leg.Position(2) - bax.Position(4) - 0.005;

    end

    if (sum(~restab.hidep) < 2) % don't show the legend
        bax.Visible = "off";
    end
end

% realign ax1 and ax2
drawnow;
ax1.Position(3) = ax2.Position(3);

opts.output = string(opts.output);
if opts.output == "" || ismissing(opts.output)
    try % if 'Pheno' col is present
        opts.output = string(matlab.lang.makeValidName(restab.Pheno(1)))+".rareForest";
    catch
        opts.output = "rareForest";
    end
end

if opts.save
    exportgraphics(ax2.Parent, opts.output + "." + opts.format, 'Resolution', opts.resolution)
end

if opts.savefig
    savefig(ax2.Parent, opts.output + ".fig")
end

end