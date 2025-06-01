function genePlotter(sumstat, p)
% plots regional/gene plots and is a subfunction of manplotter.

% hint:
% Exons = gene - introns
% CDS = gene - introns - UTRs 
% CDS = Exons - UTRs

% merged with manplotter---------------------------------------------------
arguments
    sumstat {mustBeA(sumstat, 'table')}
    p.finemap {mustBeA(p.finemap, 'table')} % fine-mapping summary stat (must have cs and pip columns, case insensitive)
    p.gene {mustBeText} = "-"
    p.region {mustBeVector, mustBeNumeric} = [nan, nan]
    p.flag {mustBeTextScalar, mustBeMember(p.flag, ["region", "gene"])} = "region" % default
    % p.annotate {mustBeTextScalar} = "SNP" % column to be treated as "annotate" column. Used only if "annotate" is absent from variables.
    
    % for  genome-wide, 5e-8 will be used; for fdr, BH adjusted p-value
    % will be used. If set to none, no line will be plotted.
    p.significance {mustBeTextScalar, mustBeMember(p.significance, ["genome-wide", "bonferroni", "fdr", "none"])} = "genome-wide"
    p.lead {mustBeTextScalar} % lead variant (diamond symbol) in the region. If left empty, the strongest association (largest -log10P) is selected as lead variant.
    p.extralead {mustBeText, mustBeVector} = "" % additional markers to be shown on the plot (will be shown as square). 
    p.causalLead (1,1) logical = false % if lead variant (diamond symbol) is also a cusal variant (large PIP), add an asterisk (*) to the lead variant's label.
    p.labelExtralead (1,1) logical = false % to label 'extralead' variants
    p.neighbour (1,1) double {mustBeNonnegative} = 1 % margin around the region, in bp, default is 1 (modification not recommended).
    p.title {mustBeTextScalar} = "" % title to be diplayed on the graph
    p.backup (1,1) logical = false % useful if called from gwascaller.
    p.format {mustBeMember(p.format, ["jpg", "png", "tif", "eps", "pdf"])} = "jpg" % format of output plots
    p.resolution (1, 1) double = 400 % resolution of output plot
    p.save (1,1) logical = true % only if plot is true
    p.savefig (1,1) logical = false % only if plot is true. Saves .fig files as well
    p.outdir {mustBeTextScalar} % output directory for backup files and figures
    p.ldmethod {mustBeTextScalar, mustBeMember(p.ldmethod, ["insample", "external"])} = "insample" % uses bgen file specified in bgenhome
    p.musthaveUTR (1,1) logical = true % to visualize only genes which have UTR (and CDS)

    % options for "insample" 'ldmethod'
    p.bgenhome {mustBeFolder} % if empty, uses default directory in getbulkgeno function
    p.parallel (1,1) logical = false % uses parallel toolbox for in-sample LD structure
    
    % aesthetic options
    p.genecolor (1, 2) {mustBeText, mustBeMember(p.genecolor, ["Coral", "OrangeRed",...
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
        "DarkGray", "Gray", "SlateGray", "DarkSlateGray", "Black"])} = ["DarkRed", "Red"]; % CDS/UTR colors
    p.boundary (1, 1) double = 3e4 % the boundary around each gene for plotting, if genes' distance is > boundary, they'll be shifted up/down-wards
    % p.textOffset (1, 1) double = 0.25 % DEPRECATED, now textOffset = boxOffset + 0.05
    p.boxOffset (1, 1) double = 0.25
    p.linewidth (1, 1) double {mustBeGreaterThan(p.linewidth, 1e-3)} = 1.2 % line width of genes to draw
    p.xfontsize (1,1) double = 12 % font size of x-axis, x-label font size is this value + 2
    p.yfontsize (1,1) double = 14 % font size of y-axis, y-label font size is this value + 3
    p.box (1,1) logical = false % box around the plot
    p.markersize (1, 1) double {mustBeGreaterThan(p.markersize, 1)} = 50 % marker size of dots
    p.genefontsize (1,1) double = 11 % gene name font size, legend font size is this value - 2
    p.topfontsize (1,1) double = 14 % font size of the top signal 
    p.map (1,1) {mustBeMember(p.map, ["jet", "turbo", "abyss", ...
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
        "flare", "crest", "colorblind", "deep", "muted"])} = "turbo" % colormap for LD structure
    p.MarkerFaceAlpha (1, 1) double {mustBeInRange(p.MarkerFaceAlpha, 0, 1)} = 0.95
    p.MarkerEdgeAlpha (1, 1) double {mustBeInRange(p.MarkerEdgeAlpha, 0, 1)} = 0.9
    p.aspectratio (1, 1) double {mustBePositive} = 4 % aspect ratio of gene to regional plot
    p.aspectratioMethod {mustBeTextScalar, mustBeMember(p.aspectratioMethod, ["auto", "manual"])} = "auto" % if "auto" aspectratio will be set based on gene rows
    p.axdistance (1, 1) double {mustBePositive} = 0.085 % margin between regional and gene plot
    p.fontname {mustBeTextScalar} = "Garamond"
    p.randSample (1,1) double = 2e5 % number of random samples to be used in LD calculation 

end
% -------------------------------------------------------------------------

% validate input arguments
close all force

% columns to be always present: chr, pos, p
sumstat.Properties.VariableNames = lower(sumstat.Properties.VariableNames);
idx = ismember({'chr', 'pos', 'p'}, sumstat.Properties.VariableNames);
if ~all(idx)
    error('genePlotter:missingArgument', 'chr, pos and p must be present!')
end

% check snp column
snpcol = find(ismember(lower(sumstat.Properties.VariableNames), ...
    ["snp", "variant_id", "variantid", "variant", "id", "annotate"]), 1);
if isempty(snpcol)
    error('genePlotter:missingArgument', 'snp column is missing')
end
sumstat.Properties.VariableNames(snpcol) = "snp";

% check if fine-mapping stats are present
if isfield(p, 'finemap')
    p.finemap.Properties.VariableNames = lower(p.finemap.Properties.VariableNames);
    fmidx = startsWith({'cs', 'pip'}, p.finemap.Properties.VariableNames);
    if ~all(fmidx)
        error('genePlotter:missingArgument', "cs and pip pattern in " + ...
            "column names must be present for fine-mapping!")
    end
    
    % snp col
    [idx1, idx2] = ismember(sumstat.snp, p.finemap.snp); idx2(idx2 < 1) = [];
    [sumstat.cs, sumstat.pip] = deal(zeros(height(sumstat), 1));
    sumstat.cs(idx1) = p.finemap.cs(idx2);
    sumstat.pip(idx1) = p.finemap.pip(idx2);
    p.finemap = true;
else
    p.finemap = false;
end

if numel(unique(sumstat.chr)) > 1
    error('genePlotter:chr', 'chr in sumstat table must be unique!')
end
p.chr = char(string(sumstat.chr(1)));

if ~isstring(sumstat.snp)
    sumstat.snp = string(sumstat.snp);
end

if strcmp(p.flag, "region")
    if all(isnan(p.region)) % region left empty, find it from pos column
        p.region = [min(sumstat.pos), max(sumstat.pos)];
    end
elseif isempty(p.gene) || p.gene == "" || ismissing(p.gene)
    error('genePlotter:missingGene', 'gene name cannot be empty with "gene" flag!')
end

% check p-values
if any(sumstat.p > 1) || any(sumstat.p < 0)
    error('genePlotter:pvalue', 'p values must be in [0 1] range')
end

if any(sumstat.p < realmin*eps)
    sumstat.p(sumstat.p < realmin*eps) = realmin*eps;
end
sumstat.p = -log10(sumstat.p);
% -------------------------------------------------------------------------
if ~isfield(p, 'outdir') || string(p.outdir) == "" || ismissing(p.outdir)
    p.outdir = pwd;
end
if ~isfolder(p.outdir), mkdir(p.outdir); end

if p.backup % generate a name for backup file
    p.backupname = "genePlotter.chr" + string(p.chr) + ".pos" + ...
        min(sumstat.pos) + "_" + max(sumstat.pos) + ".mat";
    p.backupname = fullfile(p.outdir, p.backupname);
end

% prepare LD matrix 
fprintf('getting LD matrix...')

if p.backup && isfile(p.backupname) % save for other traits (useful in phenom-wide analysis)
    ldmat = load(p.backupname, 'ldmatbackup').ldmatbackup; 
else
    if p.ldmethod == "external" % 1000 Genomes project (European subset)
        ldmat = LDlink(sumstat.snp);
    else % in-sample, uses subset of white-British from UK Biobank (as defined in UK Biobank Nature paper)
        if any(colnames(sumstat).lower == "allele0")
            sumstat = renamevars(sumstat, ["allele0", "allele1"], ["a1", "a2"]);
        end
        eids = getQCEID(3, false);
        p.randSample = min(p.randSample, numel(eids));
        eids = randsample(eids, p.randSample);
        if isfield(p, 'bgenhome')
            %@29OCT2024
            if contains(p.bgenhome.lower, "topmed")
                maf = 1e-5;
            else
                maf = 1e-3;
            end

            ldmat = getInsampleLD(sumstat.snp, string(sumstat.chr), ...
                parallel=p.parallel, bgenhome=p.bgenhome, ...
                matchid=sumstat.chr+":"+sumstat.pos+":"+sumstat.a1+":"+sumstat.a2, ...
                removenan=false, ldsample=eids, maf=maf);
        else
            ldmat = getInsampleLD(sumstat.snp, string(sumstat.chr), ...
                parallel=p.parallel, ...
                matchid=sumstat.chr+":"+sumstat.pos+":"+sumstat.a1+":"+sumstat.a2, ...
                removenan=false, ldsample=eids, maf=1e-3);
        end
        ldmat.ld = ldmat.ld.^2;
    end

    if p.backup && ~isfile(p.backupname)
        ldmatbackup = ldmat;
    end
end

if isfield(p, 'lead')
    tophit = find(sumstat.snp == p.lead, 1);
else
    [~, tophit] = max(sumstat.p);
end
tophit_ld = ldmat.snp == sumstat.snp(tophit);

p.extralead(p.extralead == "" | ismissing(p.extralead)) = [];
if ~isempty(p.extralead)
    sumstat.extralead = ismember(sumstat.snp, p.extralead);
    sumstat.extralead(tophit) = false;
else
    sumstat.extralead = false(height(sumstat), 1);
end

% get LD structure 
Cdata = -1.*ones(size(sumstat,1), 1);
if any(tophit_ld)
    [f1, f2] = ismember(sumstat.snp, ldmat.snp); f2(f2<1) = [];
    ldmat.ld(logical(eye(size(ldmat.ld, 1)))) = 0; % pairwise r2 to itself
    Cdata(f1) = ldmat.ld(tophit_ld, f2);
    Cdata(isnan(Cdata)) = -1;
else
    ldmat.ld = NaN;
end
fprintf('\b\b Done.\n')


fprintf('getting gene/region data from gnomAD...')

%@29OCT2024: topmed is on GRCh38
if isfield(p, "bgenhome") && contains(p.bgenhome.lower, "topmed")
    dataset = 'gnomad_r3';
    reference_genome = 'GRCh38';
else
    dataset = 'gnomad_r2_1';
    reference_genome = 'GRCh37';
end

if p.backup && isfile(p.backupname)
    r = load(p.backupname, 'r').r; 
else
    if strcmp(p.flag, 'gene')

        r = gnomad('gene_symbol', p.gene, 'dataset', dataset, 'reference_genome', reference_genome); % default is GRCh37
        r = r.data.gene;

        % get genes in the region
        r = gnomad('start', r.start-p.neighbour, ...
            'stop', r.stop+p.neighbour, 'chrom', r.chrom, ...
            'dataset', dataset, 'reference_genome', reference_genome);

    else % region is provided
        maxcount = 20; cnt = 1;
        while cnt <= maxcount
            r = gnomad('start', p.region(1)-p.neighbour, ...
                'stop', p.region(2)+p.neighbour, 'chrom', p.chr, ...
                'compact', true, 'dataset', dataset, ...
                'reference_genome', reference_genome);
            if isstruct(r)
                geneList = {r.data.region.genes.gene_id};
                break
            else
                if cnt == maxcount % still server doesn't allow :(
                    error('%s', r)
                end
                pause(3) % too many genes to catch?
            end
            cnt = cnt + 1;
        end

        r = ({});
        for i = 1:numel(geneList) % loop over genes in the region and fetch their transcripts info
            maxcount = 20; cnt = 1;
            while cnt <= maxcount
                rtmp = gnomad('gene_id', string(geneList{i}), ...
                    'compact', true, 'dataset', dataset, ...
                    'reference_genome', reference_genome);
                if isstruct(rtmp)
                    r{i} = rtmp.data.gene;
                    break
                else
                    if cnt == maxcount % still server doesn't allow :(
                        error('%s', rtmp)
                    end
                    pause(7) % too many genes to catch?
                end
                cnt = cnt + 1;
            end
            clear rtmp
        end
        r = cell2mat(r);
    end
    
    % check duplicated genes (same symbol but different gene id). In this case,
    % only genes returned by gnomad API based on symbol search are picked.
    dup_genes = duplicates({r.symbol});
    rmIdx = ({});
    for i = 1:numel(dup_genes)
        idx = ismember({r.symbol}, dup_genes{i});
        maxcount = 20; cnt = 1;
        while cnt <= maxcount
            try
                gene_tmp = gnomad('gene_symbol', dup_genes{i}, 'dataset', dataset, 'reference_genome', reference_genome);
                gene_tmp = {gene_tmp.data.gene.gene_id};
                break
            catch ME
                if cnt == maxcount % still server doesn't allow :(
                    error('%s', ME.message)
                end
                pause(3) % too many genes to catch?
            end
            cnt = cnt + 1;
        end
        
        rmIdx{i} = find(idx & ~ismember({r.gene_id}, gene_tmp));
    end
    r(horzcat(rmIdx{:})) = [];
end
fprintf(' Done.\n')

if p.backup && ~isfile(p.backupname)% save for other traits (useful in phenom-wide analysis)
    save(p.backupname, 'r', 'ldmatbackup'); 
end

% check if genes have CDS/utr features
if p.musthaveUTR
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
end
[~, f_sorted] = sort([r.start]);
r = r(f_sorted);

% check gene overlap ------------------------------------------------------
if numel(r) > 40
    p.boundary = max(p.boundary, 5e4);
elseif numel(r) > 20
    p.boundary = max(p.boundary, 4e4);
end

[r.row] = deal(1);
ovr = false(numel(r), numel(r)); % overlapping matrix
for i = 1:numel(r)
    for j = i+1:numel(r)
%         if (r(j).start >= r(i).start && r(j).start <= r(i).stop) || ...
%                 (r(j).stop >= r(i).start && r(j).stop <= r(i).stop)
%             if r(i).row == r(j).row
%                 r(j).row = r(i).row - 1; % overlapping genes
%             end
%             [ovr(i, j), ovr(j, i)] = deal(true);
            
        if r(j).stop+p.boundary <= r(i).start || ...
                r(j).start-p.boundary <= r(i).stop
            if r(i).row == r(j).row
                r(j).row = r(i).row - 1; % adjacent genes
            end
           [ovr(i, j), ovr(j, i)] = deal(true);
        
        end
    end
end

% rearrange genes so they have furthest distance to each
% % other
% ovr(logical(eye(size(ovr, 1)))) = true;
% meanpos = mean([r.start; r.stop], 1)';
% meandist = squareform(pdist(meanpos, 'cityblock'));
% y = [r.row]';
% mima = minmax(y');
% ct = 1; ytmp = ({});
% while numel(ovr)
%     ovr_idx = find(ovr(1, :));
%     meandist_tmp = meandist(ovr_idx, ovr_idx);
% 
%     ovr(ovr_idx, :) = []; ovr(:, ovr_idx) = [];
%     meandist(ovr_idx, :) = []; meandist(:, ovr_idx) = [];
% 
%     checky = round(linspace(mima(1), mima(2), numel(ovr_idx)));
%     [~, minidx] = sort(meandist_tmp(1, :));
%     ytmp{ct} = checky(minidx);
%     ct = ct + 1;
% end
% ytmp = horzcat(ytmp{:});
% for i = 1:numel(r); r(i).row = ytmp(i); end

% check if aspect ratio of gene/hit axes is auto
p.legfontsize = p.genefontsize - 2;
if p.aspectratioMethod == "auto"
    numrows = range([r.row]) + 1;
    if numrows < 2
        p.aspectratio = 12;
    elseif numrows < 3
        p.aspectratio = 8.2;
    elseif numrows < 5
        p.aspectratio = 6;
    elseif numrows < 6
        % p.genefontsize = 9;
        p.aspectratio = 4;
    elseif numrows < 8
        % p.genefontsize = 8.5;
        p.aspectratio = 3;
    elseif numrows < 10
        % p.genefontsize = 7;
        p.aspectratio = 2;
    else
        % p.genefontsize = 5.5;
        p.aspectratio = 1.5;
    end
end

%% draw now ---------------------------------------------------------------
figure("WindowState", "maximized"); 
hold on

% resolve gene colors (CDS/UTR)
pc = palette();
[~, idx] = ismember(p.genecolor, pc.color);
p.genecolor = pc.code(idx);

p.textOffset = p.boxOffset + 0.05;
for j = 1:numel(r) % draw genes 
    plot([r(j).start r(j).stop], [r(j).row, r(j).row], ...
        'Color', p.genecolor(1), 'linewidth', p.linewidth);

    canonicalIdx = ismember({r(j).transcripts.transcript_id}, ...
        r(j).canonical_transcript_id);
    rIN = r(j).transcripts(canonicalIdx).exons;

    if strcmp(r(j).strand, '+')
        txtarrow = r(j).symbol + ">"; % "\rightarrow";
    else
        txtarrow = "<" + r(j).symbol; % "\leftarrow"
    end

    text(r(j).start+(r(j).stop-r(j).start)/2, r(j).row+p.textOffset,...
        txtarrow, 'FontAngle', 'italic', ...
        'HorizontalAlignment', 'center',...
        'FontSize', p.genefontsize, 'Tag', r(j).symbol, ...
        'VerticalAlignment', 'baseline', 'FontName', p.fontname);

    for i = 1:numel(rIN)
        X = [rIN(i).start, rIN(i).stop, rIN(i).stop, rIN(i).start];
        Y = [r(j).row-p.boxOffset, r(j).row-p.boxOffset,...
            r(j).row+p.boxOffset, r(j).row+p.boxOffset];

        if strcmp(rIN(i).feature_type, 'CDS')
            patch('XData', X, 'YData', Y, 'FaceColor', p.genecolor(1), ...
                'FaceAlpha', 1, 'Tag', r(j).symbol, ...
                'EdgeColor', p.genecolor(1));
        elseif strcmp(rIN(i).feature_type, 'UTR')
            patch('XData', X, 'YData', Y, 'FaceColor', p.genecolor(2), ...
                'FaceAlpha', 1, 'Tag', r(j).symbol, ...
                'EdgeColor', p.genecolor(2));
        end
    end
end

ax1 = gca;
ax1.YLim(1) = ax1.YLim(1) - 0.05;
ax1.Position(1) = 0.05;
ax1.Position(2) = 0.05; ax1.Position(3) = 0.9;
ax1.Position(4) = ax1.Position(4)/p.aspectratio;

% find best YLim for ax2
te = findobj(ax1, 'Type', 'Text');
drawnow
te_extent = reshape([te.Extent], 4, [])';
[~, max_ext] = max(te_extent(:, 2) + te_extent(:, 4)); % closest to upper y-axis
ax1.YLim(2) = te(max_ext).Extent(2) + te(max_ext).Extent(4);

ax1.XAxis.Visible = 'off';
ax1.YAxis.Visible = 'off';
newpos = ax1.Position;
newpos(2) = newpos(2) + newpos(4) + p.axdistance;
newpos(4) = ax1.Position(4)*(p.aspectratio - 1);
ax2 = axes('Position', newpos);
ax2.Box = 'off';

% make sure that ax2 remains within figure margin
drawnow
ax2.Position(4) = 0.95 - ax2.Position(2);

hold on
% add fine-mapping results (if available) ---------------------------------
if p.finemap
    csgrp = unique(sumstat.cs);
    csgrp(csgrp == 0) = [];
    if isempty(csgrp)
        p.finemap = false;
    else
        % fmap_color = feval(p.map, numel(csgrp));
        fmap_color = ["Black", "RoyalBlue", "Cyan", "Gold", "DarkSlateGray"];
        [~, fmapCidx] = ismember(fmap_color, pc.color);
        fmap_color = pc.code(fmapCidx);
        if numel(csgrp) > 5
            csgrp = csgrp(1:5);
        end
        
        sc3 = (0);
        for i = 1:numel(csgrp)
            fmhits = sumstat.cs == csgrp(i);
            fmhits(union(tophit, find(sumstat.extralead))) = false;
            sc3(i) = scatter(sumstat.pos(fmhits),...
                        sumstat.p(fmhits), 'MarkerFaceColor', ...
                        'none', 'MarkerEdgeColor', fmap_color(csgrp(i)), ...
                        'MarkerEdgeAlpha', 1,...
                        'MarkerFaceAlpha', 1, ...
                        'SizeData', p.markersize*2, ...
                        'LineWidth', p.linewidth);
    %         sc3(i) = scatter(sumstat.pos(fmhits),...
    %                     sumstat.p(fmhits), p.markersize*4, sumstat.pip(fmhits), ...
    %                     'MarkerFaceColor', ...
    %                     'flat', 'MarkerEdgeColor', 'k', ...
    %                     'MarkerEdgeAlpha', 1,...
    %                     'MarkerFaceAlpha', 1);
        end
        
        if sumstat.cs(tophit) > 0 
            tophitcs = sumstat.cs(tophit);
            sc3(numel(csgrp) + 1) = scatter(sumstat.pos(tophit), ...
                sumstat.p(tophit), 'MarkerFaceColor', 'none',...
                'MarkerEdgeColor', fmap_color(tophitcs), ...
                'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 1, ...
                'Marker', 'diamond',...
                'SizeData', p.markersize*4, ...
                'LineWidth', p.linewidth);
        end
        
        % check extra lead
        if any(sumstat.cs(sumstat.extralead) > 0)
            sumstat_extra = sumstat(setdiff(find(sumstat.extralead), tophit), :);
            for i = 1:numel(csgrp)
                fmhits = sumstat_extra.cs == csgrp(i);
                if ~any(fmhits); continue; end
                sc3(numel(sc3) + 1) = scatter(sumstat_extra.pos(fmhits), ...
                    sumstat_extra.p(fmhits), 'MarkerFaceColor', 'none',...
                    'MarkerEdgeColor', fmap_color(csgrp(i)), ...
                    'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 1, ...
                    'Marker', 'square',...
                    'SizeData', p.markersize*5, ...
                    'LineWidth', p.linewidth);
            end
        end

    end

end

% create LD bins ----------------------------------------------------------
try
    CdataCols = feval(p.map, 6);
catch
    CdataCols = get_colormap(p.map, 6);
end

lbls = flip(["0.8\leqr^{2}<1.0", "0.6\leqr^{2}<0.8", ...
    "0.4\leqr^{2}<0.6", "0.2\leqr^{2}<0.4", "0.0\leqr^{2}<0.2", ...
    "no r^{2} info"]);
CdataBins = zeros(numel(Cdata), 3);
for i = 1:3
    CdataBins(Cdata == -1, i) = CdataCols(1, i);
    CdataBins(Cdata >= 0, i) = CdataCols(2, i);
    CdataBins(Cdata >= 0.2, i) = CdataCols(3, i);
    CdataBins(Cdata >= 0.4, i) = CdataCols(4, i);
    CdataBins(Cdata >= 0.6, i) = CdataCols(5, i);
    CdataBins(Cdata >= 0.8, i) = CdataCols(6, i);
end

% check not-used bins and remove them
check_CdataBins = unique(join(string(CdataBins), "-"));
check_CdataCols = join(string(CdataCols), "-");
[~, nonUsedIdx] = setdiff(check_CdataCols, check_CdataBins);
lbls(nonUsedIdx) = [];
CdataCols(nonUsedIdx, :) = [];

lbls = ["Lead variant", flip(lbls)];

[otherhits, tophit_idx] = deal(true(height(sumstat), 1));
otherhits(tophit) = false;
tophit_idx(tophit) = false;
 
sc2 = scatter(sumstat.pos(otherhits),...
    sumstat.p(otherhits), 'MarkerFaceColor', ...
    'flat', 'MarkerEdgeColor', 'k','CData', CdataBins(tophit_idx, :),...
    'MarkerEdgeAlpha', p.MarkerEdgeAlpha,...
    'MarkerFaceAlpha', p.MarkerFaceAlpha, ...
    'SizeData', p.markersize);
dtrow = dataTipTextRow('snp', sumstat.snp(otherhits));
sc2.DataTipTemplate.DataTipRows(end+1) = dtrow;
try % if beta is present
    sc2.DataTipTemplate.DataTipRows(1).Label = 'beta';
    sc2.DataTipTemplate.DataTipRows(1).Value = sumstat.beta(otherhits);
catch
    sc2.DataTipTemplate.DataTipRows(1).Label = 'pos';
end
sc2.DataTipTemplate.DataTipRows(2).Label = 'p';
sc2.DataTipTemplate.DataTipRows(2).Value = 10.^-sumstat.p(otherhits);

% add extra lead variants (to compare two GWAS together). Extralead can be
% lead variants in another GWAS
if any(sumstat.extralead)
    scatter(sumstat.pos(sumstat.extralead), sumstat.p(sumstat.extralead), ...
        'CData', CdataBins(sumstat.extralead, :), ...
        'MarkerFaceColor', 'flat', ...
        'MarkerEdgeColor', 'k',...
        'MarkerEdgeAlpha', p.MarkerEdgeAlpha, ...
        'MarkerFaceAlpha', p.MarkerFaceAlpha, 'Marker', 'square',...
        'SizeData', p.markersize.*3);
end

% add lead variant --------------------------------------------------------
sc1 = scatter(sumstat.pos(tophit), sumstat.p(tophit), 'MarkerFaceColor', 'k',...
        'MarkerEdgeColor', 'k',...
        'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', 0.9, 'Marker', 'diamond',...
        'SizeData', p.markersize.*2);

try % if beta is present
    sc1.DataTipTemplate.DataTipRows(1).Label = 'beta';
    sc1.DataTipTemplate.DataTipRows(1).Value = sumstat.beta(tophit);
catch
    sc1.DataTipTemplate.DataTipRows(1).Label = 'pos';
end
sc1.DataTipTemplate.DataTipRows(2).Label = 'p';
sc1.DataTipTemplate.DataTipRows(2).Value = 10.^-sumstat.p(tophit);


% set X-axis limits of the regional plot
ax2.XLim = [min(sumstat.pos) - 100, max(sumstat.pos) + 100];
if ax2.XLim(1) < 0
    ax2.XLim(1) = 0;
end

CdataCols = [[0 0 0]; flip(CdataCols)];
lg = (0);
for i = 1:size(CdataCols, 1)
    if i == 1
        mrkr = 'diamond';
    else
        mrkr = 'O';
    end
    lg(i) = plot(NaN, NaN, 'MarkerFaceColor', CdataCols(i,:),...
            'LineStyle', 'none', 'Marker', mrkr,...
            'MarkerEdgeColor',  'k',...
            'Color',  CdataCols(i,:), 'MarkerSize', p.markersize);
end
leg = legend(lg, lbls, 'Color', [0.93 0.93 0.93], 'AutoUpdate', 'off',...
    'FontSize', p.legfontsize, 'Location', 'northeast', ...
    'Orientation', 'Vertical', 'FontName', p.fontname);
% leg = modifyLegendPos(leg, ax2);

%--------------------------------------------------------------------------
% set tophit label
tophit_label = tophit;
if p.labelExtralead
    tophit_label = union(tophit, find(sumstat.extralead));
end
lead_variants = sumstat.snp(tophit_label);

if p.causalLead % lead variant is a putative causal variant (*)
    topvariant = sumstat.snp(tophit);
    lead_variants(ismember(lead_variants, topvariant)) = topvariant + "(*)";
end

txth = text(sumstat.pos(tophit_label), sumstat.p(tophit_label),...
    lead_variants, 'FontAngle', 'italic',...
    'HorizontalAlignment', 'center', 'Interpreter', 'none',...
    'VerticalAlignment', 'bottom', 'FontSize', p.topfontsize, ...
    'FontName', p.fontname);
for k = 1:numel(txth)
    txth(k).Units = 'normalized';
    txth(k).Position(2) = txth(k).Position(2) + 0.05*ax2.Position(4);
    drawnow;
    txth(k).Units = 'data';
    if (txth(k).Extent(2) + txth(k).Extent(4)) > ax2.YLim(2)
       ax2.YLim(2) = txth(k).Extent(2) + txth(k).Extent(4);
    end
    drawnow;
    txth(k).Units = 'normalized'; 
end

% check if legend overlaps the text ---------------------------------------
checkLeadSNPoverlap(leg, txth);

p.significance = string(p.significance);
if p.significance == "bonferroni"
    p.significance = 0.05/(numel(sumstat.p));
elseif p.significance == "genome-wide"
    p.significance = 5e-8;
elseif p.significance == "fdr"
    checkp = 10.^(-sumstat.p);
    p.significance = mafdr(checkp, 'BHFDR', true);
    p.significance = max(checkp(p.significance <= 0.05)); % max p-value passing BH FDR cutoff
    if isempty(p.significance); p.significance = 0; end
else
    p.significance = 0;
end
if any(sumstat.p >= -log10(p.significance)) && (p.significance ~= 1)
    line(ax2.XLim, [-log10(p.significance), -log10(p.significance)],...
        'LineWidth', 1, ...
        'color', 'r', 'LineStyle', '--') % Significance threshold
end
set(ax2,'box','off','LineWidth', 1.2, 'TickLength', [0.001 0], ...
    'FontSize', 12, 'FontWeight', 'bold', 'TickDir', 'out', ...
    'FontName', p.fontname)

ax2.XTickLabel = string(ax2.XTick./1e6);
ax2.XAxis.FontSize = p.xfontsize;
xlabel(['chromosome ', p.chr, ' (Mb)'], 'FontSize', p.xfontsize + 2, ...
    'FontWeight', 'bold', 'FontName', p.fontname)
ylabel('-log_{10} (\itp)', 'FontSize', p.yfontsize + 3, ...
    'FontWeight', 'bold', 'FontName', p.fontname)
ax2.YAxis.FontSize = p.yfontsize;
if p.box % box around the figure
    ax2.Box = 'on';
else
    ax2.Box = 'off';
end

% ensure that y-label is within axis borders
drawnow;
ax2.Position(1) = ax2.TightInset(1);

% set gene plot axis limits the same as gwas plot
drawnow
ax1.Position(1) = ax2.Position(1);
ax1.XLim = ax2.XLim;

% check if new ax1 limits covers gene names
tobj = findobj(ax1, 'Type', 'Text');
for i = 1:numel(tobj)
    geneIdx = strcmp({r.symbol}, replace(tobj(i).String, ["\leftarrow", "\rightarrow", ">", "<"], "")); % is this gene even visible on ax1? if not, remove it's label too
    if (r(geneIdx).stop > ax1.XLim(2) && r(geneIdx).start > ax1.XLim(2))...
            || (r(geneIdx).start < ax1.XLim(1) && r(geneIdx).stop < ax1.XLim(1))
        continue
    end
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

% check if text objects (gene names) overlap; if so, reduce their font size
checkTextObjOverlap(tobj); 

grid(ax2, 'on')
ax2.GridAlpha = 0.025;

% add fine-mapping legend -------------------------------------------------
if p.finemap
    ax3 = copyobj(ax2, gcf);
    ax3.Visible = 'off';
    leg3 = legend(ax3, sc3, 'Color', [0.93 0.93 0.93], ...
        'AutoUpdate', 'off', 'FontSize', p.legfontsize + 2, ...
        'Location', 'northwest', 'Orientation', 'Vertical');
    leg3.String = string(csgrp);
    leg3.Title.String = "Credible set";
    leg3.Title.FontName = p.fontname;

    drawnow;
    for i = 1:numel(ax3.Children)
        ax3.Children(i).Visible = 'off';
    end
    drawnow;
end

if p.title == "" || isempty(p.title)
    p.title = getRandomName("genePlotter", 3);
else
    if ~startsWith(p.title, "\it{") && ~endsWith(p.title, "}")
        p.showtitle = "\it{" + p.title + "}";
    end
    title(ax2, p.showtitle, 'Interpreter', 'tex', ...
        'FontSize', p.yfontsize + 4, 'FontName', p.fontname);
end

if p.save
    exportgraphics(ax2.Parent, fullfile(p.outdir, matlab.lang.makeValidName(p.title) ...
        + "." + p.format), 'Resolution', p.resolution)
end

if p.savefig
    savefig(ax2.Parent, fullfile(p.outdir, matlab.lang.makeValidName(p.title) + ".fig"))
end
end % END

%% subfunctions ===========================================================
function [leg, txth] = checkLeadSNPoverlap(leg, txth) 
overlap = true;
drawnow
triedWays = 0;
triedAdjustment = {'right', 'left'};
% sc = findobj(txth.Parent.Children, "Type", "Scatter"); % scatter points
while overlap
    pos1 = txth.Extent;
    pos2 = leg.Position;
    thresh = 0.08;
    
    % check x-position
    xoverlap = false;
    if pos2(1) > pos1(1)
        if (pos1(1) + pos1(3)) > pos2(1) || abs(pos1(1) + pos1(3) - pos2(1)) <= thresh
            xoverlap = true;
        end
    else
        if (pos2(1) + pos2(3)) > pos1(1) || abs(pos2(1) + pos2(3) - pos1(1)) <= thresh
            xoverlap = true;
        end
    end
    
    % check y-position
    yoverlap = false;
    if pos2(2) > pos1(2)
        if (pos1(2) + pos1(4)) > pos2(2) || abs(pos1(2) + pos1(4) - pos2(2)) <= thresh
            yoverlap = true;
        end
    else
        if (pos2(2) + pos2(4)) > pos1(2) || abs(pos2(2) + pos2(4) - pos1(2)) <= thresh
            yoverlap = true;
        end
    end
    
    if ~all([xoverlap, yoverlap]) % no overlap
        overlap = false; 
    else % change leg position
        if triedWays == 0 % first change legend location (default is northeast)
            leg.Location = 'northwest';
            triedWays = triedWays + 1;
        else % next try to change text adjustment 
            txth.HorizontalAlignment = triedAdjustment{triedWays};
            triedWays = triedWays + 1;
        end
    end
    drawnow    
end

end 

%% ------------------------------------------------------------------------
function checkTextObjOverlap(tobj)

drawnow

for i = 1:numel(tobj)
%     xex1 = tobj(i).Extent(1) + tobj(i).Extent(3);
%     x1 = tobj(i).Extent(1);
%     yex1 = tobj(i).Extent(2) + tobj(i).Extent(4);
%     y1 = tobj(i).Extent(2);
    for j = i+1:numel(tobj)
%         xex2 = tobj(j).Extent(1) + tobj(j).Extent(3);
%         x2 = tobj(j).Extent(1);
%         yex2 = tobj(j).Extent(2) + tobj(j).Extent(4);
%         y2 = tobj(j).Extent(2);

        while (((tobj(i).Extent(1) + tobj(i).Extent(3)) > tobj(j).Extent(1) ...
                && tobj(j).Extent(1) > tobj(i).Extent(1)) || ...
                ((tobj(j).Extent(1) + tobj(j).Extent(3)) > tobj(i).Extent(1) ...
                && tobj(j).Extent(1) < tobj(i).Extent(1))) && ...
                (tobj(i).Extent(2) == tobj(j).Extent(2))% x-ax overlap
            drawnow
            tobj(j).FontSize = tobj(j).FontSize - 1;
            tobj(i).FontSize = tobj(i).FontSize - 1;
        end

    end
end

end