function res = pruneFinemapRes(opts)
% the idea is to prune/merge fine mappings from FIENAMP and SuSiE based on
% a cutoff on PIP, so that:
%           1- in the first step, only PIPs above this cutoff (0.05) are
%           kept, and 
%           2- in the next step, all variants in 1st CS and the variant
%           with max PIP in other CSs will be kept. The number of variants
%           per each CS will be shown as well. If there is no 1st CS (for
%           SuSiE due to purity adjustment), top PIPs (> 0.05) in 0 CS will
%           be kept.

arguments
    opts.finemap {mustBeFile}
    opts.susie {mustBeFile}
    opts.pipthresh (1,1) double = 0.05 % remove variants below this PIP
end


if ~any(isfield(opts, {'susie', 'finemap'}))
    res = [];
    return
elseif all(isfield(opts, {'susie', 'finemap'}))
    fm = pruneFinemap(opts.finemap, opts);
    ss = pruneFinemap(opts.susie, opts);
elseif isfield(opts, 'finemap')
    res = pruneFinemap(opts.finemap, opts);
    return
else % SuSiE
    res = pruneFinemap(opts.susie, opts);
    return
end

regions = union(fm.region, ss.region);
infocols = ["SNP", "CHR", "BP", "A1", "A2", "SNPVAR", "Z", "region", "merge"];
cols = array2table(["PIP", "PIP.FINEMAP", "PIP.SuSiE"; ...
    "BETA_MEAN", "BETA.FINEMAP", "BETA.SuSiE"; "BETA_SD", ...
    "SD.FINEMAP", "SD.SuSiE"; "CS", "CS.FINEMAP", ...
    "CS.SuSiE"; "CS size", "CS.size.FINEMAP", "CS.size.SuSiE"]);
tabo = array2table(strings(1, numel(infocols) + numel(cols.(2))*2), ...
        'VariableNames', [infocols, union(cols.(2), cols.(3), 'stable').']);
tabo = convertvars(tabo, ["CHR", "BP", "SNPVAR", "Z", "merge", union(cols.(2), cols.(3)).'], @double);

res = cell(numel(regions), 1);
for i = 1:numel(regions)
    fidx = fm.region == regions(i);
    sidx = ss.region == regions(i);

    ftab = fm(fidx, :); stab = ss(sidx, :);
    snp = union(ftab.SNP, stab.SNP);
    tab = repmat(tabo, numel(snp), 1);
    if isrow(snp); snp = snp.'; end
    tab.SNP = snp;

    [fidx1, fidx2] = ismember(snp, ftab.SNP); fidx2(fidx2 < 1) = [];
    [sidx1, sidx2] = ismember(snp, stab.SNP); sidx2(sidx2 < 1) = [];

    for j = 1:height(cols)
        tab.(cols.(2)(j))(fidx1) = ftab.(cols.(1)(j))(fidx2);
        tab.(cols.(3)(j))(sidx1) = stab.(cols.(1)(j))(sidx2);
    end

    for j = 1:numel(infocols)
        tab.(infocols(j))(sidx1) = stab.(infocols(j))(sidx2);
        try % FINEMAP from polyFun doesn't have all columns
            tab.(infocols(j))(fidx1) = ftab.(infocols(j))(fidx2);
        catch 
        end
    end

    res{i} = tab;

end

res = vertcat(res{:});
res = sortrows(res, {'CHR', 'BP'});
end % END

%% subfctions =============================================================
function res = pruneFinemap(fmap, opts)
fmap = load(fmap);
fi = fieldnames(fmap);
fmap = fmap.(fi{1});

% remove variants that don't belong to any CS and have very low PIP. The
% latter filter keeps high PIPs with CS0, for which SuSiE could not create
% a 95% CS because of purity criteria. 
fmap(endsWith(fmap.CREDIBLE_SET, ":0") & fmap.PIP < opts.pipthresh, :) = [];

% unique regions
cs = unique(erase(unique(fmap.CREDIBLE_SET), ":" + digitsPattern + textBoundary("end")));
regions = array2table(split(cs, [":", "-"]), 'VariableNames', {'chr', 'start', 'end'});
regions = convertvars(regions, {'start', 'end'}, @double);
regions.chr = double(replace(regions.chr, "chr", ""));
regions.cs = cs;
regions.merge = nan(height(regions), 1);
regions = sortrows(regions, 'start');
regions = sortrows(regions, 'chr');
chr = unique(regions.chr);

grp = 1;
for i = 1:numel(chr)
    ridx = regions.chr == chr(i);
    rtab = regions(ridx, :);

    interval =  reshape(rtab{:, 2:3}',1,[]);
    ii = strfind(sign(diff(interval)), [-1 1]);
    interval([ii,ii+1]) = [];
    interval = reshape(interval, 2, [])';
    
    for j = 1:size(interval, 1)
        idx = rtab.start >= interval(j, 1) & rtab.end <= interval(j, 2);
        rtab.merge(idx) = grp;
        grp = grp + 1;
    end
    
    regions(ridx, :) = rtab;
end

% merge overlapping regions: not recommended, because may merge independent
% loci
mergeIdx = unique(regions.merge);
res = cell(numel(mergeIdx), 1);
for i = 1:numel(mergeIdx) 

    regIdx = find(regions.merge == mergeIdx(i));
    region = "chr" + regions.chr(regIdx(1)) + ":" + regions.start(regIdx(1)) + "-" + regions.end(regIdx(end));
    thisCS = regions.cs(regIdx);
    tab = fmap(startsWith(fmap.CREDIBLE_SET, thisCS), :);
    tab = sortrows(tab,'PIP','descend');
    tab.region = erase(tab.CREDIBLE_SET, ":" + digitsPattern + textBoundary);
    tab.CS = extract(tab.CREDIBLE_SET, ":" + digitsPattern + textBoundary);
    tab.CS = erase(tab.CS, ":").double;
    tab.locus = repmat(region, height(tab), 1);
    cscount = groupsummary(tab, 'CREDIBLE_SET');
    % don't merge 
    % tab.CREDIBLE_SET = extract(tab.CREDIBLE_SET, ":" + digitsPattern + textBoundary("end"));
    % tab.CREDIBLE_SET = double(erase(tab.CREDIBLE_SET, ":"));

    % max PIP per each CS and keep all SNPs in first CS
    idx = endsWith(tab.CREDIBLE_SET, ":1") & tab.PIP >= opts.pipthresh;
    if ~any(idx) % can happen with susie due to purity threshold, keep 0 then
        idx = endsWith(tab.CREDIBLE_SET, ":0") & tab.PIP >= opts.pipthresh;
    end
    tabr = tab(~idx, :);
    res{i} = groupfilter(tabr, 'CREDIBLE_SET', @(x) x == max(x), 'PIP');
    res{i} = [tab(idx, :); res{i}];
    res{i}.merge = repmat(mergeIdx(i), height(res{i}), 1);
    
    [f1, f2] = ismember(res{i}.CREDIBLE_SET, cscount.CREDIBLE_SET);
    res{i}.("CS size") = nan(height(res{i}), 1);
    res{i}.("CS size")(f1) = cscount.GroupCount(f2);
    res{i} = sortrows(res{i}, {'CS','CREDIBLE_SET'});

end
res = vertcat(res{:});

end 