function [tab, failed] = liftover(tab, opts)
% performs liftover using rtracklayer package
% Oveis Jamialahmadi, Sahlgrenska Academy, Jan 2022.
% 
% @24/01/2022: 'useEnsemble' was added to perform mapping with ENSEMBL REST
%               if liftOver fails to map some positions.
% @30/12/2022: a bug was fixed: string chr are forced to be numeric.
% @20NOV2024: 'method' option was added: ensembl, ucsc or both. If anything
%             but uscs is selected, 'useEnsemble' is not recommended,
%             because it's equivalent to 'ensembl' chain (REST API 'map'
%             method uses the same approach). If provided, ignores 'chain'
%             argument.

arguments
    tab {mustBeA(tab, 'table')}
    opts.chain {mustBeFile} = fullfile(fileparts(which('liftover.m')), 'chain', 'GRCh37_to_GRCh38.chain')
    opts.useEnsemble (1,1) logical = true % to use ENSEMBL REST API for variants liftOver fails to map
    opts.backup (1,1) logical = false % save liftovers locally for future use 
    opts.method {mustBeTextScalar, mustBeMember(opts.method, ["ensembl", "ucsc", "both"])} % overrides 'chain'
end

if opts.backup % under development
    opts.backupdir = fullfile(fileparts(which('liftover.m'), 'backup'));
    if ~isfolder(opts.backupdir)
        mkdir(opts.backupdir)
    end
    opts.backupfiles = string({dir(fullfile(opts.backupdir, '*.mat')).name}.');
    opts.backupfiles = opts.backupfiles(startsWith(opts.backupfiles, 'liftover.'));
    opts.backupfiles = fullfile(opts.bakupdir, opts.backupfiles);
end

file = getRandomName("lift", 3);
cols = lower(colnames(tab));
posIdx = find(ismember(cols, {'pos', 'bp', 'pos.grch37', 'pos.hg19', 'genpos', 'pos37'}));
chrIdx = find(ismember(cols, {'chr', 'ch', 'chrom'}));
pos = tab.(posIdx(1)); chr = double(tab.(chrIdx(1)));
save(file + ".mat", 'pos', 'chr')

%@20NOV2024: 'method'
if isfield(opts, "method")
    chfiles = getfilenames(fullfile(fileparts(which('liftover.m')), 'chain'), "chain", fullpath=true).chain;
    assert(~isempty(chfiles), "liftover:no chain file was found in chain subfolder of liftover function!")
    
    % chfileE = chfiles(chfiles.contains("GRCh" + digitsPattern + "_to_GRCh" + digitsPattern));
    chfileE = chfiles(chfiles.contains("GRCh37_to_GRCh38"));
    chfileU = chfiles(chfiles.contains("hg19ToHg38.over"));

    if opts.method == "ensembl"
        assert(~isempty(chfileU), "liftover:no ensembl chain file was found in chain subfolder of liftover function!")
        opts.chain = chfileE;
    elseif opts.method == "ucsc"
        assert(~isempty(chfileU), "liftover:no ucsc chain file was found in chain subfolder of liftover function!")
        opts.chain = chfileU;
    else
        assert(~isempty(chfileU), "liftover:no ensembl chain file was found in chain subfolder of liftover function!")
        assert(~isempty(chfileU), "liftover:no ucsc chain file was found in chain subfolder of liftover function!")
        opts.chain = [chfileE; chfileU];
    end

end

opts.chain = string(replace(opts.chain, "\", "/"));
r = string;
r(1) = 'if("rtracklayer" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rtracklayer")}';
r(numel(r) + 1) = "library(rtracklayer)";

if numel(opts.chain) > 1
    for k = 1:numel(opts.chain)
        r(numel(r) + 1) = "chainObject" + k + " <- import.chain('" + opts.chain(k) + "')";
    end
else
    r(numel(r) + 1) = "chainObject <- import.chain('" + opts.chain + "')";
end

r(numel(r) + 1) = "mat <- R.matlab::readMat('" + file + ".mat')";
r(numel(r) + 1) = "grObject <- GRanges(seqnames=paste0('chr', mat$chr), ranges=IRanges(start=mat$pos, end=mat$pos))";

if numel(opts.chain) > 1
    for k = 1:numel(opts.chain)
        r(numel(r) + 1) = "mapped" + k +" <- as.data.frame(liftOver(grObject, chainObject" + k + "))";
        r(numel(r) + 1) = "mapped" + k + " <- mapped" + k + "[, c('group', 'start')]";
    end
    r(numel(r) + 1) = "mapped = list(" + join("mapped" + (1:numel(opts.chain)), ",") + ")";
    r(numel(r) + 1) = "mapped <- Reduce(function(x, y) merge(x, y, by = 'group', all = TRUE), mapped)";

else
    r(numel(r) + 1) = "mapped <- as.data.frame(liftOver(grObject, chainObject))";
    r(numel(r) + 1) = "mapped <- mapped[, c('group', 'start')]";
end
r(numel(r) + 1) = "write.table(mapped, file = '" + file + ".txt', quote = F, row.names = F)";


MATLAB2Rconnector(file + ".r", 'delr', true, 'code', r');

if isempty(gcp("nocreate"))
    map = readtable(file + ".txt", 'VariableNamingRule', 'preserve');
else
    map = tabularTextDatastore(file + ".txt", TextType="string", ...
        VariableNamingRule="preserve", TreatAsMissing="NA");
    map = gather(tall(map));
end

if numel(opts.chain) > 1
    cols = colnames(map);
    cols = cols(cols.startsWith("start"));
    assert(numel(cols) == 2, "something went wrong with rtracklayer liftOver function (less than 2 mapped cols returned!)")

    % keep the consistent cols
    map = renamevars(map, cols, ["en", "uc"]);
    map.start(:) = nan;
    idx = map.en == map.uc;
    map.start(idx) = map.en(idx);
    map.method(:) = "Ensembl&UCSC";

    % find inconsistent cols and prefer Ensembl over UCSC
    nan_idx_en = ismissing(map.en);
    nan_idx_uc = ismissing(map.uc);
    fprintf("\t%d failed liftovers for Ensembl\n", sum(nan_idx_en))
    fprintf("\t%d failed liftovers for UCSC\n", sum(nan_idx_uc))
    fprintf("\t%d failed liftovers for Ensembl and UCSC\n", sum(nan_idx_en & nan_idx_uc))
    nan_idx = nan_idx_en | nan_idx_uc;

    inc_idx = (map.en ~= map.uc) & (~nan_idx);
    fprintf("\t%d inconsistent liftovers between Ensembl and UCSC\n", sum(inc_idx))
    fprintf("\t\tkept Ensembl locations for inconsistent maps\n")
    assert(all(ismissing(map.start(inc_idx))))
    map.start(inc_idx) = map.en(inc_idx);
    map.method(inc_idx) = "Ensembl(conflict)";


    % fill the rest with whatever exclusively found by UCSC or Ensembl
    exc_idx = nan_idx_en & ~nan_idx_uc;
    fprintf("\t%d sites were filled with UCSC exclusively\n", sum(exc_idx))
    assert(all(ismissing(map.start(exc_idx))))
    map.start(exc_idx) = map.uc(exc_idx);
    map.method(exc_idx) = "UCSC";

    exc_idx = ~nan_idx_en & nan_idx_uc;
    fprintf("\t%d sites were filled with Ensembl exclusively\n", sum(exc_idx))
    assert(all(ismissing(map.start(exc_idx))))
    map.start(exc_idx) = map.en(exc_idx);
    map.method(exc_idx) = "Ensembl";

    map(:, ["uc", "en"]) = [];
    nan_idx = ismissing(map.start);
    if any(nan_idx)
        fprintf("\t%d sites failed in the end\n", sum(nan_idx))
    end

end

if numel(opts.chain) > 1
    opts.multichain = true;
else
    opts.multichain = false;
end

delete(file + ".mat")
delete(file + ".txt")

posCol = string(tab.Properties.VariableNames(posIdx(1)));
tab.(posCol + "_old") = tab.(posCol);

if opts.multichain
    assert(~any(colnames(tab) == "liftover"), "cannot add liftover column to input data")
    tab.liftover(:) = "";
    tab.liftover(tab.liftover == "") = missing;
end

if size(tab, 1) ~= size(map, 1) % some unmapped variants
    tab.(posCol) = nan(size(tab, 1), 1);
    if isempty(map) % nothing could be mapped
        failed = tab;
    else
        tab.(posCol)(map.group) = map.start;

        if opts.multichain
            tab.liftover(map.group) = map.method;
        end

        nanidx = isnan(tab.(posCol));
        failed = tab(nanidx, :);
    end
else
    failed = [];
    tab.(posCol) = map.start;

    if opts.multichain, tab.liftover = map.method; end
end

if opts.useEnsemble && ~isempty(failed)
    fprintf("\tTrying ENSEMBL API for %d failed variants\n", height(failed))
    for i = 1:numel(failed.(posCol + "_old"))
        try
            tmp = EnsemblREST(failed.(chrIdx)(i) + ":" + failed.(posCol + "_old")(i),...
                'map', 'refGenome', '37');
            failed.(posCol)(i) = tmp.start;
            if opts.multichain, failed.liftover(i) = "EnsemblAPI"; end
            fprintf("\t\t%d of %d: success\n", i, height(failed))
        catch % failed to map
            % fprintf("\t\t%d of %d:failed\n", i, height(failed))
        end
    end

    % fill tab table
    tab(nanidx, :) = failed; 
    failed(~isnan(failed.(posCol)), :) = [];
    if isempty(failed); failed = []; end

    if ~isempty(failed) % try for the second time (API may fail)
        nanidx = isnan(tab.(posCol));
        for i = 1:numel(failed.(posCol + "_old"))
            try
                tmp = EnsemblREST(failed.(chrIdx)(i) + ":" + failed.(posCol + "_old")(i),...
                    'map', 'refGenome', '37');
                failed.(posCol)(i) = tmp.start;
                if opts.multichain, failed.liftover(i) = "EnsemblAPI"; end
            catch % failed to map
            end
        end

        % fill tab table
        tab(nanidx, :) = failed; 
        failed(~isnan(failed.(posCol)), :) = [];
        if isempty(failed); failed = []; end
    end

end


end % END