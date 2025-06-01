function regionalPlotter(opts)

arguments
    % either the name of summary stat file or a mat file of regions around
    % each index variant. In case of a mat file, it must be a struct with
    % two fields: gss and index. in either case when 'usebackup' is true
    % ('gss' can be left empty in this case) the function first tries to
    % find it in 'outdir'.
    opts.gss 

    % either a table or mat file of index variants for 'gss'. If finds
    % index variants in the 'gss' mat file, this is not mandatory. Note
    % that is 'usebackup' is true, this field won't be used.
    opts.index 
    opts.finemap % results of fine-mapping (either a table or a mat file)
    opts.pipfilter (1,1) double = nan % cutoff to filter PIPs below this. default is nan: keeps max PIP per each CS
    opts.win (1,1) double = 500 % in kbp. Window around the index variants 
    opts.workers (1,1) = nan
    opts.parallel (1,1) logical = true % use tall/gather for reading summary stat file
    opts.outdir = fullfile(pwd, "RegionalPlots")
    opts.savegss (1,1) logical = true % to save index/gss around each index variant to a file in 'outdir' 
    opts.usebackup (1,1) logical = false % to use saved gss/index/finemap data. This flag overrides 'gss', 'index' but not 'finemap'.
    opts.findlocus (1,1) logical = true % to find the nearest gene for index variants
    opts.refGenome {mustBeTextScalar, mustBeMember(opts.refGenome, ["GRCh37", "GRCh38"])} =  "GRCh37" % note: opentargets only retruns data from GRCh38 (for mapVariant2gene)

    % genePlotter
    opts.significance {mustBeTextScalar, mustBeMember(opts.significance, ["genome-wide", "bonferroni", "fdr", "none"])} = "genome-wide" 
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "jpg" % format of output plots
    opts.resolution (1, 1) double = 300 % resolution of output plot
    opts.fontname {mustBeTextScalar} = "Garamond"
    opts.musthaveUTR (1,1) logical = true % to visualize only genes which have UTR (and CDS)
end

% check parallel pool
if opts.parallel && isempty(gcp('nocreate'))
    if isnan(opts.workers)
        % number of physical cores usually gives the best efficiency for
        % regression analyses. Although this has not been thoroughly
        % tested.
        maxNumCompThreads('automatic'); % reset (only if was changed in current session)
        opts.workers = maxNumCompThreads; 
    end
    parc = parcluster('Processes');
    parc.NumWorkers = opts.workers;
    parpool(parc); 
end

if ~isfolder(opts.outdir), mkdir(opts.outdir); end

gssfile = fullfile(opts.outdir, "gss_index_regional.mat");
if ~isfile(gssfile) && ~isfield(opts, 'gss')
    error("cannot find gss and index files within 'outdir' dir!")
elseif isfile(gssfile)
    disp("found gss_index_regional.mat file!")
    opts.gss = gssfile;
end

% use backed up files already saved to the 'outdir'
if opts.usebackup && endsWith(opts.gss, ".mat")
    gss = load(opts.gss, "gss").gss;
    opts.index = load(opts.gss, "index").index;

    if ~isfield(opts, 'finemap')
        try
            opts.finemap = load(opts.gss, "finemap").finemap;
        catch
        end
    end
end

if isfield(opts, 'index') && ~istable(opts.index)
    index_tmp = load(opts.index);
    fis = string(fieldnames(index_tmp));
    opts.index = index_tmp.(fis(1));
end

% add locus (nearest gene) to index table
locusCol = contains(lower(colnames(opts.index)), ["gene", "locus"]);
if any(locusCol)
    locusCol = find(locusCol, 1);
    opts.index.locus = opts.index.(locusCol);
end

if opts.findlocus
    cols = createGWASheader(colnames(opts.index)); % create a unified header
    snpcol = find(cols == "SNP");
    if ~any(colnames(opts.index) == "locus")
        vep = mapVariant2gene(opts.index, verbose=false, refGenome=opts.refGenome);
        [f1, f2] = ismember(opts.index.(snpcol), vep.id);
        opts.index.locus(:) = "-";
        opts.index.locus(f1) = vep.nearestGene(f2(f1));
    end
end

% fetch fine-mapping results ----------------------------------------------
if isfield(opts, 'finemap') && ~istable(opts.finemap)
    opts.finemap = load(opts.finemap);
    fis = string(fieldnames(opts.finemap));
    opts.finemap = opts.finemap.(fis(1));
end

% fetch gwas summary stats for all loci -----------------------------------
if ~endsWith(opts.gss, ".mat")
    if ~isfield(opts, 'index')
        error("index is missing!")
    end
    opts.index.Properties.VariableNames = createGWASheader(colnames(opts.index));
    [gss, ~, index] = readGWASfile(opts.gss, index=opts.index, ...
        parallel=opts.parallel, light=false, full=true, n=true, ...
        win=opts.win, legacy=false);

    if isfield(opts, 'finemap')
        finemap = opts.finemap;
    end

    if opts.findlocus
        cols = createGWASheader(colnames(opts.index)); % create a unified header
        snpcol = find(cols == "SNP");
        [f1, f2] = ismember(index.snp, opts.index.(snpcol));
        index.locus(:) = "-";
        index.locus(f1) = opts.index.locus(f2(f1));
    end

    opts.index = index;

    if opts.savegss
        save(fullfile(opts.outdir, "gss_index_regional.mat"), "gss", "index", "finemap")
    end
end

% match fine-mapped variants to gss variants
if isfield(opts, 'finemap')
    opts.finemap.id = opts.finemap.CHR + "." + opts.finemap.BP + "." + ...
        opts.finemap.A1 + "." + opts.finemap.A2;
    gss.id = gss.chr + "." + gss.pos + "." + gss.allele1 + "." + gss.allele0;
    [f1, f2] = ismember(gss.id, opts.finemap.id); 
    if any(~f1) % flip A1/A2
        gss.id = gss.chr + "." + gss.pos + "." + gss.allele0 + "." + gss.allele1;
        [f1f, f2f] = ismember(gss.id, opts.finemap.id);
        if sum(f1f) > sum(f1)
            f1 = f1f; f2 = f2f;
        end
    end
    gss = gss(f1, :);
    opts.finemap = opts.finemap(f2(f1), ["SNP", "CREDIBLE_SET", "PIP"]);
    opts.finemap.Properties.VariableNames = ["SNP", "cs", "pip"];
end

% loop over each region ---------------------------------------------------
index = opts.index;
for k = 1:height(index)

    idx = gss.chr == index.chr(k) & gss.pos >= (index.pos(k) - opts.win*1e3)...
        & gss.pos <= (index.pos(k) + opts.win*1e3);
    
    tab = gss(idx, :);
    if isempty(tab), continue; end
    region = tab.chr(1) + "_" + min(tab.pos) + "_" + max(tab.pos);
    if opts.findlocus
        region = index.locus(k) + "_" + region;
    end
    
    gopts = struct;
    if isfield(opts, 'finemap')
        fm = opts.finemap(idx, :);

        if isnan(opts.pipfilter)
            fm = groupfilter(fm, "cs", @(x) x == max(x), "pip");
        else
            % single PIP cut-off 
            fm_max = groupfilter(fm, "cs", @(x) x == max(x), "pip");
            fm(fm.pip < opts.pipfilter, :) = [];
            if isempty(fm), fm = fm_max; end
        end

        gopts.finemap = fm;

        % check if there are more than 1 putative causal variant (PIP >= cutoff)
        idx = fm.pip >= opts.pipfilter;
        if sum(idx) >= 1
            tmp = tab(ismember(tab.snp, fm.SNP(idx)), :);
            % tmp = sortrows(tmp, "p", "ascend");
            % gopts.extralead = tmp.snp(2:end);
            gopts.extralead = setdiff(tmp.snp, index.snp(k));
            gopts.labelExtralead = true;
            
            idx = fm.SNP == index.snp(k);
            if any(idx) % lead variant is amongst putative causal variants
                gopts.causalLead = true;
            end
        end
    end
    
    gopts.ldmethod="insample";
    gopts.backup=true; 
    gopts.parallel=opts.parallel;
    if opts.findlocus
        gopts.title=index.locus(k);
    else
        gopts.title=" ";
    end
    gopts.save=true;
    gopts.savefig=true;
    gopts.resolution=opts.resolution;
    gopts.significance=opts.significance;
    gopts.outdir=fullfile(opts.outdir, region);
    gopts.format=opts.format;
    gopts.lead = index.snp(k);
    gopts.fontname = opts.fontname;
    gopts.musthaveUTR = opts.musthaveUTR;
    gopts = namedargs2cell(gopts);
    genePlotter(tab, gopts{:})
    
end

end