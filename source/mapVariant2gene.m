function [vep, bgi] = mapVariant2gene(raw, opts)
% find nearest coding gene to a variant and fetches the most severe
% consequence from VEP
% 
% @06APR2023: we keep protein coding genes for 'nearestGene' calculations.
%            Now, protein coding genes across GRCh37 and 38 are kept. This
%            was done because of an observation that some genes were
%            pseudogenes in 37 but protein coding in 38.

arguments
    raw {mustBeNonempty}
    opts.eqtl (1,1) logical = false
    opts.veponly (1,1) logical = false % only vep table returned from callVEP function
    opts.bgionly (1,1) logical = false % only bgi struct
    opts.refGenome {mustBeTextScalar, mustBeMember(opts.refGenome, ["GRCh37", "GRCh38"])} =  "GRCh37" % note: opentargets only retruns data from GRCh38 
    opts.verbose (1,1) logical = true
    opts.bgenhome {mustBeFolder}

    % GTF files: if available, the function uses them for finding gene
    % symbols, otherwise, uses Ensembl REST API
    opts.gtf37 {mustBeFile} = fullfile(fileparts(which("mapVariant2gene.m")), "GTF", "Homo_sapiens.GRCh37.87.gtf.gene.mat")
    opts.gtf38 {mustBeFile} = fullfile(fileparts(which("mapVariant2gene.m")), "GTF", "Homo_sapiens.GRCh38.107.gtf.gene.mat")
end

% for gnomAD
if opts.refGenome == "GRCh37"
    opts.gdataset = "gnomad_r2_1";
else
    opts.gdataset = "gnomad_r3";
end

% get BIM struct ----------------------------------------------------------
if istable(raw)
    cols = lower(raw.Properties.VariableNames);
    snpcol = find(ismember(cols, ["annotate", "snp", "id", "variant", "variantid", "variant_id"]), 1);
    chrcol = find(ismember(cols, ["chr", "chrom", "chromosome"]), 1);
    snp.snp = raw.(snpcol); snp.chr = string(raw.(chrcol));
elseif isstring(raw)
    snp.snp = raw;
    chr = EnsemblREST(snp.snp, 'vep', 'refGenome', erase(opts.refGenome, "GRCh"));
    if isstruct(chr)
        chr = string(chr.seq_region_name);
    else
        chr = cellfun(@(x) string({x.seq_region_name, x.id}), chr, 'uni', false);
        chr = vertcat(chr{:});
        [~, idx] = ismember(snp.snp, chr(:, 2)); idx(idx < 1) = [];
        chr = chr(idx, 1);
    end
    if isempty(chr); error('mapVariant2gene::cannot find chromosomes!'); end
    snp.chr = chr;
else
    error('mapVariant2gene::unknown input!')
end
bgi = bgireader(snp, 'bgenhome', opts.bgenhome, 'verbose', opts.verbose);
bgi = table2struct(struct2table(vertcat(bgi{:})), 'ToScalar', true);
strFields = ["chr", "a1", "a2", "pos", "snp"];
for i = 1:numel(strFields)
    bgi.(strFields(i)) = string(bgi.(strFields(i))); 
end

if isstring(bgi.chr)
    bgi.chr = regexprep(bgi.chr, "^0", "");
end

if opts.bgionly
    vep = [];
    return
end

% % VEP annotation ----------------------------------------------------------
vep = callVEP(bgi, 'method', 'vep', 'backup', true, ...
    'genomeRef', erase(opts.refGenome, "GRCh"), 'verbose', opts.verbose);

if opts.veponly
    return
end

vep_genes = vep.gene(vep.gene ~= "" & ~contains(vep.gene, ';')); 
if ~isempty(vep_genes)
    refTag = erase(opts.refGenome, "GRCh");
    if isfield(opts, "gtf"+refTag)
        genes = load(opts.("gtf"+refTag), "gene_name", "gene_biotype");
        genes = struct2table(genes);
        genes.Properties.VariableNames = regexprep(colnames(genes), "^gene_", "");
        genes(~ismember(genes.name, vep_genes), :) = [];

        % @06APR2023: see notes above
        otherRef = setdiff(["37", "38"], erase(opts.refGenome, "GRCh"));
        geneso = load(opts.("gtf"+otherRef), "gene_name", "gene_biotype");
        geneso = struct2table(geneso);
        geneso.Properties.VariableNames = regexprep(colnames(geneso), "^gene_", "");
        geneso(~ismember(geneso.name, vep_genes), :) = [];
        geneso(~ismember(geneso.name, genes.name), :) = []; % only keep the same genes in the otherRef file

        genes = unique([genes; geneso], "rows");
        dups = duplicates(genes.name);
        genes.rem(:) = false;
        for k = 1:numel(dups)
            idx = find(genes.name == dups(k));
            idx = idx(genes.biotype(idx) ~= "protein_coding");
            if isempty(idx) || numel(idx) > 1, idx = 1; end % otherwise pick one 
            genes.rem(idx) = true;
        end
        genes(genes.rem, :) = [];
        genes.rem = [];
    else
        genes = EnsemblREST(unique(vep_genes), 'lookup', 'geneSymbol', true, ...
            'refGenome', refTag);
    end
    vep_genes_rem = setdiff(vep_genes, genes.name);

    otherRef = setdiff(["37", "38"], erase(opts.refGenome, "GRCh"));
    if ~isempty(vep_genes_rem) % check other refGenome
        if ~isfield(opts, "gtf"+otherRef)
            genes = [genes; EnsemblREST(vep_genes_rem, 'lookup', ...
                'geneSymbol', true, 'refGenome', otherRef)];
        end
    end
    ncd_idx = ismember(vep.gene, ...
        genes.name(genes.biotype ~= "protein_coding")) | ...
        contains(vep.gene, ';') | vep.gene == "" | ismissing(vep.gene);
    
else % all snps are intergenic
    ncd_idx = true(height(vep), 1); 
end

if any(ncd_idx)
    keepbgi = bgi;
    for i = 1:numel(strFields)
        bgi.(strFields(i))(~ncd_idx) = []; 
    end

    % check GTF file
    refTag = erase(opts.refGenome, "GRCh");
    if isfield(opts, "gtf"+refTag)
        gtf = struct2table(load(opts.("gtf"+refTag)));
        gtf(gtf.gene_biotype ~= "protein_coding", :) = [];
    end
   
    % find nearest gene from open target for non-coding genes from VEP --------
    [vep.nearestGene, vep.gnomAD_variant] = deal(strings(height(vep), 1));
    for j = 1:numel(bgi.snp)
        vep_idx = find(vep.(1) == bgi.snp(j));

        if isfield(opts, "gtf"+refTag) %    use GTF file
            idx1 = gtf.seqname == string(bgi.chr(j));
            tmpgtf = gtf(idx1, :);
            tmpdi = double(bgi.pos(j));
            di = min(abs(tmpgtf.start - tmpdi), abs(tmpgtf.stop - tmpdi));
            [~, idx2] = min(di);
            ens = tmpgtf(idx2, :);
            if ~isempty(ens), vep.nearestGene(vep_idx) = ens.gene_name(1); end
        elseif opts.refGenome == "GRCh37"
            ens = EnsemblREST(bgi.chr(j) + ":" + bgi.pos(j), ...
                "nearest", "refGenome", "37");
            if ~isempty(ens), vep.nearestGene(vep_idx) = ens.name(1); end

        else % for GRCh38 use opentargets (doesn't support GRCh37)
            % step 1: use VEP table input IDs
            t = opentargets(vep.id(vep_idx));
            if ~isempty(t.data.search.variants)
                t = struct2table(t.data.search.variants, 'AsArray', true);
                if double(bgi.pos(j)) == t.positionB37(1)
                    nearestGene = {t.nearestCodingGene.symbol};
                    vep.nearestGene(vep_idx) = string(nearestGene(1));
                    continue
                end
            end
            
            % step 2: use VEP table annotated IDs
            t = opentargets(vep.snp(vep_idx));
            if ~isempty(t.data) && ~isempty(t.data.search.variants)
                t = struct2table(t.data.search.variants, 'AsArray', true);
                if double(bgi.pos(j)) == t.positionB37(1)
                    nearestGene = {t.nearestCodingGene.symbol};
                    vep.nearestGene(vep_idx) = string(nearestGene(1));
                    continue
                end
            end
    
            % step 3: use gnomAD 
            gnmd = gnomad('variantId', join([bgi.chr(j), ...
                string(bgi.pos(j)), bgi.a1(j), bgi.a2(j)], "-"), ...
                'dataset', opts.gdataset, 'reference_genome', opts.refGenome);
            if isempty(gnmd.data.variant)
                gnmd = gnomad('variantId', join([bgi.chr(j), ...
                    string(bgi.pos(j)), bgi.a2(j), bgi.a1(j)], "-"), ...
                    'dataset', opts.gdataset, 'reference_genome', opts.refGenome);
                if isempty(gnmd.data.variant)
                    continue
                end
            end
            
            t = opentargets(string(gnmd.data.variant.rsid));
            if ~isempty(t.data.search.variants)
                t = struct2table(t.data.search.variants, 'AsArray', true);
                if double(bgi.pos(j)) == t.positionB37(1)
                    vep.gnomAD_variant(vep_idx) = string(gnmd.data.variant.rsid);
                    nearestGene = {t.nearestCodingGene.symbol};
                    vep.nearestGene(vep_idx) = string(nearestGene(1));
                    continue
                end
            end
        end
    
    end

    fill_idx = vep.nearestGene == "";
    vep.nearestGene(fill_idx) = vep.gene(fill_idx);
    bgi = keepbgi;

else % all VEP genes are of protein_coding type
    vep.nearestGene = vep.gene;
    vep.gnomAD_variant = strings(height(vep), 1);
end

if opts.eqtl % get eQTL data for each variant
    vep.eqtl = cell(height(vep), 1);
    for i = 1:height(vep)
        if opts.verbose
            fprintf('fetching eQTL for variant %d of %d\n', i, height(vep))
        end
        eqtl = EnsemblREST(vep.id(i), 'eqtl', 'refGenome', erase(opts.refGenome, "GRCh"));
        if isempty(eqtl)
            eqtl = LDexpress(vep.id(i), 'verbose', false);
        end
        if isempty(eqtl)
            snp = split(vep.snp(i), ';'); snp(~startsWith(snp, 'rs')) = [];
            for j = 1:numel(snp)
                eqtl = EnsemblREST(snp(j), 'eqtl', 'refGenome', erase(opts.refGenome, "GRCh"));
                if isempty(eqtl)
                    eqtl = LDexpress(snp(j), 'verbose', false);
                end
                if ~isempty(eqtl)
                    vep.eqtl{i} = eqtl;
                    break
                end
            end

            if isempty(eqtl)
                if vep.gnomAD_variant(i) ~= ""
                    eqtl = EnsemblREST(vep.gnomAD_variant(i), ...
                        'eqtl', 'refGenome', erase(opts.refGenome, "GRCh"));
                    if isempty(eqtl)
                        eqtl = LDexpress(vep.gnomAD_variant(i), 'verbose', false);
                    end

                    if ~isempty(eqtl)
                        vep.eqtl{i} = eqtl;
                        continue
                    end
                end
            else
                continue
            end
        else
            vep.eqtl{i} = eqtl;
            continue
        end
    end
end

vep.hgvsc = arrayfun( @(x) join(split(x, ';'), newline), vep.hgvsc, 'uni', false);
vep.hgvsp = arrayfun( @(x) join(split(x, ';'), newline), vep.hgvsp, 'uni', false);

end