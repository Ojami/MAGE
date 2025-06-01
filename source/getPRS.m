function prs = getPRS(tab, opts)
% calculates a polygenic risk score for variants with their weights in
% 'tab'. 'tab' must have the following columns (or alternatively use
% options):
%   -chr: chromosome 
%   -snp: variant ids as they appear in the input bgi/bim files
%   -beta: weight of each variant in 'snp'
%   -ea: effect allele for beta
%   -nea: non-effect allele
% 
% Oveis Jamialahmadi, University of Gothenburg. August 2024.

arguments
    tab {mustBeA(tab, "table")} 
    opts.chr {mustBeTextScalar} = "chr"
    opts.snp {mustBeTextScalar} = "snp" 
    opts.beta {mustBeTextScalar} = "beta"
    opts.ea {mustBeTextScalar} = "ea"
    opts.nea {mustBeTextScalar} = "nea"

    % options for fetching imputed  genotypes
    opts.bgenhome {mustBeFolder}
    opts.parallel (1,1) logical = false

    opts.verbose (1,1) logical = true
end

prs = struct;
% verify columns ----------------------------------------------------------
cols = colnames(tab);
fi = ["chr", "snp", "beta", "ea", "nea"];
for k = 1:numel(fi)
    assert(any(cols == opts.(fi(k))), "column %s was not found!", opts.(fi(k)))
end

tab_raw = tab; % for output weights keeping also original table
tab = tab(:, [opts.chr, opts.snp, opts.beta, opts.ea, opts.nea]);
tab = renamevars(tab, 1:5, ["chr", "snp", "beta", "ea", "nea"]);
tab.chr = string(tab.chr);

% fetch variants 
geno = getbulkgeno(tab.snp, tab.chr, parallel=opts.parallel, ...
        datatype="double", verbose=false, merge=true, infoscore=1e-6);


% match geno with variants in tab -----------------------------------------
tab.id = tab.chr + ":" + tab.snp + ":" + tab.ea + ":" + tab.nea;
geno.idf = geno.chr + ":" + geno.snp + ":" + geno.a2 + ":" + geno.a1;
geno.idr = geno.chr + ":" + geno.snp + ":" + geno.a1 + ":" + geno.a2;

non_matched_idx = ~(ismember(geno.idf, tab.id) | ismember(geno.idr, tab.id));
fi = setdiff(fieldnames(geno), ["eid", "bed"]);
if any(non_matched_idx)
    if opts.verbose, fprintf("%d non-matching variants were removed from geno file.\n", sum(non_matched_idx)); end

    for k = 1:numel(fi)
        geno.(fi(k))(non_matched_idx) = [];
    end
    geno.bed(:, non_matched_idx) = [];
end

% realign geno calls 
idx = ~ismember(geno.idf, tab.id);
if any(idx)
    geno.bed(:, idx) = 2 - geno.bed(:, idx);
    a1 = geno.a2(idx);
    geno.a2(idx) = geno.a1(idx);
    geno.a1(idx) = a1;
    geno.idf = geno.chr + ":" + geno.snp + ":" + geno.a2 + ":" + geno.a1;
end

assert(all(ismember(geno.idf, tab.id)), "Something went wrong!")
[f1, f2] = ismember(tab.id, geno.idf); f2(f2<1) = [];

% calculate and return the PRS
prs.eid = geno.eid;
prs.value = geno.bed(:, f2)*tab.beta(f1);

tab = tab(f1, :);
tab_raw.id = tab_raw.(opts.chr) + ":" + tab_raw.(opts.snp) + ...
    ":" + tab_raw.(opts.ea) + ":" + tab_raw.(opts.nea);
[~, idx] = ismember(tab.id, tab_raw.id);
tab_raw = tab_raw(idx, :);
tab_raw.id = [];

prs.weights = tab_raw;
clear geno

end % END