function call_geneWrapper

% an example of how to call geneWrapper function. To use this example,
% unzip the files in misc.7z
% 12JULY2024

% schemes are different annotation cateogries using different algorithms.
% scheme 3 and 4 are the same, with the difference that qced variants were
% not filtered (i.e. those from DNAnexus) in scheme 4.

% @16JULY2024: based on some benchmarks on different algorithms, it turned
% out that CADD >= 20 has the highest loading/coeff in PCA on whole exome
% matrix composing of all selected algorithms compared to a more
% conservative 30 cut-off. Also, PrimateAI doesn't seem to add much to
% neither of top PCs. So, this is the latest scheme: scheme 5
% 
% @07APR2025: "vep.missense" now is used over "vep.general" to restrict the
% down-stream analyses to missnese+lof variants only (for the sake of
% consistency with other studies/reports).

clc
hm = fileparts(which("call_geneWrapper.m"));
qcpath = "X:\WES\WES.500K\plink.qc"; % QC files are in this directory
genes = load("protein_coding_genes_grch38.mat").genes; % set of genes to keep
qc = getfilenames(fullfile(hm, "wes500.variant.qc"), "txt", fullpath=true).txt; % provided by DNANexus
wes = fullfile(hm, "annotation", "WES.500K.txt"); % see runVEP

% define categories
categories = struct;
method = "revel(0.5),cadd(20),sift,polyphen(s),LRT," + ...
    "MutationTaster,M-CAP,AlphaMissense";
categories.name = ["consensus"; "relaxed"; "consensus4"; "AlphaMissense"];
categories.method = [repmat(method, numel(categories.name)-1, 1); "AlphaMissense"];
categories.rule = ["and"; "or"; "4"; "and"];

method = "revel(0.5),cadd(20)";
categories.name = [categories.name; "REVEL5ORCAD20"];
categories.method = [categories.method; method];
categories.rule = [categories.rule; "or"];
% 
% method = "revel(0.5),cadd(30),M-CAP,AlphaMissense";
% categories.name = [categories.name; "RCMA3"];
% categories.method = [categories.method; method];
% categories.rule = [categories.rule; "3"];

categories = struct2table(categories);
geneWrapper(qcpath, wes, geneset=genes, exclude=qc, categories=categories, ...
    model="vep.missense")
% geneWrapper(freq, wes, geneset=genes, categories=categories)

end % END