function res = pops(opts)
% @06APR2024: a wrapper for Polygenic Priority Score (PoPS)
% https://www.nature.com/articles/s41588-023-01443-6
% https://github.com/FinucaneLab/gene_features/tree/master
% https://github.com/FinucaneLab/pops
% 
% @05SEP2025: used munged_features from FLAMES pipeline (full set)
%             -gene_annot.txt from FLAMES does not cover all genes,
%             therefore we generate a temp one based on ENSEMBL GRCh37.87
arguments
    opts.magma {mustBeFolder} % path to MAGMA files: genes.raw and genes.out must be there
    opts.popswd {mustBeFolder} = fileparts(which("pops.m"))
    opts.featurewd {mustBeFolder} = fullfile(fileparts(which("pops.m")), "pops_features_full_FUMA_compatible") % downloaded from FLAMES
    opts.out_path {mustBeTextScalar} = fullfile(pwd, "out")
    opts.outprefix {mustBeTextScalar} = "pops_results"
    opts.num_feature_chunks (1,1) double = 116
    % opts.gene_annot_path {mustBeFile} = fullfile(fileparts(which("pops.m")), "pops_features_full", "gene_annot.txt")
    % opts.feature_mat_prefix {mustBeTextScalar} = fullfile(fileparts(which("pops.m")), "pops_features_full", "munged_features", "pops_features")
    % opts.control_features {mustBeFile} = fullfile(fileparts(which("pops.m")), "pops_features_full", "control.features")
end

if ~isfolder(opts.out_path)
    mkdir(opts.out_path)
end

opts.out_path_wsl = makeWSLpath(opts.out_path);
if ~endsWith(opts.out_path_wsl, "/")
    opts.out_path_wsl = opts.out_path_wsl + "/";
end

% % create gene annotation file
% g = struct2table(load("Homo_sapiens.GRCh37.87.gtf.gene.mat"));
% g(g.gene_biotype ~= "protein_coding", :) = [];
% g.TSS = g.start;
% idx = g.strand == "-";
% g.TSS(idx) = g.stop(idx);
% g = g(:, ["gene_id", "gene_name", "seqname", "start", "stop", "TSS"]);
% g.Properties.VariableNames = ["ENSGID", "NAME", "CHR", "START", "END", "TSS"];
% mag = readtable(fullfile(opts.magma, "magma.genes.out"), "TextType", "string", "FileType", "text");
% g(~ismember(g.ENSGID, mag.GENE), :) = [];
% writetable(g, fullfile(opts.featurewd, "gene_annot.txt"), "Delimiter", "\t")
% dos2unix(fullfile(opts.featurewd, "gene_annot.txt"))

ch = ["featurewd", "magma"];
for k = 1:numel(ch)
    opts.(ch(k)) = makeWSLpath(opts.(ch(k)));
    if ~endsWith(opts.(ch(k)), "/")
        opts.(ch(k)) = opts.(ch(k)) + "/";
    end
end

% " --gene_annot_path " + makeWSLpath(fullfile(opts.popswd, "gene_annot.txt")) + ...
cmd = "python3 " + makeWSLpath(fullfile(opts.popswd, "pops.py")) + ...
    " --gene_annot_path " + opts.featurewd + "gene_annots.txt" + ...
    " --feature_mat_prefix " + opts.featurewd + "features_munged/pops_features " + ...
    "--num_feature_chunks " + opts.num_feature_chunks + ...
    " --magma_prefix " + opts.magma + "magma --verbose " + ...
    "--out_prefix " + opts.out_path_wsl + opts.outprefix  + ...
    " --control_features_path " + opts.featurewd + "control.features";
runbash(cmd, "pops_bash", "parallel", false, "verbose", true, "conda", "pops");

% delete(fullfile(opts.popswd, "gene_annot.txt"))

res = readtable(fullfile(opts.out_path, opts.outprefix + ".preds"), ...
    "TextType", "string", "VariableNamingRule", "preserve", "FileType", "text");



end