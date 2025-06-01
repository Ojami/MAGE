function res = pops(opts)
% @06APR2024: a wrapper for Polygenic Priority Score (PoPS)
% https://www.nature.com/articles/s41588-023-01443-6
% https://github.com/FinucaneLab/gene_features/tree/master
% https://github.com/FinucaneLab/pops

arguments
    opts.magma {mustBeFolder} % path to MAGMA files: genes.raw and genes.out must be there
    opts.popswd {mustBeFolder} = fileparts(which("pops.m"))
    opts.featurewd {mustBeFolder} = fullfile(fileparts(which("pops.m")), "features") % downloaded from https://www.finucanelab.org/data along with control features
    opts.out {mustBeTextScalar} = fullfile(pwd, "out")
    opts.outprefix {mustBeTextScalar} = "pops_results"
    opts.mungeonly (1,1) logical = false % only do the munging of features and return
    opts.chunk (1,1) double = 1
end

if ~isfolder(opts.out)
    mkdir(opts.out)
end

if ~isfolder(fullfile(opts.out, "features_munged"))
    mkdir(fullfile(opts.out, "features_munged"))
end

% check if already munged features are there
npy = getfilenames(fullfile(opts.out, "features_munged"), "npy").npy;
if ~isempty(npy)
    opts.ignoremunge = true;
    opts.chunk = numel(npy);
else
    opts.ignoremunge = false;
end


fi = ["featurewd", "popswd", "out", "magma"];
for k = 1:numel(fi)
    opts.(fi(k) + "wsl") = makeWSLpath(opts.(fi(k)));
    if ~endsWith(opts.(fi(k) + "wsl"), "/")
        opts.(fi(k) + "wsl") = opts.(fi(k) + "wsl") + "/";
    end
end

% create gene annotation file
g = struct2table(load("Homo_sapiens.GRCh37.87.gtf.gene.mat"));
g(g.gene_biotype ~= "protein_coding", :) = [];
g.TSS = g.start;
idx = g.strand == "-";
g.TSS(idx) = g.stop(idx);
g = g(:, ["gene_id", "gene_name", "seqname", "start", "stop", "TSS"]);
g.Properties.VariableNames = ["ENSGID", "NAME", "CHR", "START", "END", "TSS"];
mag = readtable(fullfile(opts.magma, "magma.genes.out"), "TextType", "string", "FileType", "text");
g(~ismember(g.ENSGID, mag.GENE), :) = [];
writetable(g, fullfile(opts.out, "gene_annot.txt"), "Delimiter", "\t")
dos2unix(fullfile(opts.out, "gene_annot.txt"))


% munge
if ~opts.ignoremunge
    cmd = "python3 " + opts.popswdwsl + "munge_feature_directory.py " + ...
            "--gene_annot_path " + opts.outwsl + "gene_annot.txt" + ...
            " --feature_dir " + opts.featurewdwsl + ...
            " --save_prefix " + opts.outwsl + "features_munged/pops_features --max_cols 5000 --nan_policy mean";
    runbash(cmd, "pops_bash", "parallel", false, "verbose", true);
end

if opts.mungeonly
    res = table;
    return
end

cmd = "python3 " + opts.popswdwsl + "pops.py " + ...
    "--gene_annot_path " + opts.outwsl + "gene_annot.txt " + ...
    "--feature_mat_prefix " + opts.outwsl + "features_munged/pops_features " + ...
    "--num_feature_chunks " + opts.chunk + ...
    " --magma_prefix " + opts.magmawsl + "magma --verbose " + ...
    "--out_prefix " + opts.outwsl + opts.outprefix  + ...
    " --control_features_path " + opts.popswdwsl + "control.features.txt ";
runbash(cmd, "pops_bash", "parallel", false, "verbose", true);


res = readtable(fullfile(opts.out, opts.outprefix + ".preds"), ...
    "TextType", "string", "VariableNamingRule", "preserve", "FileType", "text");

end