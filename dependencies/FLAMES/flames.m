function flames(opts)
% A wrapper for FLAMES gene prioritization pipeline:
% https://github.com/Marijn-Schipper/FLAMES

arguments
    opts.flames {mustBeFile} = fullfile(fileparts(which("flames.m")), "flames.py") 
    opts.a {mustBeFolder} % = fullfile(fileparts(which("flames.m")), "Annotation_data") 
    opts.o {mustBeTextScalar} = pwd() 
    opts.p {mustBeFile} % /path/to/PoPS preds file
    opts.m {mustBeFile} % MAGMA *.genes.out file
    opts.mt {mustBeFile} % MAGMA *.gsa.out file
    opts.id {mustBeFile} % indexfile of credible sets (see makeFinemaptable function)
    opts.conda {mustBeTextScalar} = "flames"
end

if ~isfolder(opts.o), mkdir(opts.o); end


% python FLAMES.py annotate \
% -o {DESIRED_OUTPUT_DIRECTORY} \
% -a {PATH_TO_THE_DOWNLOADED_ANNOTATION_DATA_DIRECTORY} \
% -p {DESIRED_POPS_OUTPUT_PREFIX}.preds \
% -m {DESIRED_ZSCORE_FILENAME}.genes.out \
% -mt {DESIRED_TISSUE_RELEVANCE_FILENAME}.gsa.out \
% -id {PATH_TO_INDEXFILE} 
% 
% check if annotations already exist
apth = readtable(opts.id, TextType="string", Delimiter="\t");
apth = apth.Annotfiles;
apth = arrayfun(@(x)makeWSLpath(x, true), apth);

run_flames_ann = ~all(isfile(apth));

if run_flames_ann
    cmd = "python3 " + makeWSLpath(opts.flames) + " annotate" + ...
        " -a " + makeWSLpath(opts.a) + ...
        " -o " + makeWSLpath(opts.o) + ...
        " -p " + makeWSLpath(opts.p) + ...
        " -m " + makeWSLpath(opts.m) + ...
        " -mt " + makeWSLpath(opts.mt) + ...
        " -id " + makeWSLpath(opts.id);
    
    runbash(cmd, "flames_bash", parallel=false, ...
        verbose=true, conda=opts.conda);
end

% run FLAMES scroing step
% python FLAMES.py FLAMES \
% -id {INDEX_FILE_NCLUDING_COLUMN Annotfiles} \
% -o {DESIRED_OUTPUT_DIRECTORY}

if isfile(fullfile(opts.o, "FLAMES_scores.raw"))
    fprintf("Skipped FLAMES, already exists: %s\n", ...
        fullfile(opts.o, "FLAMES_scores.raw"))
    return
end
cmd = "python3 " + makeWSLpath(opts.flames) + " FLAMES " + ...
    " -o " + makeWSLpath(opts.o) + ...
    " -id " + makeWSLpath(opts.id);
runbash(cmd, "flames_bash", parallel=false, verbose=true, conda=opts.conda);

end % END