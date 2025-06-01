function bim = bimreader(bimfile, opts)
% par of bedreader function for reading bim files, allowing the selection
% of 'loader' and 'columns' to be read.
% Oveis Jamialahmadi, University of Gothenburg, Apr 2022.
% 
% @19MAY2025: 'extract' option was added to extract only rows containing
%             variants in 'extract'.

arguments
    bimfile {mustBeFile}
    opts.cols {mustBeVector, mustBeMember(opts.cols, ["chr", "snp", "cmg", "pos", "a1", "a2"])} = ["chr", "snp", "pos", "a1", "a2"]
    opts.loader {mustBeTextScalar, mustBeMember(opts.loader, ["matlab", "bfilereader"])} = "matlab"; % "matlab" is usually faster
    opts.struct (1,1) logical = true
    opts.extract {mustBeVector}
end

% is it a bim file?
if ~endsWith(bimfile, ".bim")
    error('bimreader:wrongInput', 'not a bim file!')
end

% a1: ALT, a2: REF in PLINK 2
defcols = {'chr', 'snp', 'cmg', 'pos', 'a1', 'a2'};

if isfield(opts, "extract")
    filterFileRows(bimfile, opts.extract, "2", infile=false, keep=true);
    opts.extract = [];
    [pth, name, ext] = fileparts(bimfile);
    bimfile = fullfile(pth, name + "_filtered" + ext);
    opts.delbim = true; % remove this temp file when finished
else
    opts.delbim = false;
end

if opts.loader == "matlab"
    try
        bim = tabularTextDatastore(bimfile,'FileExtensions', '.bim', ...
            'TextType', 'string', 'NumHeaderLines', 0);
    catch % only 1 variant is present
        bim = tabularTextDatastore(bimfile,'FileExtensions', '.bim','TextType',...
            'string','Delimiter', char(9), 'ReadVariableNames', false, ...
            'NumHeaderLines', 0);
    end

    bim.VariableNames = defcols;
    bim.SelectedVariableNames = opts.cols;
    idx = find(ismember(bim.SelectedVariableNames, ["a1", "a2"]));
    bim.SelectedFormats(idx) = repmat({'%q'}, 1, numel(idx));

    if isempty(gcp("nocreate"))
        bim = readall(bim);
    else % parpool is active
        matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
        bim = gather(tall(bim));
        matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);
    end
    
else
    colidx = find(ismember(defcols, opts.cols));
    bim = bfilereader(bimfile, 'extractCol', colidx, 'verbose', 'off', 'return', 'rawTable');
    bim.Properties.VariableNames = defcols(colidx);
    bim = bim(:, opts.cols);
    bim.pos = double(bim.pos);
end

if opts.struct
    bim = table2struct(bim, 'ToScalar', true);
end

if opts.delbim
    delete(bimfile);
end

end % END