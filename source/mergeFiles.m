function mergeFiles(files, opts)
% an updated version of 'mergeStatFiles' function. 
% Oveis Jamialahmadi, University of Gothenburg, August 2024.

arguments
    files {mustBeFile}
    opts.output {mustBeTextScalar} = "merged_file.txt" % output file path/name: fullfile(pwd, "merged.txt")
    opts.cols {mustBeText, mustBeVector} % columns to merge. If empty, uses all columns
    opts.delimiter {mustBeTextScalar} = "\t" % delimiter of output file
    opts.sort (1,1) logical = true % natsort files?
end

% check output
opts.output = string(opts.output);
[pth, opts.output, ext] = fileparts(opts.output);
opts.output = opts.output + ext;
if pth == ""
    opts.output = fullfile(pwd, opts.output);
else
    if ~isfolder(pth), mkdir(pth); end
    [~, pth] = fileattrib(pth);
    opts.output = fullfile(pth.Name, opts.output);
end
opts.output = makeWSLpath(opts.output);

% find header
hdr = bfilereader(files(1), "summary", "firstline");

% check input files
if opts.sort, files = natsort(files); end
for k = 1:numel(files)
    file = files(k);
    [~, file] = fileattrib(file);
    file = file.Name;
    files(k) = makeWSLpath(file);
end

cmd = string;
cmd(1) = "# Output file name";
cmd(numel(cmd) + 1) = 'output_file="' + opts.output + '"';

if isfield(opts, "cols")
    idx = find(ismember(hdr, opts.cols));
else
    idx = 1:numel(hdr);
end
cmd(numel(cmd) + 2) = "# Columns to include (1-based, comma-separated)";
cmd(numel(cmd) + 1) = 'columns="' + join(string(idx), ",") + '"';

cmd(numel(cmd) + 1) = 'output_delimiter="' + opts.delimiter + '"';

cmd(numel(cmd) + 1) = 'files=("' + files.join('" "') + '")';

cmd(numel(cmd) + 1) = "# Convert columns to awk's format";
cmd(numel(cmd) + 1) = "awk_columns=$(echo $columns | awk -F, '{for (i=1;i<=NF;i++) printf ""%s%s"", $i, (i<NF ? OFS : """")}')";

cmd(numel(cmd) + 1) = "# Process files with awk";
cmd(numel(cmd) + 1) = "awk -v cols=""$awk_columns"" -v OFS=""$output_delimiter"" '";
cmd(numel(cmd) + 1) = "    BEGIN { split(cols, c, "" "") } ";
cmd(numel(cmd) + 1) = "   FNR == 1 && NR != 1 { next }  # Skip header of all but first file";
cmd(numel(cmd) + 1) = "    {";
cmd(numel(cmd) + 1) = '       for (i in c) printf "%s%s", (i>1?OFS:""), $c[i]';
cmd(numel(cmd) + 1) = '       print ""';
cmd(numel(cmd) + 1) = "    }";
cmd(numel(cmd) + 1) = "' ""${files[@]}"" > ""$output_file""";

runbash(cmd, "mergeFiles", parallel=false, wait=true, verbose=true);


end %END
