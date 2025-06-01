function joinFileCols(file, cols, opts)
% joinFileCols merges the file columns, 'cols', by joining them via 'join'
% option (default is ":") and creates a new file of joined values. This
% function's mainly use is to create unique GWAS ids based on
% chr:pos:sort(a1,a2), where a1 and a2 are sorted naturally (inspired by
% LDSC munge function). columns to be sorted should be added to 'sort'
% option.
% 
% To compile the c code:  gcc -o joinFileCols joinFileCols.c -lz

% Oveis Jamialahmadi, Univeristy of Gothenburg. November 2024.

arguments
    file {mustBeFile}
    cols {mustBeText, mustBeVector} % column names in 'file' to be joined: order matters
    opts.sep {mustBeTextScalar} % columns delimiter: determined automatically by bfilereader if left empty (recommended)
    opts.verbose (1,1) logical = true
    opts.gzip {mustBeMember(opts.gzip, ["gzip", "pigz"])} = "pigz"
    opts.sort {mustBeText, mustBeVector} % must be a subset of 'cols'
    opts.suffix {mustBeTextScalar} = "joined"
    opts.wrapper {mustBeMember(opts.wrapper, ["c", "awk"])} = "c" % default is C which is faster

end

file = string(file);

[~, file] = fileattrib(file);
file = string(file.Name);

opts.cols = cols; % keep col names
if isfield(opts, "sort")
    assert(all(ismember(opts.sort, cols)), "'sort' must be a subset of 'cols'!")
end

% check if cols exist in header
[cc, sep] = bfilereader(file, summary="firstline");
if ~isfield(opts, "sep"), opts.sep = sep.sep; end
assert(isempty(setdiff(cols, cc)), "Error:input column cannot be found in the header!")

[~, col_index]= ismember(cols, cc);

filein = "'" + makeWSLpath(file) + "'";
cmd = string;

% if wrapper is c: use the binary c code
if opts.wrapper == "c"
    exe = string(regexprep(which("joinFileCols.m"), ".m$", ""));
    assert(isfile(exe), "cannot find the binary file joinFileCols!")
    exe = "'" + makeWSLpath(exe) + "'";
    cmd = exe + " " + filein + " " + opts.suffix + ' "' + ...
        join(string(col_index), ",") + '"';

    if isfield(opts, "sort")
        [~, sort_col_index ]= ismember(opts.sort, cc);
        cmd = cmd + ' "' + join(string(sort_col_index), ",") + '"';
    end
    runbash(cmd, parallel=false, wait=true, verbose=opts.verbose);

    return
end

cmd(numel(cmd) + 1) = "# Input arguments";
cmd(numel(cmd) + 1) = "input_file=" + filein;
cmd(numel(cmd) + 1) = "suffix='" + opts.suffix + "'";

% # Columns to merge (1-based indices, comma-separated)
cmd(numel(cmd) + 1) = "cols=""" + join(string(col_index), ",") + """";

% # Columns to be naturally sorted (subset of cols, 1-based indices, comma-separated)
if isfield(opts, "sort")
    [~, sort_col_index ]= ismember(opts.sort, cc);
    cmd(numel(cmd) + 1) = "sort_cols=""" + join(string(sort_col_index), ",") + """";
else
    % nothing to sort
    cmd(numel(cmd) + 1) = "sort_cols=""""";
end


% # Detect if the input file is gzipped
cmd(numel(cmd) + 1) = 'if [[ $(head -c 2 "$input_file" | xxd -p) == "1f8b" ]]; then'; % for BGZF files
cmd(numel(cmd) + 1) = "  is_gzip=true";
cmd(numel(cmd) + 1) = "else";
cmd(numel(cmd) + 1) = "  is_gzip=false";
cmd(numel(cmd) + 1) = "fi";

% # Define output file name
cmd(numel(cmd) + 1) = 'base_name="${input_file%.*}"';
cmd(numel(cmd) + 1) = 'extension="${input_file##*.}"';
cmd(numel(cmd) + 1) = 'output_file="${base_name}_${suffix}.${extension}"';

% # Process the file
cmd(numel(cmd) + 1) = "process_file() {";
cmd(numel(cmd) + 1) = '  local infile="$1"';
cmd(numel(cmd) + 1) = '  local outfile="$2"';
cmd(numel(cmd) + 1) = "";
cmd(numel(cmd) + 1) = "  awk -F'[ \t]+' -v cols=""$cols"" -v sort_cols=""$sort_cols"" '";
cmd(numel(cmd) + 1) = "  BEGIN {";
cmd(numel(cmd) + 1) = "      # Parse column indices";
cmd(numel(cmd) + 1) = '      n = split(cols, merge_indices, ",")';
cmd(numel(cmd) + 1) = '      m = (sort_cols != "" ? split(sort_cols, sort_indices, ",") : 0)';
cmd(numel(cmd) + 1) = "";
cmd(numel(cmd) + 1) = "      # Prepare sort lookup for fast access";
cmd(numel(cmd) + 1) = "      if (m > 0) {";
cmd(numel(cmd) + 1) = "          for (j = 1; j <= m; j++) {";
cmd(numel(cmd) + 1) = "              sort_lookup[sort_indices[j]] = 1";
cmd(numel(cmd) + 1) = "          }";
cmd(numel(cmd) + 1) = "      }";
cmd(numel(cmd) + 1) = "  }";
cmd(numel(cmd) + 1) = "  {";
cmd(numel(cmd) + 1) = "      # Initialize arrays for sorted and non-sorted columns";
cmd(numel(cmd) + 1) = '      sorted_values = ""';
cmd(numel(cmd) + 1) = '      non_sorted_values = ""';
cmd(numel(cmd) + 1) = "";
cmd(numel(cmd) + 1) = "      # Process each column";
cmd(numel(cmd) + 1) = "      for (i = 1; i <= n; i++) {";
cmd(numel(cmd) + 1) = "          col_idx = merge_indices[i]";
cmd(numel(cmd) + 1) = "";
cmd(numel(cmd) + 1) = "          if (m > 0 && col_idx in sort_lookup) {";
cmd(numel(cmd) + 1) = '              sorted_values = (sorted_values == "" ? $col_idx : sorted_values ":" $col_idx)';
cmd(numel(cmd) + 1) = "          } else {";
cmd(numel(cmd) + 1) = '              non_sorted_values = (non_sorted_values == "" ? $col_idx : non_sorted_values ":" $col_idx)';
cmd(numel(cmd) + 1) = "          }";
cmd(numel(cmd) + 1) = "      }";

cmd(numel(cmd) + 1) = "      # Perform natural sorting on sorted_values";
cmd(numel(cmd) + 1) = '      if (m > 0 && sorted_values != "") {';
cmd(numel(cmd) + 1) = '          split(sorted_values, values, ":")';
cmd(numel(cmd) + 1) = "          asort(values)";
cmd(numel(cmd) + 1) = "          sorted_values = values[1]";
cmd(numel(cmd) + 1) = "          for (k = 2; k <= length(values); k++) {";
cmd(numel(cmd) + 1) = '              sorted_values = sorted_values ":" values[k]';
cmd(numel(cmd) + 1) = "          }";
cmd(numel(cmd) + 1) = "      }";
cmd(numel(cmd) + 1) = "";
cmd(numel(cmd) + 1) = "      # Combine non-sorted and sorted columns";
cmd(numel(cmd) + 1) = '      final = (non_sorted_values == "" ? sorted_values : (non_sorted_values ":" sorted_values))';
cmd(numel(cmd) + 1) = '      sub(/:$/, "", final)  # Remove trailing colon';
cmd(numel(cmd) + 1) = "      print final";
cmd(numel(cmd) + 1) = "  }";
cmd(numel(cmd) + 1) = "  ' ""$infile"" > ""$outfile""";
cmd(numel(cmd) + 1) = "}";

cmd(numel(cmd) + 1) = "# Handle gzipped or regular file";
cmd(numel(cmd) + 1) = "if $is_gzip; then";
cmd(numel(cmd) + 1) = "  temp_file=$(mktemp)";
cmd(numel(cmd) + 1) = "  " + opts.gzip + ' -cd "$input_file" > "$temp_file"';
cmd(numel(cmd) + 1) = '  process_file "$temp_file" "$output_file"';
cmd(numel(cmd) + 1) = '  rm -f "$temp_file"';
cmd(numel(cmd) + 1) = "  " + opts.gzip + ' "$output_file"';
cmd(numel(cmd) + 1) = "else";
cmd(numel(cmd) + 1) = '  process_file "$input_file" "$output_file"';
cmd(numel(cmd) + 1) = "fi";

% echo "Merged file saved to: $output_file"
cmd = cmd';

runbash(cmd, "joinFileCols", "wait", true, "parallel", false, ...
    "verbose", opts.verbose);


end % END

