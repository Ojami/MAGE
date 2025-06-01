function filterFileRows(file, val, col, opts)
% deleteRowsFromFile remove rows which the values in 'val' match values in
% file's column 'col'. The filtered file will be created with 

% Oveis Jamialahmadi, Univeristy of Gothenburg. November 2024.

arguments
    file {mustBeFile}
    val {mustBeText} % values to be matched: either string vector or a file name
    col {mustBeTextScalar} % column name in 'file'
    opts.sep {mustBeTextScalar} % columns delimiter: determined automatically by bfilereader if left empty (recommended)
    opts.infile (1,1) logical = false % to modify the original file (default). If false, writes to a new file with a "_filtered" suffix
    opts.verbose (1,1) logical = false
    opts.gzip {mustBeMember(opts.gzip, ["gzip", "pigz"])} = "pigz"
    opts.keep (1,1) logical = false % false: removes matching rows, true: keeps (only) matching rows
    
    %@25APR2025: comment was added to skip the comment lines
    opts.comment {mustBeTextScalar} = "#@!~&^%#$&()^%^$%#$@#%&^&"
end

tmp_name = getRandomName("tmp.txt", 5);

if isfile(val)
    copyfile(val, tmp_name);
else
    writematrix(val, tmp_name, QuoteStrings="none", FileType="text");
end

file = string(file);

[~, file] = fileattrib(file);
file = string(file.Name);

% check if cols exist in header
[cc, sep] = bfilereader(file, summary="firstline", comment=opts.comment);

if ~isfield(opts, "sep"), opts.sep = sep.sep; end

if all(~ismissing(double(col))) % col index 1-based
    col_index = col;
else
    assert(all(ismember(col, cc)), "Error:input column cannot be found in the header!")
    col_index = find(cc == col);
end

filein = "'" + makeWSLpath(file) + "'";

% cmd = string;
% cmd(1) = "dos2unix " + tmp_name;
% cmd(2) = "if file " + filein + " | grep " + '-q "gzip compressed"; then';
% cmd(3) = "    gzip -cd " + filein + " | awk -F'" + opts.sep + ...
%     "' -v col=" + col_index + " -v OFS='" + opts.sep + ...
%     "' 'NR==FNR {remove[$1]; next} !(($col) in remove)' " + ...
%     tmp_name + " - | gzip > filtered_output.gz";
% cmd(4) = "else";
% cmd(5) = "    cat " + filein + " | awk -F'" + opts.sep + ...
%     "' -v col=" + col_index + " -v OFS='" + opts.sep + ...
%     "' 'NR==FNR {remove[$1]; next} !(($col) in remove)' " + ...
%     tmp_name + " - > filtered_output.txt";
% cmd(6) = "fi";
% cmd = cmd';
% 
% runbash(cmd, "deleteRowsFromFile", "wait", true, "parallel", false);

cmd = string;
cmd(1) = "dos2unix " + tmp_name;
cmd(numel(cmd) + 1) = "# Input arguments";
cmd(numel(cmd) + 1) = "input_file=" + filein;
cmd(numel(cmd) + 1) = "pattern_file=" + tmp_name; 
cmd(numel(cmd) + 1) = "infile=" + string(opts.infile);
cmd(numel(cmd) + 1) = "keep=" + string(opts.keep);
cmd(numel(cmd) + 1) = "comment='" + string(opts.comment) + "'";
cmd(numel(cmd) + 1) = "col=" + col_index;

cmd(numel(cmd) + 1) = "# helper: first non-comment line = header";
cmd(numel(cmd) + 1) = "get_header() {";
cmd(numel(cmd) + 1) = '   grep -v "^${comment}" | head -n 1';
cmd(numel(cmd) + 1) = "}";

cmd(numel(cmd) + 1) = "# helper: strip comment lines";
cmd(numel(cmd) + 1) = "strip_comments() {";
cmd(numel(cmd) + 1) = '    grep -v "^${comment}"';
cmd(numel(cmd) + 1) = "}";

% cmd(numel(cmd) + 1) = 'if file "$input_file" | grep -q "gzip compressed"; then';
cmd(numel(cmd) + 1) = 'if [[ $(head -c 2 "$input_file" | xxd -p) == "1f8b" ]]; then'; % for BGZF files
cmd(numel(cmd) + 1) = "  is_gzip=true";
cmd(numel(cmd) + 1) = "else";
cmd(numel(cmd) + 1) = "  is_gzip=false";
cmd(numel(cmd) + 1) = "fi";

% # Process the file
cmd(numel(cmd) + 1) = "if $is_gzip; then";
  % # For gzip file
cmd(numel(cmd) + 1) = '  header=$(pigz -cd "$input_file" | get_header)';
cmd(numel(cmd) + 1) = "";
cmd(numel(cmd) + 1) = "  " + opts.gzip + ' -cd "$input_file" \';
cmd(numel(cmd) + 1) = "      | strip_comments \";
cmd(numel(cmd) + 1) = "      | tail -n +2 \";
cmd(numel(cmd) + 1) = "      | awk -F'" + opts.sep + "' -v OFS='" + opts.sep + "' -v col=""$col"" -v keep=""$keep"" \";
cmd(numel(cmd) + 1) = "          'NR==FNR {remove[$1]; next}";
cmd(numel(cmd) + 1) = '           (keep=="true"  && (($col) in remove)) ||';
cmd(numel(cmd) + 1) = "           (keep==""false"" && !(($col) in remove))' \";
cmd(numel(cmd) + 1) = "          ""$pattern_file"" - > temp_filtered.txt";


% cmd(numel(cmd) + 1) = "  " + opts.gzip + ' -cd "$input_file" | ' + "awk -F'" + opts.sep + ...
%     "' -v OFS='" + opts.sep + "' -v col=" + col_index + ...
%     " -v keep=""$keep"" 'NR==FNR {remove[$1]; next} " + ...
%     "(keep == ""true"" && (($col) in remove)) || " + ...
%     "(keep == ""false"" && !(($col) in remove))' " + ...
%     """$pattern_file"" - > temp_filtered.txt";
if opts.keep
     cmd(numel(cmd) + 1) = '    { printf "%s\n" "$header"; cat temp_filtered.txt; } \';
     cmd(numel(cmd) + 1) = "      > temp_filtered_hd.txt && mv temp_filtered_hd.txt temp_filtered.txt";
    % cmd(numel(cmd) + 1) = "{ " + opts.gzip + " -cd ""$input_file"" | head -n 1; cat temp_filtered.txt; } > temp_filtered_hd.txt && mv temp_filtered_hd.txt temp_filtered.txt";
end

cmd(numel(cmd) + 1) = "else";

% # For plain text file
cmd(numel(cmd) + 1) = '    header=$(strip_comments < "$input_file" | head -n 1)';
cmd(numel(cmd) + 1) = '    strip_comments < "$input_file" \';
cmd(numel(cmd) + 1) = "      | tail -n +2 \";
cmd(numel(cmd) + 1) = "      | awk -F'" + opts.sep + "' -v OFS='" + opts.sep + "' -v col=""$col"" -v keep=""$keep"" \";
cmd(numel(cmd) + 1) = "          'NR==FNR {remove[$1]; next}";
cmd(numel(cmd) + 1) = '           (keep=="true"  && (($col) in remove)) ||';
cmd(numel(cmd) + 1) = "           (keep==""false"" && !(($col) in remove))' \";
cmd(numel(cmd) + 1) = "          ""$pattern_file"" - > temp_filtered.txt";


% cmd(numel(cmd) + 1) = "awk -F'" + opts.sep + ...
%     "' -v col=" + col_index + " -v OFS='" + opts.sep + ...
%     "' -v keep=""$keep"" 'NR==FNR {remove[$1]; next} " + ...
%     "(keep == ""true"" && (($col) in remove)) || " + ...
%     "(keep == ""false"" && !(($col) in remove))' " + ...
%     """$pattern_file"" ""$input_file"" > temp_filtered.txt";

if opts.keep
    cmd(numel(cmd) + 1) = '    { printf "%s\n" "$header"; cat temp_filtered.txt; } \';
    cmd(numel(cmd) + 1) = "      > temp_filtered_hd.txt && mv temp_filtered_hd.txt temp_filtered.txt";
    % cmd(numel(cmd) + 1) = "{ head -n 1 ""$input_file""; cat temp_filtered.txt; } > temp_filtered_hd.txt && mv temp_filtered_hd.txt temp_filtered.txt";
end
cmd(numel(cmd) + 1) = "fi";

% # Handle output based on 'infile' flag
cmd(numel(cmd) + 1) = 'if [[ "$infile" == "true" ]]; then';
  % # Replace original file
cmd(numel(cmd) + 1) = "  if $is_gzip; then";
cmd(numel(cmd) + 1) = "    " + opts.gzip + ' < temp_filtered.txt > "$input_file"';
cmd(numel(cmd) + 1) = "  else";
cmd(numel(cmd) + 1) = '    mv temp_filtered.txt "$input_file"';
cmd(numel(cmd) + 1) = "  fi";
cmd(numel(cmd) + 1) = "  echo ""Original file '$input_file' has been replaced.""";
cmd(numel(cmd) + 1) = "else";
  % # Create a new file with "_filtered" suffix
cmd(numel(cmd) + 1) = '  output_file="${input_file%.*}_filtered.${input_file##*.}"';
cmd(numel(cmd) + 1) = "  if $is_gzip; then";
cmd(numel(cmd) + 1) = "    " + opts.gzip + ' < temp_filtered.txt > "$output_file"';
cmd(numel(cmd) + 1) = "  else";
cmd(numel(cmd) + 1) = '    mv temp_filtered.txt "$output_file"';
cmd(numel(cmd) + 1) = "  fi";
cmd(numel(cmd) + 1) = "  echo ""Filtered file created: '$output_file'""";
cmd(numel(cmd) + 1) = "fi";

% # Clean up
cmd(numel(cmd) + 1) = "rm -f temp_filtered.txt";
cmd = cmd';

runbash(cmd, "deleteRowsFromFile", "wait", true, "parallel", false, ...
    "verbose", opts.verbose, "delete", true);
delete(tmp_name)

end % END

%% future development: numerical filtering based on multiple cols
% #!/bin/bash
% # Input arguments
% input_file='/path/to/your/input_file.txt'
% pattern_file='patterns.txt'   # File with strings for text matching
% text_filter=true              # Enable text filtering
% numerical_filter=false        # Enable numerical filtering (set to true when needed)
% 
% # Numerical filtering variables
% declare -A conditions
% conditions[1]="> 0.5"         # Format: column_number="operator threshold"
% conditions[2]="< 100"         # Add more columns as needed
% logical_operator="AND"        # "AND" for intersection, "OR" for union of criteria
% 
% # Check if input is compressed
% if [[ $(head -c 2 "$input_file" | xxd -p) == "1f8b" ]]; then
%   is_gzip=true
% else
%   is_gzip=false
% fi
% 
% # Function for text filtering
% text_filtering() {
%   awk -F'	' -v OFS='	' -v col=4 -v keep=true '
%   NR==FNR {patterns[$1]; next}
%   (keep == "true" && ($col in patterns)) || (keep == "false" && !($col in patterns)) {
%     print $0
%   }' "$pattern_file" "$1" > temp_filtered.txt
% }
% 
% # Function for numerical filtering
% numerical_filtering() {
%   awk -F'	' -v OFS='	' -v logical_operator="$logical_operator" -v header=1 '
%   BEGIN {
%     # Parse conditions
%     num_conditions = ARGV[1]
%     delete ARGV[1]
%     for (i = 1; i <= num_conditions; i++) {
%       split(ARGV[i], cond, " ");
%       col[i] = cond[1];
%       op[i] = cond[2];
%       val[i] = cond[3];
%       delete ARGV[i];
%     }
%   }
%   NR == 1 {print; next} # Always include the header
%   {
%     match_count = 0;
%     for (i = 1; i <= num_conditions; i++) {
%       value = $col[i];
%       if ((op[i] == ">" && value > val[i]) ||
%           (op[i] == "<" && value < val[i]) ||
%           (op[i] == "==" && value == val[i]) ||
%           (op[i] == ">=" && value >= val[i]) ||
%           (op[i] == "<=" && value <= val[i]) ||
%           (op[i] == "!=" && value != val[i])) {
%         match_count++;
%       }
%     }
%     # Apply logical operator
%     if ((logical_operator == "AND" && match_count == num_conditions) ||
%         (logical_operator == "OR" && match_count > 0)) {
%       print $0;
%     }
%   }' "${#conditions[@]}" "${!conditions[@]}" "$1" > temp_filtered.txt
% }
% 
% # Determine which filtering to apply
% if $text_filter && ! $numerical_filter; then
%   if $is_gzip; then
%     pigz -cd "$input_file" | text_filtering
%   else
%     text_filtering "$input_file"
%   fi
% elif $numerical_filter && ! $text_filter; then
%   if $is_gzip; then
%     pigz -cd "$input_file" | numerical_filtering
%   else
%     numerical_filtering "$input_file"
%   fi
% else
%   echo "Error: Enable either text_filter=true or numerical_filter=true, but not both."
%   exit 1
% fi
% 
% # Save output
% output_file="${input_file%.*}_filtered.${input_file##*.}"
% if $is_gzip; then
%   pigz < temp_filtered.txt > "$output_file"
% else
%   mv temp_filtered.txt "$output_file"
% fi
% echo "Filtered file created: '$output_file'"
% 
% # Clean up
% rm -f temp_filtered.txt
