function addIDtoGWASfile(files, cols, opts)
% addIDtoGWASfile adds a new column by merging columns in 'cols' to the
% file, and writes it to a new file with a suffix in 'suffix'. This is
% useful for meta-analysis and other purposes where unique ID column
% ('newcol') is required to avoid confusion (e.g. using SNP column with
% duplicated names).

% Oveis Jamialahmadi, Univeristy of Gothenburg. August 2024.

arguments
    files {mustBeFile}
    cols {mustBeText, mustBeVector} % columns to join
    opts.sep {mustBeTextScalar} = ":" % separator to merge 'cols'
    opts.delim {mustBeTextScalar} = "\t" % columns delimiter
    opts.newcol {mustBeTextScalar} = "MERGED_ID" % new column by mergeing 'cols'
    opts.suffix {mustBeTextScalar} = "merged" % suffix to add to new edited file: file.suffix
end


files = string(files);
for k = 1:numel(files)
    file = files(k);
    
    % check if cols exist in header
    cc = bfilereader(file, summary="firstline");
    assert(all(ismember(cols, cc)), "Error:some columns cannot be found in the header!")
    assert(all(~ismember(opts.newcol, cc)), "newcol cannot be an already existing column!")

    % construct the awk command 
    [~, idx] = ismember(cols, cc); idx(idx < 1) = [];
    mcmd = join("$" + idx, ' "' + opts.sep + '" ');
    
    filein = "'" + makeWSLpath(file) + "'";

    if file.endsWith(".gz")
        cmd = "zcat " + filein + " | ";
        file = regexprep(file, ".gz$", "");
        gzipped = true;
    else
        gzipped = false;
        cmd = "";
    end

    [~, ~, ext] = fileparts(file);
    fileout = regexprep(file, ext + "$", "") + "." + opts.suffix + ext;
    fileout = "'" + makeWSLpath(fileout);
    
    if gzipped
        cmd = cmd + "awk 'BEGIN {OFS=""" + opts.delim + """} {if (NR==1) " + ...
            "{print $0, """ + opts.newcol + """} else " + ...
            "{print $0, " + mcmd + "}}' | gzip > " + fileout + ".gz'";
    else
        cmd = "awk 'BEGIN {OFS=""" + opts.delim + """} {if (NR==1) " + ...
            "{print $0, """ + opts.newcol + """} else " + ...
            "{print $0, " + mcmd + "}}' " + filein + " > " + fileout + "'";
    end
    runbash(cmd, "addIDtoGWASfile", "wait", true);

    % normalize delimiter
    if ~gzipped
        cmd = "awk 'BEGIN {OFS=""" + opts.delim + """} " + ...
            "{gsub(/[ \t]+/, """ + opts.delim + """); print}' " + ...
            fileout + "' > " + fileout + "_TMP' && mv " + fileout + ...
            "_TMP' " + fileout + "'";

        % cmd = "sed -i 's/ \+/ \t/g' " + fileout + "'";
        runbash(cmd, "addIDtoGWASfile", "wait", true);
    end
end

end % END
