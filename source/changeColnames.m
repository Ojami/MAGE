function changeColnames(files, oldcols, newcols, opts)
% changeColnames rename column in the 'cols' string vector in the header of
% input file in 'files'.

% Oveis Jamialahmadi, Univeristy of Gothenburg. August 2024.

arguments
    files {mustBeFile}
    oldcols {mustBeText, mustBeVector} % columns to be renamed
    newcols {mustBeText, mustBeVector} % new column names
    opts.sep {mustBeTextScalar} = "\t" % columns delimiter
    opts.infile (1,1) logical = true % to modify the original file (default). If false, writes to a new file with a ".newcols" suffix
    opts.suffix {mustBeTextScalar} = "newcols" % only if infile is false
    opts.gzip {mustBeMember(opts.gzip, ["gzip", "pigz"])} = "pigz"

end

assert(numel(oldcols) == numel(newcols), "New and old columns must have the same size!")
opts.suffix = "." + regexprep(opts.suffix, '^\.+', '');

files = string(files);

for k = 1:numel(files)
    file = files(k);
    [~, file] = fileattrib(file);
    [file, files(k)] = deal(string(file.Name));
    
    % check if cols exist in header
    [cc, sep] = bfilereader(file, summary="firstline");
    sep = sep.sep;
    assert(all(ismember(oldcols, cc)), "Error:some columns cannot be found in the header!")
    
    % check that newcols differ from other cols (no duplicates)
    cc_check = setdiff(cc, oldcols);
    assert(~all(ismember(newcols, cc_check)), "duplicated columns found! choose new column names!")

    % construct the awk command 
    [f1, f2] = ismember(cc, oldcols);
    cc(f1) = newcols(f2(f1));
    mcmd = join(cc, sep);

    filein = "'" + makeWSLpath(file) + "'";

    if file.endsWith(".gz")
        cmd = opts.gzip + " -dc " + filein + " | ";
        file = regexprep(file, ".gz$", "");
        gzipped = true;
        opts.suffix2 = opts.suffix;
    else
        gzipped = false;
        cmd = "";
        if opts.infile
            opts.suffix2 = "";
        else
            opts.suffix2 = opts.suffix;
        end
    end

    if opts.infile && ~gzipped
        iopt = " -i ";
        [~, ~, ext] = fileparts(file);
        fileout = regexprep(file, ext + "$", "") + opts.suffix2 + ext;
    else
        iopt = "";
        [~, file, ext] = fileparts(file);
        fileout = file + opts.suffix2 + ext;
    end
    fileout = "'" + makeWSLpath(fileout);

    
    if gzipped
        cmd = cmd + "sed " + iopt + "'1s/.*/" + mcmd + "/' " + " | " + opts.gzip + " > " + fileout + ".gz'";
    else
        if opts.infile
            cmd = "sed " + iopt + "'1s/.*/" + mcmd + "/' " + filein;
        else
            cmd = "sed " + iopt + "'1s/.*/" + mcmd + "/' " + filein + " > " + fileout + "'";
        end
    end
    runbash(cmd, "changeColnames", "wait", true);

    if gzipped
        fileout = fileout + ".gz";
        fileout = fileout.erase("'");
        if ispc, fileout = makeWSLpath(fileout, true); end
        if opts.infile
            movefile(fileout, files(k))
        end
    end

end

end % END