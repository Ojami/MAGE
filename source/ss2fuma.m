function ss2fuma(sumstat, opts)
% converts summary statistics file for FUMA web service. 
% The output of this function should be uploaded on https://fuma.ctglab.nl/
% 
% Oveis Jamialahmadi, University of Gothenburg, March 2023.

arguments
    sumstat {mustBeFile}
    opts.output {mustBeTextScalar}
    opts.leadsnps {mustBeA(opts.leadsnps, 'table')} % lead snp file
    opts.map {mustBeA(opts.map, 'table')} % mapping between old and new variable names with two columns: 1-old and 2-new
    opts.metal (1,1) logical = false % if true deletes some extra METAL specific columns, and writes a light sumstat file
    opts.threads (1,1) double = 35 % only when metal is true
end

t1 = tic;
tmpname = getRandomName("fuma", 5);
sumstat = string(sumstat);
[pth, name] = fileparts(sumstat);
if pth == "", pth = pwd; end

if isfield(opts, 'output')
    [pth, name, ext] = fileparts(string(opts.output));

    if pth == "", pth = pwd; end

    if ~isfolder(pth)
        mkdir(pth)
    end
    output = fullfile(pth, name + ext);

else
    output = fullfile(pth, name + "_fuma");
end

if isfield(opts, 'leadsnps')
    tab = opts.leadsnps;
    h = colnames(tab);
    snpidx = find(ismember(h.lower, ["snp", "varid", "rsid", "rsids", "id"]), 1);
    chridx = find(contains(h.lower, "chr"), 1);
    posidx = find(ismember(h.lower, ["bp", "position", "pos"]) | contains(h.lower, "genpos"), 1);
    tab = tab(:, [snpidx, chridx, posidx]);
    writetable(tab, fullfile(pth, name + "_leadsnps.txt"), Delimiter="\t");
end

hdr = bfilereader(sumstat, "header", true, "summary", "firstline");

if any(contains(hdr.lower, "log10p")) % regenie
    pcol = find(contains(hdr.lower, "log10p"));

     % convert log10 P -> P
    runbash("awk 'NR!=1 {$"+ pcol(1) + "=10^-$" + pcol(1) + ...
        "}1' " + makeWSLpath(sumstat) + " > " + makeWSLpath(output), ...
        tmpname, "wait", true);
    sumstat = output;
end

% change column names
if isfield(opts, 'map')
    hdr_new = hdr;
    for k = 1:height(opts.map)
        [idx1, idx2] = ismember(hdr.lower, lower(opts.map.old));
        hdr_new(idx1) = opts.map.new(idx2);
    end
    hdr = hdr_new;
else
    hdr = createGWASheader(hdr);
end

if opts.metal
    opts.hdr = hdr;
    parseMETAL2FUMA(sumstat, output, opts);
else
    if sumstat == output
        cmd = "sed -i '1s/.*/" + hdr.join(" ") + "/' '" + makeWSLpath(sumstat) + ...
            "'";
    else
        cmd = "sed '1s/.*/" + hdr.join(" ") + "/' '" + makeWSLpath(sumstat) + ...
            "'  > '" + makeWSLpath(output) + "'";
    end
    runbash(cmd, tmpname, "wait", true);
end

% compress it
gzip(output)
delete(output)

fprintf('summary stat file %s is ready (done in %.0f sec)\n', output, toc(t1))

end %END

%% subfunctions ===========================================================
function parseMETAL2FUMA(sumstat, output, opts)

if isfile(output)
    return
end

if isempty(gcp("nocreate"))
    parpool("Processes", opts.threads);
end

ds = tabularTextDatastore(sumstat, TextType="string", VariableNamingRule="preserve");
delim = ds.Delimiter;

ds = tall(ds);
keep_hdr = setdiff(opts.hdr, ["FreqSE", "MinFreq", "MaxFreq", "Direction", ...
        "HetISq", "HetChiSq", "HetDf", "HetPVal"], "stable");
ds.Properties.VariableNames = opts.hdr;

ds = ds(:, keep_hdr);
fastWriteTable(ds, output=output, delimiter=delim)
dos2unix( output, verbose=true);

end % END