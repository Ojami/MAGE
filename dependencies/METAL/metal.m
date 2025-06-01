function metal(files, opts)
%@05AUG2024: a wrapper for METAL meta-analysis package.
% https://genome.sph.umich.edu/wiki/METAL_Documentation

arguments
    files {mustBeFile}
    opts.metal {mustBeFile} = regexprep(which("metal.m"), ".m$", "") % metal executable

    % column of input files: should have the same number and in the orde as
    % 'files'. If is scalar, it assumes all input files have the same
    % column names.
    opts.snp {mustBeText, mustBeVector} % SNP column
    opts.ea {mustBeText, mustBeVector} % effect allele column
    opts.eaf {mustBeText, mustBeVector} % effect allele frequency column  
    opts.nea {mustBeText, mustBeVector} % other allele 
    opts.effect {mustBeText, mustBeVector} % effect column (beta)
    opts.p {mustBeText, mustBeVector}
    opts.n {mustBeText, mustBeVector} % WEIGHT 
    opts.se {mustBeText, mustBeVector}
    opts.chr {mustBeText, mustBeVector}
    opts.pos {mustBeText, mustBeVector}

    % other options
    opts.SEPARATOR {mustBeMember(opts.SEPARATOR, ["WHITESPACE", "COMMA", "TAB"])} = "WHITESPACE"    
    opts.HETEROGENEITY (1,1) logical = true % To allow for heterogeneity
    opts.SCHEME {mustBeMember(opts.SCHEME, ["STDERR", "SAMPLESIZE"])} = "STDERR"
    opts.GENOMICCONTROL {mustBeMember(opts.GENOMICCONTROL, ["OFF", "ON"])} = "ON"
    opts.OVERLAP {mustBeMember(opts.OVERLAP, ["OFF", "ON"])} = "OFF" % Sample Overlap Correction: valid only with SCHEME SAMPLESIZE
    opts.VERBOSE {mustBeMember(opts.VERBOSE, ["OFF", "ON"])} = "OFF"
    opts.AVERAGEFREQ {mustBeMember(opts.AVERAGEFREQ, ["OFF", "ON"])} = "ON"
    opts.MINMAXFREQ {mustBeMember(opts.MINMAXFREQ, ["OFF", "ON"])} = "ON"
    opts.TRACKPOSITIONS {mustBeMember(opts.TRACKPOSITIONS, ["OFF", "ON"])} = "ON"
    opts.OUTFILE {mustBeTextScalar} = "METAANALYSIS" % suffix is .tbl
    opts.LOGPVALUE {mustBeMember(opts.LOGPVALUE, ["OFF", "ON"])} = "OFF"
    opts.EFFECT_PRINT_PRECISION (1,1) double = 7
    opts.STDERR_PRINT_PRECISION (1,1) double = 7

    % filtering: can be a vector if more filtering criteria are needed
    opts.ADDFILTER {mustBeText, mustBeVector} % e.g. ADDFILTER=["N > 1000", "MAF > 0.01", "MARKER_ID IN (rs1234,rs123456,rs123)"]

end

if opts.SCHEME ~= "SAMPLESIZE"
    opts.OVERLAP = "OFF"; % only valid for this scheme
end

% check column options
mp = struct;
mp.n1 = ["snp", "eaf", "effect", "p", "n", "se", "chr", "pos"]';
mp.n2 = ["MARKER", "FREQLABEL", "EFFECT", "PVALUE", "WEIGHT", "STDERR", "CHROMOSOME", "POSITION"]';

% write the scripts
sr = string;
sr(1) = "#METAL wrapper";
if isfield(opts, "n")
    sr(2) = "CUSTOMVARIABLE N";
    sr(3) = "LABEL N as " + opts.n(1);
end
op = ["GENOMICCONTROL", "SEPARATOR", "SCHEME", "OVERLAP", "VERBOSE", ...
    "AVERAGEFREQ", "MINMAXFREQ", "TRACKPOSITIONS", "LOGPVALUE", ...
    "EFFECT_PRINT_PRECISION", "STDERR_PRINT_PRECISION"];
for k = 1:numel(op)
     sr(numel(sr) + 1) = op(k) + " " + opts.(op(k));
end

% check custrom filters
if isfield(opts, "ADDFILTER")
    for k = 1:numel(opts.ADDFILTER)
        sr(numel(sr) + 1) = "ADDFILTER " + opts.ADDFILTER(k);
    end
end

inopts = opts;
inopts = rmfield(inopts, setdiff(fieldnames(inopts), union(mp.n1, ["ea", "nea"])));
fi = string(fieldnames(inopts));

% check size of inputs
for k = 1:numel(fi)
    if numel(files) ~= numel(inopts.(fi(k))) % use the value for all inputs
        inopts.(fi(k)) = repmat(inopts.(fi(k)), numel(files), 1);
    end
end

for k = 1:numel(files)
    sr(numel(sr) + 2) = "# FILE " + k;
    
    visitedEA = false;
    for j = 1:numel(fi)
        if any(fi(j) == ["ea", "nea"])
            if ~visitedEA
                sr(numel(sr) + 1) = "ALLELE " + inopts.ea(k) + " " +inopts.nea(k);
                visitedEA = true;
            end
        else
            cc = mp.n2(mp.n1 == fi(j));
            sr(numel(sr) + 1) = cc + " " + inopts.(fi(j))(k);
        end
    end
    
    sr(numel(sr) + 1) = "PROCESS " + makeWSLpath(files(k));

end

sr(numel(sr) + 2) = "OUTFILE " + opts.OUTFILE + char(32) + ".txt";

if opts.HETEROGENEITY
    sr(numel(sr) + 1) = "ANALYZE HETEROGENEITY";
else
    sr(numel(sr) + 1) = "ANALYZE";
end

sr(numel(sr) + 1) = "QUIT";

writelines(sr', "METAL_SCRIPT.txt")
dos2unix(fullfile(pwd, "METAL_SCRIPT.txt"), "verbose", true)

cmd = makeWSLpath(opts.metal) + " " +  ...
    makeWSLpath(fullfile(pwd, "METAL_SCRIPT.txt"));
runbash(cmd, "META_SCRIPT", "verbose", true);
if isfile("METAL_SCRIPT.txt"), delete("METAL_SCRIPT.txt"); end
if isfile(opts.OUTFILE + "1.txt")
    movefile(opts.OUTFILE + "1.txt", opts.OUTFILE + ".txt")
    movefile(opts.OUTFILE + "1.txt.info", opts.OUTFILE + ".txt.info")
end

% To apply genomic control to the meta-analysis results, just perform an
% initial meta-analysis and then load the initial set of results into METAL
% to get final, genomic control adjusted results.
if opts.GENOMICCONTROL == "ON"
    % write the scripts
    sr = string;
    sr(1) = "#METAL wrapper";
    if isfield(opts, "n")
        sr(2) = "CUSTOMVARIABLE N";
        sr(3) = "LABEL N as " + opts.n(1);
    end
    
    for k = 1:numel(op)
         sr(numel(sr) + 1) = op(k) + " " + opts.(op(k));
    end
    
    mp = struct;				
    mp.n1 = ["MarkerName", "Freq1", "Effect", "P-value", "StdErr", "Chromosome", "Position"]';
    mp.n2 = ["MARKER", "FREQLABEL", "EFFECT", "PVALUE", "STDERR", "CHROMOSOME", "POSITION"]';
    for k = 1:numel(mp.n1)
        sr(numel(sr) + 1) = mp.n2(k) + " " + mp.n1(k);
    end
    sr(numel(sr) + 1) = "ALLELE Allele1 Allele2";
    sr(numel(sr) + 1) = "PROCESS " + makeWSLpath(opts.OUTFILE + ".txt");

    sr(numel(sr) + 2) = "OUTFILE " + opts.OUTFILE + " .txt";
    sr(numel(sr) + 1) = "ANALYZE";
    sr(numel(sr) + 1) = "QUIT";
    writelines(sr', "METAL_SCRIPT.txt")
    dos2unix(fullfile(pwd, "METAL_SCRIPT.txt"), "verbose", true)
    
    cmd = makeWSLpath(opts.metal) + " " +  ...
        makeWSLpath(fullfile(pwd, "METAL_SCRIPT.txt"));
    runbash(cmd, "METAL_SCRIPT", "verbose", true);

    if isfile(opts.OUTFILE + "1.txt")
        movefile(opts.OUTFILE + "1.txt", opts.OUTFILE + ".GC.txt")
        movefile(opts.OUTFILE + "1.txt.info", opts.OUTFILE + ".GC.txt.info")
    end

end
if isfile("METAL_SCRIPT.txt"), delete("METAL_SCRIPT.txt"); end



end % END