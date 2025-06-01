function geno = readvcf(region, opts)

% uses bcftools to read variants from a VCF file. If not indexed, the
% function first indexes the input vcf file.

% Oveis Jamialahmadi, 30 April 2024. University of Gothenburg.

arguments 
    region {mustBeText, mustBeVector} % region(s) should be fetched in format of chrN:POS1-POS2
    opts.vcf {mustBeFile, mustBeVector} % input vcf file if not provided looks in 'home' directory to find the pattern
    opts.home {mustBeFolder} % directory where VCF/index files exist
    opts.parallel (1,1) logical = true % run jobs in parallel
    opts.dtype {mustBeMember(opts.dtype, ["single", "double"])} = "double"
    opts.fields {mustBeMember(opts.fields, ["GT", "DS", "HDS", "GP"])} = ["GT", "DS", "HDS", "GP"] % fields to read
end

geno = struct;
assert(any(isfield(opts, ["vcf", "home"])), "readvcf::at least one of 'vcf' or 'home' should be provided!")

% find vcf file if not provided
if ~isfield(opts, "vcf")
    files = getfilenames(opts.home, "vcf.gz").vcf_gz;
    if isempty(files)
        files = getfilenames(opts.home, "vcf").vcf;
    end
    assert(~isempty(files), "readvcf::there is no VCF file in 'home' directory!")
    files = fullfile(opts.home, files);
    
    patt = replace(strshare(files, "pattern", true), "%d", "%s");
    for k = 1:numel(region)
        chr = regexprep(extractBefore(region(k), ":"), "^chr", "", "ignorecase");
        opts.vcf(k, 1) = compose(patt.replace("\", "/"), chr);
    end
end

opts.out = getRandomName("vcfreader");
opts.outwsl = makeWSLpath(fullfile(pwd, opts.out));

% check if index files are present, if not first index
vcf_files = unique(opts.vcf);
cmd = strings(numel(vcf_files), 1);
for k = 1:numel(vcf_files)

    if ~any(isfile(vcf_files(k) + [".tbi", ".csi"]))
        name = makeWSLpath(vcf_files(k).replace(["\", "/"], filesep));
        cmd(k) = bcftools + " index " + name + " -o " + name + ".csi";
    end

end

cmd(cmd == "") = [];
if ~isempty(cmd)
    runbash(cmd, "readvcf", "parallel", opts.parallel, "verbose", true);
end


cmd = strings(numel(region), 1);
for k = 1:numel(region)
    cmd(k) = bcftools + " view " + makeWSLpath(opts.vcf(k).replace(["\", "/"], filesep)) + ' "' + region(k) + '" > "' + opts.outwsl + k + '.txt"';
end

% get sample names
cmd(numel(cmd) + 1) = bcftools + " query -l " + makeWSLpath(opts.vcf(k).replace(["\", "/"], filesep)) + ' >"' + opts.outwsl + '.sample.txt"';
runbash(cmd, "readvcf", "parallel", opts.parallel, "verbose", true);
geno.eid = readmatrix(opts.out + ".sample.txt", "FileType", "text", "NumHeaderLines", 0, "OutputType", "string");
try
    ds = tabularTextDatastore(opts.out + (1:numel(region)) + ".txt", "TextType", ...
        "string", "CommentStyle", "##", "VariableNamingRule", "preserve", ...
        "ReadVariableNames", true, "NumHeaderLines", 0);
    ds = readall(ds);
catch
    ds = cell(numel(region), 1);
    for k = 1:numel(ds)
        ds{k} = readtable(opts.out + k + ".txt", "CommentStyle", "##", ...
            "ReadVariableNames", true, "VariableNamingRule", "preserve", ...
            "NumHeaderLines", 0);
    end
    ds(cellfun(@isempty, ds)) = [];
    ds = vertcat(ds{:});
end

arrayfun(@delete, opts.out + (1:numel(region)) + ".txt")
delete(opts.out + ".sample.txt")

if isempty(ds), return; end

% append variant info to meta struct
cols = colnames(ds);
mcols.name = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"];
mcols.tag = ["chr", "pos", "snp", "ref", "alt", "QUAL", "FILTER"];
for k = 1:numel(mcols.name)
    if any(cols == mcols.name(k))
        geno.(mcols.tag(k)) = ds.(mcols.name(k));
        ds.(mcols.name(k)) = [];
    end
end

% extract imputation score
c = cols(cols.lower == "info");
geno.R2 = nan(height(ds), 1);
if ~isempty(c)
    idx = contains(ds.(c), "R2=");
    geno.R2(idx) = double(extractBetween(ds.(c)(idx), (";"|textBoundary("start")) + "R2=", (";"|textBoundary("end"))));
    ds.(c) = [];
end

% find format
if ~any(colnames(ds) == "FORMAT"), return; end
fr = ds.FORMAT(1).split(":");

% GT - Estimated most likely genotype.
% DS - Estimated alternate allele dosage [P(0/1)+2*P(1/1)].
% HDS - Estimated phased haploid alternate allele dosage.
% GP - Estimated Posterior Genotype Probabilities P(0/0), P(0/1) and P(1/1).
dd = table(["GT", "HDS", "GP", "DS"]', [2, 2, 3, 1]');
dd.Properties.VariableNames = ["tag", "span"];
[~, idx] = ismember(fr, dd.tag); idx(idx < 1) = [];
dd = dd(idx, :);
dd.span = [cumsum(dd.span)-dd.span+1, cumsum(dd.span)];
dd(~ismember(dd.tag, opts.fields), :) = [];
fr = dd.tag;

ds.FORMAT = [];
ds = ds(:, geno.eid);
ds = rows2vars(ds);
ds = ds{:, 2:end};
ds = double(split(ds, [":", "|", ",", "/"]));
if opts.dtype ~= "double"
    ds = feval(opts.dtype, ds);
end
ct = 1;
while true
    if size(ds, 3) > 1
        bed = squeeze(ds(:, 1, :));
        ds(:, 1, :) = [];
    else
        bed = ds; ds = [];
    end
    
    for k = 1:height(dd)
        geno.(fr(k)){ct, 1} = bed(:, dd.span(k, 1):dd.span(k, 2));
    end
    clear bed

    ct = ct + 1;
    if isempty(ds), break; end

end

for k = 1:numel(fr)
    geno.(fr(k)) = horzcat(geno.(fr(k)){:});
end
% 
% geno.desc = table(["GT", "DS", "HDS", "GP"]', ...
%     ["Estimated most likely genotype", "Estimated alternate allele dosage [P(0/1)+2*P(1/1)]", ...
%     "Estimated phased haploid alternate allele dosage", "Estimated Posterior Genotype Probabilities P(0/0), P(0/1) and P(1/1)"]');

end % END