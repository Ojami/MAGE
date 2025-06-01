function out = vcftrimmer(infiles, opts)
% vcftrimmer trims the input VCF (variant call format) and removes genotype
% calls to prepare it for variant annotation software (e.g. VEP or
% SnpEff). It can also converts a BIM/BGI file to a VCF file suitable for
% VEP/SnpEff.
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, University of Gothenburg, March 2022.

arguments
    infiles {mustBeFile, mustBeVector}
    opts.parallel (1,1) logical = false
    opts.workers (1,1) double = 10
    opts.ref_first (1,1) logical = false % if true, treats the first allele column in BIM (or allele1 in BGI) as REF, and the other allele as ALT. In case of UKBB PLINK, should be false
end

% check inputs
if ~isstring(infiles); infiles = string(infiles); end
if ~all(endsWith(infiles, [".bim", ".bgi", ".vcf"]))
    error('vcftrimmer:invalidDataType', 'only bim/big/vcf files are allowed!')
end

% check file types
[~, ~, ext] = fileparts(infiles);
opts.bgi = ext == ".bgi";
opts.bim = ext == ".bim";
opts.vcf = ext == ".vcf";
if any(opts.vcf)
    % vcftrimmerParallel is outdated and nonsensical.
    fprintf('vcftrimmer: using vcf files is deprecated, use bim/bgi files\n')
    if ~opts.parallel % cannot perform without parallel toolbox
        pans = input("'parallel' flag must be true with vcf files!\n" + ...
            "[1]ignore vcf files\n[2]use parallel toolbox (default:2):");
        if pans == 1
            infiles(opts.vcf) = [];
            opts.bgi(opts.vcf) = [];
            opts.bim(opts.vcf) = [];
            opts.vcf = false(numel(infiles), 1);
        else
            opts.parallel = true;
        end
    end
end

if opts.parallel && isempty(gcp('nocreate'))
    parc = parcluster('local');
    parc.NumWorkers = opts.workers;
    parpool(parc);
end

out = strings(numel(infiles), 1);
for i = 1:numel(infiles)
    fprintf("(%d of %d)-creating vcf for: %s\n", i, numel(infiles), infiles(i))
    if opts.vcf(i)
        vcftrimmerParallel(infiles(i))
    elseif opts.bim(i)
        out(i, 1) = bim2vcf(infiles(i), opts);
    elseif opts.bgi(i)
        out(i, 1) = bgi2vcf(infiles(i), opts);
    end
end


end

%% subfunctions -----------------------------------------------------------
function out = bgi2vcf(fileName, opts)
conn = sqlite(fileName);

tic
getBGI = fetch(conn, 'SELECT * FROM Variant');
close(conn)
getBGI(:, [4, 7, 8]) = [];
if istable(getBGI), getBGI = getBGI.Variables; end

if ~opts.ref_first % allele1 is ALT
    tmp = getBGI(:, 4);
    getBGI(:, 4) = getBGI(:, 5);
    getBGI(:, 5) = tmp;
    clear tmp
end

% Remove starting 0 from chromosome
getBGI(:, 1) = regexprep(getBGI(:, 1), '^0', '');

N = size(getBGI, 1);
getBGI = [getBGI, repmat(".",N,2),  repmat(["PR","GT"],N,1)];
getBGI = join(getBGI, char(9));
toc

vcfDate = string(datetime("now", "Format", "MM/dd/uuuu"));
vcfHdr = ["##fileformat=VCFv4.2";"##fileDate=" + vcfDate;
"##source=MATLAB vcfTrimmer bim2vcf output";
"##contig=<ID=BGI,length=000000>";
"##INFO=<ID=PR,Number=0,Type=Flag,Description=""Provisional reference allele, may not be based on real reference genome"">";
"##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">"
"#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT"];


% fprintf('Writing data to a new VCF file...\n')
[pth, out] = fileparts(fileName);
out = fullfile(pth, out + ".t.vcf");
fileID = fopen(out,'w');
fprintf(fileID,'%s\n',vcfHdr(:));
fprintf(fileID,'%s\n',getBGI(1:end-1));
fprintf(fileID,'%s',getBGI(end));
fclose(fileID);

fclose('all');

end

%% ------------------------------------------------------------------------
function out = bim2vcf(fileName, opts)

try
    getbim = tabularTextDatastore(fileName,'FileExtensions', '.bim', ...
        'TextType', 'string', 'NumHeaderLines', 0);
catch % only 1 variant is present
    getbim = tabularTextDatastore(fileName,'FileExtensions', '.bim','TextType',...
        'string','Delimiter', char(9), 'ReadVariableNames', false, ...
        'NumHeaderLines', 0);
end
getbim.SelectedVariableNames(3) = [];
getbim = readall(getbim);

% at least in UKBB, PLINK BIM files are formatted as chr, id, cm, pos, alt,
% ref
if opts.ref_first % last column is ALT
    getbim = [getbim.(1), getbim.(3), getbim.(2), getbim.(4), getbim.(5),...
        repmat(".",size(getbim,1),2),repmat(["PR","GT"],size(getbim,1),1)];
else
    getbim = [getbim.(1), getbim.(3), getbim.(2), getbim.(5), getbim.(4),...
        repmat(".",size(getbim,1),2),repmat(["PR","GT"],size(getbim,1),1)];
end
getbim = join(getbim,char(9));

vcfDate = string(datetime("now", "Format", "MM/dd/uuuu"));
vcfHdr = ["##fileformat=VCFv4.2";"##fileDate=" + vcfDate;
"##source=MATLAB vcfTrimmer largebim2vcf output";
"##contig=<ID=allchr,length=000000>";
"##INFO=<ID=PR,Number=0,Type=Flag,Description=""Provisional reference allele, may not be based on real reference genome"">";
"##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">"
"#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT"];

[pth, out] = fileparts(fileName);
out = fullfile(pth, out + ".t.vcf");
fileID = fopen(out,'w');
fprintf(fileID,'%s\n',vcfHdr(:));
fprintf(fileID,'%s\n',getbim(1:end-1));
fprintf(fileID,'%s',getbim(end));
fclose(fileID);
fclose('all');

end

%% ------------------------------------------------------------------------
function vcftrimmerParallel(fileName)
% DEPRECATED
fid = fopen(fileName, 'rt');
fprintf('Finding headers....\n')
goME = true;
while goME
    tline = textscan(fid, '%c%*[^\n]', 1);
    if char(tline) ~= '#'
        fseek(fid,tellMeLine,'bof');
        goME = false;
    else
        tellMeLine = ftell(fid)+1;
    end
end
fclose(fid);

m = memmapfile(fileName);
headerStr = split(string(char(m.Data(1:tellMeLine)')),newline);
clear m
headerStr = strrep(headerStr,char(13),'');

headerStrMain = headerStr(end);
if headerStrMain == ""
    headerStrMain = headerStr(end-1);
end
headerStrMain = split(headerStrMain,char(9));
headerStrMain = headerStrMain(cellfun(@isempty,regexp(headerStrMain(1:20),'\d+_\d+')));
% headerStrMain(end) = [];

headerStr_new = headerStr;
headerStr_new(3) = headerStr_new(3) + " MATLAB vcfTrimmer output";
headerStr_new(end) = join(headerStrMain);

skipLines = numel(headerStr) - 1;

fprintf('Getting datastore....\n')
tic
warning off
ds = datastore(fileName,'FileExtensions','.vcf','Type','tabulartext','NumHeaderLines',skipLines);
getFmt = numel(headerStrMain); %count spaces
warning on
toc

% fmt = [ repmat('%s',1,getFmt), repmat('%*s',1,sum(uint8(headerStr{end}) == 9) - getFmt), '%*[^\n]' ];
% fmt = [ repmat('%s',1,getFmt), '%*[^\n]' ];
fprintf('Running mapreduce....\n')
tic
mapRes = mapreduce(ds, @(data,info,kvs)vcftrimmerMapper(data,info,kvs,getFmt), @vcftrimmerReducer);
mapRes = readall(mapRes);
toc

fprintf('Getting data....\n')
[~,offsetIDXsorted] = sort(mapRes.Key);
mapRes = mapRes.Value;
mapRes = mapRes(offsetIDXsorted);
mapRes = vertcat(mapRes{:});
mapRes = join(mapRes,char(9));

%Delete unnecessary MAT files
matFiles = dir('*.mat');
matFiles = {matFiles.name};
matFiles = matFiles(~cellfun(@isempty,regexp(matFiles,'^result_')));
delete(matFiles{:})

fprintf('Writing data to a new VCF file...\n')
[pth, outfileName, fileNameTAG] = fileparts(fileName);
outfileName = [outfileName,'.t',fileNameTAG];
fileID = fopen(fullfile(pth, outfileName),'w');
fprintf(fileID,'%s\n',headerStr_new(:));
fprintf(fileID,'%s\n',mapRes(1:end-1));
fprintf(fileID,'%s',mapRes(end));
fclose(fileID);
fprintf('VCF file was successfully written to %s\n',outfileName)
fclose('all');
end

%--------------------------------------------------------------------------
function vcftrimmerMapper(data, intermKeys, intermKVStore, getFmt)
 data = string(data.(data.Properties.VariableNames{1}));
 trimmedSTR = arrayfun(@(x)textscan(x,[repmat('%s',1,getFmt),'%*[^\n]']),data,'UniformOutput',false);
 trimmedSTR = cellfun(@(x)string(x),trimmedSTR,'UniformOutput',false);
 trimmedSTR = vertcat(trimmedSTR{:});
 offsetData = intermKeys.Offset;
%  
%  dataTrimmed.offset = offsetData;
%  dataTrimmed.str = trimmedSTR;
% add(intermKVStore, "dataTrimmed", dataTrimmed)

add(intermKVStore, offsetData, trimmedSTR)
end

function vcftrimmerReducer(intermKey, intermValsIter, outKVStore)
  trimmedSTR = []; 
%   offsetData = [];
  while(hasnext(intermValsIter))
%       dataTrimmed = getnext(intermValsIter);
%       offsetData = [offsetData; dataTrimmed.offset];
%       trimmedSTR = [trimmedSTR; dataTrimmed.str];
      trimmedSTR = [trimmedSTR; getnext(intermValsIter)];
  end
  
  add(outKVStore,intermKey,trimmedSTR);
%   add(outKVStore, offsetData, trimmedSTR);
% 
%   add(outKVStore, "trimmedSTR", trimmedSTR);
%   add(outKVStore, "offsetData", offsetData);
end
%--------------------------------------------------------------------------