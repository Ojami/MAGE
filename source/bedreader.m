function [getBED, getBIM, getFAM] = bedreader(fileName, inputSNPs, opts)
% getBED: a nXm matrix, which n corresponds to number of individuals (in
% getFAM) and m corresponds to variants (m = 2* variant numbers). It
% returns data as:
%     - 0,0 for homozygote major 
%     - 0,1 for heterozygote
%     - 1,1 for homozygote minor
%     - -2,-2 for missing genotype

% inputSNPs = "rs11591147"; % only if snp-major, otherwise must be sample ID 
% @06/01/2022: some bugs were fixed.
% Oveis Jamialahmadi, University of Gothenburg, Jan 2019.

arguments
    fileName {mustBeFile}
    inputSNPs {mustBeVector, mustBeText}
    opts.bim (1,1) logical = true
    opts.fam (1,1) logical = true
    opts.parallel (1,1) logical = false
end

fclose('all');

if numel(string(fileName)) > 1
    error('bedreader only accepts one bed file!')
else
    fileName = string(fileName);
end

if ~isstring(inputSNPs)
    inputSNPs = string(inputSNPs);
end

fileName = regexprep(fileName, '.bed$', '');
bedfile = fileName + ".bed";
bimfile = fileName + ".bim";
famfile = fileName + ".fam";

% sanity check -------------------------
% check snp/sample-major file
magicNumber = ["01101100"; "00011011"]; % BED magic number
fid = fopen(char(bedfile), 'r');
bedType = fread(fid, 3);
fclose(fid);
if bedType(3) == 0
    mode = "individual";
elseif bedType(3) == 1
    mode = "snp";
else
    error("unrecognized mode flag!")
end

bedType = join(string(de2bi(bedType(1:2), 8, 'left-msb')),'');
if ~all(magicNumber == bedType)
    error('input file is not a BED file!')
end

% read BIM file ------------------------
tic
getBIM = bimreader(bimfile);

if mode == "snp"
    inSnpid = find(ismember(getBIM.snp, inputSNPs));
    if isempty(inSnpid)
        error('bedreader failed to find input rsIDs in provided bim file!')
    end
end

% Read BIM file
if opts.bim
    if mode ~= "snp"
        getBIM = getBIM.snp;
    else
        fi = fieldnames(getBIM);
        for i = 1:numel(fi)
            getBIM.(fi{i}) = getBIM.(fi{i})(inSnpid); 
        end
    end
end
toc
if opts.fam
    % Red FAM file
    getFAM = famreader(famfile);
    
    if mode == "individual"
        inSnpid = find(ismember(getFAM, double(string(inputSNPs)))); % assuming double (may not always be the case)
        if isempty(inSnpid)
            error('bedreader failed to find input sample in provided fam file!')
        end
        getFAM = getFAM(inSnpid);
    end
end


nstart0 = 4; % Starting byte for BED file

if mode == "individual"
    snplength = ceil(numel(getBIM.snp)/4);
elseif mode == "snp"
    snplength = ceil(numel(getFAM)/4); % number of bytes per SNP
end

% Every byte encodes up to four genotypes (2 bits per genotype).
getBED = bedchunkMex(char(bedfile), nstart0, snplength, inSnpid, opts.parallel);
if mode == "individual"
    getBED(numel(getBIM.snp)+1:end,:) = [];
elseif mode == "snp"
    getBED(numel(getFAM)+1:end,:) = [];
end
fclose('all');
end

%%  subfunctions ==========================================================
% function getBIM = bimreader(fileName)
% try
%     getBIM = tabularTextDatastore(fileName,'FileExtensions', '.bim', ...
%         'TextType', 'string', 'NumHeaderLines', 0);
% catch % only 1 variant is present
%     getBIM = tabularTextDatastore(fileName,'FileExtensions', '.bim','TextType',...
%         'string','Delimiter', char(9), 'ReadVariableNames', false, ...
%         'NumHeaderLines', 0);
% end
% 
% getBIM.SelectedVariableNames(3) = [];
% getBIM = readall(getBIM);
% getBIM.Properties.VariableNames = {'chr', 'snp', 'pos', 'a1', 'a2'}; % a1: ALT, a2: REF in PLINK 2
% getBIM = table2struct(getBIM, 'ToScalar', true);
% 
% % getBIM = bfilereader(fileName, 'extractCol', 2, 'verbose', 'off', 'parallel',true);
% % getBIM.Properties.VariableNames{1} = 'snp';
% end

%% 
function getFAM = famreader(fileName)
fid = fopen(fileName, 'r');
line1 = fgetl(fid);
frewind(fid);
fmt = [repmat('%*s', 1, numel(split(string(line1)))), '%*[^\n]'];
fmt(2) = [];
getFAM = textscan(fid, fmt, 'CollectOutput', true);
getFAM = double(string(getFAM{1}));
fclose(fid);
end