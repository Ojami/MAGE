function callData = bgenreader(bgen_file, opts)
% INPUTS:
%   - bgen_file: a character or string array for name (and path if not in
%       pwd) of the input bgen file (.bgen). The corresponding bgi file
%       with the same name also must accompany the bgen_file.
%       NOTE: if .sample file corresponding to the input bgen file cannot
%       be found, the function will search for all available sample files
%       in the current directory.
% OPTIONAL INPUTS:
%   - opts.varInfo: a string array containing rsid or position of
%       input query variants; can be mixture of both (id and position interval). For
%       instance: ["rs372654813", "107289436-107289437"]
%       Default: empty; extracting all variants (not recommended)
%   - opts.gpType: 
%         - 'dosage' (DEPRECATED)
%         - 'GT', hard genotype call (same as PLINK2
%           --import-dosage-certainty 0.9)
%         - 'GP', returns a matrix with 2 columns, corresponding to Aa and
%           aa respectively. (default)
%   - opts.verbose: printing info (enabled: true, disables: false). Default:
%              true.
%   - opts.numWorkers: Number of threats for parallel computing. Default: [].
%   - opts.save: a logical indicating bim/bed/fam to be saved (true).
%               Default: false
%   - opts.datatype: 'single' or 'double' opts.datatype for genotype calls.
% 
% OUTPUT:
%   -callData: a struct with fields: bim, bed, fam, I_A (imputation
%   information score).
% 
% Oveis Jamialahmadi, GU, Jan. 2020.
% 
% @02/17/2022: several improvments have been made. 
% @09AUG2024: now supports BGEN v1.3 with zstd compression.

arguments
    bgen_file {mustBeTextScalar, mustBeFile} % path to bgen file
    opts.varInfo {mustBeVector, mustBeText} 
    opts.gpType {mustBeTextScalar, mustBeMember(opts.gpType, ["GP", "GT"])} = "GP" % "dosage" is deprecated
    opts.verbose (1,1) logical = false
    opts.numWorkers (1,1) double = 0 
    opts.datatype {mustBeTextScalar, mustBeMember(opts.datatype, ["single", "double"])} = "single"
    opts.save (1,1) logical = false % deprecated
    opts.tall (1,1) logical = false % read memory friendly genotype calls as tall arrays. Only works in parallel.

    % sample size. If provided, function skips reading header/sample
    % block(s). This option has been implemented for getbulkgeno fnction.
    % Use it with caution, since bgenheader ensures the integrity of bgen
    % file. 
    opts.N (1,1) double

    % indices of samples to be read (a logical vector). This option has
    % been implemented for getbulkgeno fnction. 'N' must be provided as
    % well in this case. Use of logical vector protects the function
    % against unwarranted inputs (i.e. user must have samples in hand in
    % advance and after subsetting, feed the function with a logical
    % vector, instead of samples themselves). Both 'N' and 'samples' can be
    % read using bgenheader function or directly from sample files (as
    % implemented in getbulkgeno).
    opts.samples logical %  {mustBeVector} 
end

fclose('all');

% Check .bgen extension
if ~endsWith(bgen_file,'.bgen')
    bgen_file = bgen_file + ".bgen";
end
% Check bgi file presence
bgi_file = bgen_file + ".bgi";
if ~isfile(bgi_file)
    error('Missing bgi file!')
end

if ~isfield(opts, 'N')
    opts.N = nan;
end

if opts.numWorkers > 0
    if isempty(gcp('nocreate')) && ~isinf(opts.numWorkers)
        distcomp.feature('LocalUseMpiexec', false);
        parpool('local', opts.numWorkers); % Use parallel computations
    end
    parFlag = true;
else
    parFlag = false;
end

if isempty(gcp('nocreate')), opts.tall = false; end

% ----------------------------- Header block ------------------------------
 [N, ~, opts.bversion] = bgenheader(bgen_file, 'sample', false, 'verbose', opts.verbose);

if ~isnan(opts.N)
    N = opts.N;
end

% samples to be used
if ~isfield(opts, 'samples')
    opts.samples = true(N, 1); % fetch all samples
% else % check if sex information are present in 2nd column
%     if size(opts.sample, 2) > 1
%         opts.sexinfo = opts.samples(:, 2);
%         opts.samples = opts.samples(:, 1);
%     end
end

if sum(opts.samples(:, 1)) > N
    error('the length of sample indices cannot be larger than sample size!')
end

if ~all(opts.samples(:, 1)) % subsetting for genotype calls
    N_bed = sum(opts.samples(:, 1));
else
    N_bed = N;
end

% ------------------------ Fetch data from bgi file -----------------------
conn = sqlite(bgi_file);
if isempty(opts.varInfo) % get all variants
    if opts.verbose
        warning(['Reading all varaints may not be applicable',...
            ' due to large memory needed!'])
        fprintf('Fetching bgi file. It may take few minutes...\n')
    end
    bgi_data = fetch(conn, 'SELECT * FROM Variant');
    close(conn);
else
    if opts.verbose
        fprintf('Fetching input variants from bgi file...\n')
    end
    if ~isstring(opts.varInfo)
        opts.varInfo = string(opts.varInfo);
    end
    % find position based queries: they must be in format of pos1-pos2
    checkPos = arrayfun(@(x)split(x,'-'), opts.varInfo, 'uni', false);
    checkPosNumel = cellfun(@numel, checkPos); 
    checkPosIdx1 = checkPosNumel == 2; % condition 1: numel must be 2
    checkPosDouble = cellfun(@double, checkPos, 'uni', false);
    checkPosIdx2 = cell2mat(cellfun(@(x)all(~isnan(x)), checkPosDouble, 'uni', false)); % condition 2: must be double
    checkPosIdx = checkPosIdx1 & checkPosIdx2;

    opts.varInfo_rs = opts.varInfo(~checkPosIdx);
    opts.varInfo_pos = checkPos(checkPosIdx);
    START_str = "SELECT * FROM Variant WHERE ";
    if ~isempty(opts.varInfo_rs)
        opts.varInfo_rs = "('" + join(opts.varInfo_rs,"','") + "')";
        sqlquery_rs = "rsid IN "+opts.varInfo_rs;
    else
        sqlquery_rs = "";
    end
    
    if ~isempty(opts.varInfo_pos)
        for k = 1:numel(opts.varInfo_pos)
            if k > 1
                sqlquery_pos = sqlquery_pos + "(position BETWEEN " + opts.varInfo_pos{k}(1) + " AND " + opts.varInfo_pos{k}(end) + ")";
            else
                sqlquery_pos = "(position BETWEEN " + opts.varInfo_pos{k}(1) + " AND " + opts.varInfo_pos{k}(end) + ")";
            end
            if k < numel(opts.varInfo_pos)
                sqlquery_pos = sqlquery_pos + " OR ";
            end
        end
    else
        sqlquery_pos = "";
    end
    if isempty(opts.varInfo_rs) || isempty(opts.varInfo_pos)
        OR_str = "";
    else
        OR_str = " OR ";
    end
    bgi_data = fetch(conn, START_str+sqlquery_rs+OR_str+sqlquery_pos);
    close(conn)
end

if isempty(bgi_data)
    error('Input variants cannot be found in %s!', bgi_file)
end

% % calculate bgen_var_len matrix: [chr, rsid, a1, a2, pos]
% chr_len_no_zero = cellfun(@(x)numel(x), regexprep(bgi_data(:, 1), '^0', ''));
% bgen_var_len = [cellfun(@(x)numel(x),bgi_data(:,[1,3,5,6])),...
%     arrayfun(@(x)ceil(log10(max(1,double(x)+1))),[bgi_data{:, 2}]')];

% Convert to struct
bgi_hdrs = {'chr', 'pos', 'snp', 'num_allele', 'a1', 'a2',...
    'file_start_position', 'size_in_bytes'};
try % < R2022a
    bgi_data = cell2struct(bgi_data', bgi_hdrs);
catch % >= R2022a
    bgi_data.Properties.VariableNames = bgi_hdrs;
    bgi_data = table2struct(bgi_data, "ToScalar", false);
end

% ----------------------------- Get variant block -------------------------
callData = struct;

if opts.datatype == "single"
    callData.bed = zeros(N_bed, 1, 'single');
else
    callData.bed = zeros(N_bed, 1);
end
if opts.tall, callData.bed = tall(callData.bed); end
callData.I_A = (0); I_A = (0); SWAP_FLAG = (0);

if opts.verbose; tic; end
if parFlag
    bed = ({});
    gpType = opts.gpType;
    datatype = opts.datatype;
    verbose = opts.verbose;
    samples = opts.samples;
    istall = opts.tall;
    bversion = opts.bversion;

    parfor ii = 1:numel(bgi_data)
        fid = fopen(bgen_file, 'r');
        fseek(fid, bgi_data(ii).file_start_position, 'bof');
        g = fread(fid, bgi_data(ii).size_in_bytes, '*uint8', 'l');
        fclose(fid);
        [bed{ii}, I_A(ii, 1), SWAP_FLAG(ii, 1)] = ...
            decompressCalls(g, N, 2, gpType, datatype, samples, ...
            verbose, istall, bversion);
    end

    callData.bed = horzcat(bed{:});
    callData.I_A = I_A;
    clear bed I_A
else
    for ii = 1:numel(bgi_data)
        fid = fopen(bgen_file, 'r');
        fseek(fid, bgi_data(ii).file_start_position, 'bof');
        g = fread(fid, bgi_data(ii).size_in_bytes, '*uint8', 'l');
        fclose(fid);
        [checkbed, callData.I_A(ii, 1), SWAP_FLAG(ii, 1)] = ...
            decompressCalls(g, N, 2, opts.gpType, opts.datatype, ...
            opts.samples, opts.verbose, opts.tall, opts.bversion);
        try
            if strcmp(opts.gpType, 'GP') || strcmp(opts.gpType, 'GT')
                callData.bed(:, 2*ii-1:2*ii) = checkbed;
            else
                callData.bed(:, ii) = checkbed;
            end
        catch % haplotype? diploid and biallelic only
            callData.bed(:, 4*ii-3:4*ii) = checkbed;
        end
        clear checkbed
    end
end

if opts.verbose
    gettime = toc;
    fprintf('Elapsed time is %.5g seconds.\n', gettime)
end
% Read sample file --------------------------------------------------------
% Check sample file
fam = strrep(bgen_file, '.bgen', '.sample');
if ~isstring(fam)
    fam = string(fam);
end
if ~exist(fam, 'file')
    fam = dir('*.sample');
    fam = {fam.name}';   
end
    
if isempty(fam)
    error('could not find any sample file in the current directory!\n')
end
if numel(fam) > 1
    fprintf('Select the sample file:\n')
    continueLoop = true;
    while continueLoop
        try
            for jj = 1:numel(fam)
                fprintf('%d-%s\n', jj, fam{jj})
            end
            selSampleFile = input('? ');
            fam = fam(selSampleFile);
            continueLoop = false;
        catch
            fprintf('Wrong answer! Try again!\n')
        end
    end
end
fid = fopen(fam{1}, 'r');
callData.fam = textscan(fid, '%s %*[^\n]');
fclose(fid);
callData.fam = callData.fam{1};
if N > numel(callData.fam)
    error('%s does not correspond to sample numbers in the bgen file!',...
        callData.fam)
end
diff_N = numel(callData.fam) - N;
eids = double(string(callData.fam(diff_N+1:numel(callData.fam))));
if all(isnan(eids))
    eids = string(callData.fam(diff_N+1:numel(callData.fam)));
else
    eids = int32(eids);
end
callData.fam = eids;
if N_bed < N % keep only samples used in genotype calls
    callData.fam = callData.fam(opts.samples(:, 1));
end
% if any(isnan(callData.fam))
%     error('NaN was found in sample! Check you sample file!')
% end
callData.bim.chr = string({bgi_data.chr})';
callData.bim.pos = double([bgi_data.pos]');
callData.bim.snp = string({bgi_data.snp})';
SWAP_FLAG = logical(SWAP_FLAG);
if all(~SWAP_FLAG) % No need to swap alt/ref alleles
    callData.bim.a1 = categorical(string({bgi_data.a1}))';
    callData.bim.a2 = categorical(string({bgi_data.a2}))';
else
    DEF_REF_AL = categorical(string({bgi_data.a1}))';
    DEF_ALT_AL = categorical(string({bgi_data.a2}))';
    SWAP_REF = DEF_REF_AL(SWAP_FLAG);
    DEF_REF_AL(SWAP_FLAG) = DEF_ALT_AL(SWAP_FLAG);
    DEF_ALT_AL(SWAP_FLAG) = SWAP_REF;
    callData.bim.a1 = DEF_REF_AL; 
    callData.bim.a2 = DEF_ALT_AL;
end
callData.SWAP_FLAG = SWAP_FLAG;
%--------------------------------------------------------------------------
% Export data [DEPRECATED]
if opts.save
    if strcmp(opts.gpType, 'dosage')
        if size(callData.bed, 2) == 1 % Single variant
            callData.dosageFlag = true;
            save(callData.bim(1).snp+'.dos.mat', '-struct', 'callData')
        end
    else % hardcall
        if size(callData.bed, 2) == 2 % Single variant
            save(callData.bim(1).snp, '-struct', 'callData')
        end
    end
end

end % END


%% Subfunctions ===========================================================
function [getG, I_A, SWAP_FLAG] = decompressCalls(genData, N, alleleN, ...
    gpType, datatype, samples, verbose, istall, bversion)

genData = innerzlibdecomp(decodeBgenBlock(genData), bversion);

% note that decompressCallsM_mex only supports 'single' data type and 8
% bits probability (see variable B)
% [getG, I_A, SWAP_FLAG] = decompressCallsM_mex(genData, N, alleleN, gpType, verbose);
% return

ii = 1;
N_v = typecast(genData(ii:ii + 3), 'uint32'); % Number of individuals 
ii = ii + 4;
alleleN_v = typecast(genData(ii:ii + 1), 'uint16');% number of alleles
% ii = ii + 2;
if N_v ~= N || alleleN_v ~= alleleN
    error('Number of individuals or alleles discrepancy in decompressed calls!')
end
genData(1:6) = [];
% ii = 1;
% P_min = genData(ii); % The minimum ploidy Pmin of samples in the row
% P_max = genData(ii + 1); % The maximum ploidy Pmax of samples in the row
genData(1:2) = [];

% sample subsetting
if ~all(samples(:, 1))
    N = sum(samples(:, 1));
    samplesIdx = [samples(:, 1); true; true; repelem(samples(:, 1), alleleN)]; % ploidy + 2(trailing block) + trailing
    if numel(samplesIdx) ~= numel(genData)
        error('something went wrong with sample indexing!')
    end
    genData = genData(samplesIdx); % only keep these samples
end

genData_trailing = genData(N + 1:end);

% Ploidy and missingness should be checked --------------------------------
% genData = [];

genData = genData(1:N)'; % ploidy and missingness
% Note:  Ploidy      = least significant 6 bits (right-msb)
%        Missingness = most significant bit (left-msb), if 1 => missing
% genDataPloidy = de2bi(genData, 6, 'right-msb');
% genDataMissing = de2bi(genData, 8, 'left-msb'); genDataMissing = logical(genDataMissing(:,end));
% genDataMissing = dec2bin(genData, 8);
genDataMissing = int2bit(genData, 8);
genDataMissing = genDataMissing(8:8:end) > 0;
%--------------------------------------------------------------------------

phaseFlag = genData_trailing(1); 
% If Phased=1 the row stores one probability per allele
% (other than the last allele) per haplotype (e.g. to represent phased data).
% If Phased=0 the row stores one probability per possible genotype 
% (other than the 'last' genotype where all alleles are the last allele), 
% to represent unphased data.

if ~phaseFlag
    B = double(genData_trailing(2));  % number of bits used to store each probability in this row
    genData_trailing(1:2) = [];
    genData_trailing = typecast(genData_trailing, "uint"+B);
    numG = numel(genData_trailing) / alleleN;
    if datatype == "single"
        getG = zeros(numG, alleleN + 1, 'single');
    else
        getG = zeros(numG, alleleN + 1);
    end

    for ii = 1:alleleN
        if datatype == "single"
            getG(1:numG, ii) = single(genData_trailing(ii:alleleN:end)) ./ (2^B - 1);
        else
            getG(1:numG, ii) = double(genData_trailing(ii:alleleN:end)) ./ (2^B - 1);
        end
    end
    % Colex order for 2 alleles: 11,12,22
    getG(1:numG, alleleN + 1) = 1 - sum(getG,2);
elseif phaseFlag == 1
    if verbose, fprintf('phased data.\n'); end
    B = double(genData_trailing(2));  % number of bits used to store each probability in this row
    genData_trailing(1:2) = [];
    genData_trailing = typecast(genData_trailing, "uint"+B);
    hapl = alleleN;
    if hapl ~= 2
        error('bgenreader cannot handle ploidy > 2 for phased data!')
    end
    numG = numel(genData_trailing) / hapl; % only biallelic tested
    haplReal = 2*hapl;   
    if datatype == "single"
        getG = zeros(numG, haplReal, 'single'); 
    else
        getG = zeros(numG, haplReal); 
    end
    haplFlag = true; ii = 1; jj = 1;
    while haplFlag % for haplotypes: P11 P12(absent) P21 P22(absent)
        if mod(jj, 2)
            getG(1:numG, jj) = double(genData_trailing(ii:hapl:end)) ./ (2^B - 1);
            ii = ii + 1;
        else
            getG(1:numG, jj) = 1 - getG(1:numG, jj-1);
        end
        jj = jj + 1;
        if jj > haplReal
            haplFlag = false;
        end
    end
else
    error('Wrong phaseFlag value!')
end

if phaseFlag
    if datatype == "single"
        I_A = single(NaN);
    else
        I_A = NaN;
    end
    SWAP_FLAG = false;
    if istall, getG = tall(getG); end
    return
end


E = getG(:,2) + 2.*getG(:,3); %expected genotype
F = getG(:,2) + 4.*getG(:,3);
coeff = 1;

diff_EF = F - (E.^2);
Teta = sum(E(:))/(2*N); %MLE for the population minor allele frequency
if Teta == 0 || Teta == 1
    if datatype == "single"
        I_A = single(1);
    else
        I_A = 1;
    end
else
    % I_A = (sum(E.^2)/numel(E) - mean(E)^2)/(mean(E).*(1-mean(E)/2)); % MaCH
    I_A = coeff*(1 - sum(diff_EF(:))./(2*N*Teta*(1-Teta))); % IMPUTE
end

% Note: only applicable when alleleN = 2
% Alternale allele dosage (based on PLINK2)
getG_CHECK = getG(:,2) + 2.*getG(:,3);
getG_SUM = sum(getG, 2);
missing_value = getG_SUM < 1;
if strcmp(gpType, 'dosage')
    if (sum(getG_CHECK(:))/(2*numel(getG_CHECK))) > 0.5
        getG = getG(:,2) + 2.*getG(:,1);
        SWAP_FLAG = true;
    else
        getG = getG(:,2) + 2.*getG(:,3);
        SWAP_FLAG = false;
    end
    getG(missing_value) = -2;
    getG(genDataMissing) = -2;

    if istall, getG = tall(getG); end
    return
elseif strcmp(gpType, 'GP') % Genotype probabilities
    % This can be used to run different regression models as described in
    % ProbABEL paper https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-134
    if (sum(getG_CHECK(:))/(2*numel(getG_CHECK))) > 0.5
        getG = [getG(:, 2),getG(:, 1)];
        SWAP_FLAG = true;
    else
        getG = getG(:, 2:3);
        SWAP_FLAG = false;
    end
    getG(missing_value) = -2;
    getG(genDataMissing) = -2;

    if istall, getG = tall(getG); end
    return
end

if (sum(getG_CHECK(:))/(2*numel(getG_CHECK))) > 0.5
    genotype_dosage = getG(:,2) + 2.*getG(:,1);
    SWAP_FLAG = true;
else
    genotype_dosage = getG(:,2) + 2.*getG(:,3);
    SWAP_FLAG = false;
end
d_11 = (genotype_dosage <= 0.1); % [0,0.1]
d_12 = ((genotype_dosage >= 0.9) & (genotype_dosage <= 1.1)); % [0.9,1.1]
d_13 = (genotype_dosage >= 1.9); % [1.9,2.0]
d_miss = ~logical(d_11 + d_12 + d_13);

% Note: from GTOOL: In the GEN format each SNP is represented as a set of 
% three probabilities which correspond to the allele pairs AA,AB,BB. If the
% largest of the probabilities is over the threshold specified by
% --threshold, then the genotype in the PED file is expressed as the
% corresponding allele pair. The genotypes are expressed as pairs of 1,2,0
% where 1 corresponds to alleleA from the GEN file and 2 corresponds to 
% alleleB. If none of the probabilities are over the calling threshold then 
% the pair is unknown, 0 0.

% Note: bcftools get max prob as call: https://github.com/samtools/bcftools/blob/b935f56a842383479e72ca9eed845411f90f116b/plugins/tag2tag.c

% call threshol
call_threshold = 0.9;
[nr, nc] = find(getG >= call_threshold);
f_11 = nr(nc == 1);
f_12 = nr(nc == 2);
f_22 = nr(nc == 3);

% Option 1: PLINK2 + --import-dosage-certainty 0.9
if datatype == "single"
    getG = -2*ones(size(getG, 1), 2, 'single');
else
    getG = -2*ones(size(getG, 1), 2);
end

if SWAP_FLAG
    getG(f_11, 1) = 1; getG(f_11, 2) = 1;
    getG(f_12, 1) = 0; getG(f_12, 2) = 1;
    getG(f_22, 1) = 0; getG(f_22, 2) = 0;
else
    getG(f_11, 1) = 0; getG(f_11, 2) = 0;
    getG(f_12, 1) = 0; getG(f_12, 2) = 1;
    getG(f_22, 1) = 1; getG(f_22, 2) = 1;
end
getG(d_miss, 1) = -2; getG(d_miss, 2) = -2;

if istall, getG = tall(getG); end

end

%% ------------------------------------------------------------------------
function [data, id_info] = decodeBgenBlock(data)
% implemented for BGEN layout 2 v1.2/v1.3
offset = 1;

for i = 1:3
    sz = typecast(data(offset:offset + 1), 'uint16');
    offset = offset + 2;
    id_info = cast(data(offset:offset + sz - 1), 'char');
    offset = offset + numel(id_info);
end

offset = offset + 4;

num_allele = typecast(data(offset:offset+1), 'uint16');
offset = offset + 2;

for i = 1:num_allele
    a_sz = typecast(data(offset:offset + 3), 'uint32');
    offset = offset + 4;
    alleles = cast(data(offset:offset + a_sz - 1), 'char');
    offset = offset + numel(alleles);
end

data(1: offset + 7) = [];

end

%% ------------------------------------------------------------------------
function genData =  innerzlibdecomp(genData, bversion)
%@09AUG2024: added zstd decompression support
if bversion == 1
    buffer = java.io.ByteArrayOutputStream();
    zlib = java.util.zip.InflaterOutputStream(buffer);
    zlib.write(genData, 0, numel(genData));
    zlib.close();
    genData = typecast(buffer.toByteArray(), 'uint8');
elseif bversion == 2
    genData = zipmat(genData, 0, 'zstd', 1, 0, 1);
end

end

