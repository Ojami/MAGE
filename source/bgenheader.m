function [N, samples, CompressedSNPBlocks] = bgenheader(bgen_file, opts)
% This is a subfunction of bgenreader, and reads header block of bgen file.
% Can read sample block (if present) as well.
% 
% Oveis Jamialahmadi, Sahlgrenska academy, Feb 2022.
% @09AUG2024: now supports also BGEN v1.3 with zstd compression.

arguments
    bgen_file {mustBeTextScalar, mustBeFile}
    opts.sample (1,1) logical = false % returns sample file as well
    opts.verbose (1,1) logical = false
end

if ~isfield(opts, 'sample')
    opts.sample = false;
end

if ~opts.sample; samples = string; end

fid = fopen(bgen_file, 'r');
hdrData = fread(fid, 4, 'uint32', 'l');
% varoffset = hdrData(1);
% L_H =  hdrData(2); % Length of header block
% M = hdrData(3); % Number of variants
N = hdrData(4); % Number of samples
dec_magic = fread(fid, 4, 'uint8', 'l');
magic_number = char(dec_magic'); % BGEN magic number
if ~strcmp(magic_number, 'bgen')
    if sum(dec_magic)
        fclose(fid);
        error('ERROR: magic number cannot be verified!')
    end
end
% if L_H ~= 20
%     free_data_area =  fread(fid,L_H - 20, '*char', 'l')';
% else
%     free_data_area = [];
% end
header_flags = fread(fid, 1, 'uint32', 'l');
header_flags = int2bit(header_flags, 32, false); % de2bi(header_flags, 32);
CompressedSNPBlocks = bit2int(header_flags(1:2), 2 ,0); % double(bi2de(header_flags(1:2)));
if opts.verbose
    fprintf("SNP blocks layout is %d\n", bit2int(header_flags(3:6), 4 ,0))
end

if CompressedSNPBlocks == 0
    if opts.verbose
        fprintf('SNP block probability data is not compressed\n')
    end
elseif CompressedSNPBlocks == 1 
    if opts.verbose
        fprintf('SNP block probability data is compressed using zlib\n')
    end
elseif CompressedSNPBlocks == 2
    if opts.verbose
        fprintf('SNP block probability data is compressed using zstd\n')
    end
else
    fclose(fid);
    error('Unexpected value for CompressedSNPBlocks')
end
Layout = bit2int(header_flags(3:6), 4 ,0); % double(bi2de(header_flags(3:6)));
if Layout ~= 2
    fclose(fid);
    error('bgenreader can only read layout 2!\n')
end

SampleIdentifiers = bit2int(header_flags(32), 1 ,0); % double(bi2de(header_flags(32))); 
if SampleIdentifiers == 0
    if opts.verbose
        fprintf('Sample identifiers are not stored in this file\n')
    end
    samples = missing;
elseif SampleIdentifiers == 1
    if opts.verbose
        fprintf('Sample identifier block follows the header\n')
    end
    if opts.sample % read sample block
        sblock = fread(fid, fread(fid, 1, 'uint32', 'l'), '*uint8', 'l');
        samples = parseSblock(sblock, N, fid);
    end
else
    fclose(fid);
    error('Unexpected value for SampleIdentifiers')
end
fclose(fid);
end % END

%% subfunctions ===========================================================
function samples = parseSblock(sblock, N, fid)
idx = 1;
Ns = typecast(sblock(idx:idx + 3), 'uint32');
if Ns ~= N
    fclose(fid);
    error('Mismatched sample size between bgen header and sample blocks!')
end
idx = idx + 4;

samples = strings(N, 1);
% uint16 max is 65535, so unless there is a better solution, I convert it
% to double to avoid errors.
for i = 1:N
    Lsi = double(typecast(sblock(idx:idx + 1), 'uint16'));
    idx = idx + 2;
    samples(i) = char(sblock(idx:idx + Lsi - 1))';
    idx = idx + Lsi;
end

end