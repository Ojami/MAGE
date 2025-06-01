function res = LDexpress(snps, opts)
% API access to ldepxress tool using LDlinkR package
% Search if a list of variants (or variants in LD with those variants) is
% associated with gene expression in multiple tissue types
% documentation: https://cran.r-project.org/web/packages/LDlinkR/index.html
% 
% Oveis Jamialahmadi, Sahlgrenska Akademy, June 2021.
% 
% @23/12/2021: 'verbose' flag was added.

arguments
    snps {mustBeNonzeroLengthText, mustBeText} % donno exactly what's the limit (300 maybe)
    opts.pop {mustBeMember(opts.pop, {'ALL';'AFR';'YRI';'LWK';'GWD';'MSL';'ESN';'ASW';'ACB';'AMR';'MXL';'PUR';'CLM';'PEL';'EAS';'CHB';'JPT';'CHS';'CDX';'KHV';'EUR';'CEU';'TSI';'FIN';'GBR';'IBS';'SAS';'GIH';'PJL';'BEB';'STU';'ITU'})} = 'GBR' % can be accessed via list_pop in R
    opts.tissue {mustBeTextScalar} = "ALL" % can be accessed via list_getex_tissues() for full list
    opts.r2d {mustBeMember(opts.r2d, ["r2", "d"])} = "r2"
    opts.r2d_threshold {mustBeInRange(opts.r2d_threshold, 0, 1)} = 0.3; % default in LDlinkR is 0.2, this was changed to be more precise and fast
    opts.p_threshold {mustBeInRange(opts.p_threshold, 0, 1)} = 0.1;
    opts.win_size {mustBeInRange(opts.win_size, 0, 1e6)} = 5e5
    opts.token {mustBeTextScalar} 
    opts.file {mustBeTextScalar} = "LDexpress.output"
    opts.verbose (1,1) logical = true
end

% generate unique random string for parallel calls
charSet = ['a':'z',upper('a':'z'),'0':'9'];
randStr = string(charSet(randi(numel(charSet),1, 10)));
maxloop = 5;
while exist("LDexpress."+randStr+".r", 'file')
    randStr = string(charSet(randi(numel(charSet),1, 10)));
    maxloop = maxloop + 1;
    if maxloop > 5
        error('cannot generate unique random string for LDexpress :(')
    end
end

currDir = strrep(pwd, '\', '/');
Rcode = string;
Rcode(1, 1) = "setwd('" + currDir + "')";
Rcode(2, 1) = "library(LDlinkR)";
Rcode(4, 1) = "out <- LDexpress(snps = " + getvector(snps) + ...
    "pop = " + getvector(opts.pop) + ...
    "tissue = " + getvector(opts.tissue) + ...
    "r2d = " + getrvar(opts.r2d) + ", r2d_threshold = " + opts.r2d_threshold + ...
    ", p_threshold = " + opts.p_threshold + ...
    ", win_size = " + opts.win_size + ...
    ", token = " + getrvar(opts.token) + ...
    ", file = " + getrvar(opts.file) + ")";
writematrix(Rcode, "LDexpress"+randStr+".r", 'QuoteStrings', false, 'FileType','text')
[timeM2R, failedM2R] = MATLAB2Rconnector( "LDexpress"+randStr+".r",...
    'delr', true, 'printStdout', opts.verbose);

if failedM2R
    res = [];
    return
end

res = readtable(opts.file, 'FileType', 'text',...
    'VariableNamingRule', 'preserve', 'NumHeaderLines', 0, ...
    'TextType', 'string');
res = sortrows(res, 'P_value', 'ascend');
delete(opts.file)
if opts.verbose; fprintf('Elapsed time: %.2f seconds.\n', timeM2R); end

end % END

%% subfunctions
function r = getvector(r, endFlag)
% creates R vector format
if nargin < 2
    endFlag = false;
end

r = 'c("' + join(string(r), '","') + '")';
if ~endFlag
    r = r + ",";
end
end

function r = getrvar(r)
% adds quotes for strings
if ~isnumeric(r)
    r = '"' + string(r) + '"';
end
end