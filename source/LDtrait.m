function res = LDtrait(snps, opts)
% API access to LDtrait tool: https://ldlink.nci.nih.gov/?tab=apiaccess
% snps can be as rsid or chrN:POS
% Oveis Jamialahmadi, Sahlgrenska Akademy, April 2023.


arguments
    snps {mustBeNonzeroLengthText, mustBeText}
    opts.genome_build {mustBeMember(opts.genome_build, ["grch37", "grch38", "grch38_high_coverage"])} = "grch37"
    opts.pop {mustBeMember(opts.pop, {'ALL';'AFR';'YRI';'LWK';'GWD';'MSL';'ESN';'ASW';'ACB';'AMR';'MXL';'PUR';'CLM';'PEL';'EAS';'CHB';'JPT';'CHS';'CDX';'KHV';'EUR';'CEU';'TSI';'FIN';'GBR';'IBS';'SAS';'GIH';'PJL';'BEB';'STU';'ITU'})} = 'GBR' % can be accessed via list_pop in R
    opts.r2_d {mustBeMember(opts.r2_d, ["r2", "d"])} = "r2"
    opts.r2_d_threshold {mustBeInRange(opts.r2_d_threshold, 0, 1)} = 0.1; 
    opts.window {mustBeInRange(opts.window, 0, 1e6)} = 5e4
    opts.token {mustBeTextScalar}
    
    opts.timeout (1,1) double = 5e4
    opts.verbose (1,1) logical = true
end

getURL = "https://ldlink.nci.nih.gov/LDlinkRest/ldtrait?token="+opts.token;
data = struct('snps', join(snps, newline), 'pop', join(opts.pop, "+"), ...
    'r2_d', string(opts.r2_d), 'r2_d_threshold', string(opts.r2_d_threshold), ...
    'window', string(opts.window));
headers = {'Content-Type' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', opts.timeout);

t1 = tic;
res = webwrite(getURL, data, options);
res = string(splitlines(res));
res(res == "") = [];
res = split(res, char(9));
hdr = res(1, :); 
res(1, :) = [];
res = array2table(res, VariableNames=hdr);
res = convertvars(res, ["D'", "R2", "P-value"], @double);
timeSpent = toc(t1);
if opts.verbose; fprintf('Elapsed time: %.2f seconds.\n', timeSpent); end

end % END
