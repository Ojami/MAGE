function res = LDmatrix(snps, opts)
% API access to LDmatrix tool: https://ldlink.nci.nih.gov/?tab=apiaccess
% snps can be as rsid or chrN:POS
% Oveis Jamialahmadi, Sahlgrenska Akademy, April 2023.

arguments
    snps {mustBeNonzeroLengthText, mustBeText}
    opts.genome_build {mustBeMember(opts.genome_build, ["grch37", "grch38", "grch38_high_coverage"])} = "grch37"
    opts.pop {mustBeMember(opts.pop, {'ALL';'AFR';'YRI';'LWK';'GWD';'MSL';'ESN';'ASW';'ACB';'AMR';'MXL';'PUR';'CLM';'PEL';'EAS';'CHB';'JPT';'CHS';'CDX';'KHV';'EUR';'CEU';'TSI';'FIN';'GBR';'IBS';'SAS';'GIH';'PJL';'BEB';'STU';'ITU'})} = 'GBR' % can be accessed via list_pop in R
    opts.r2_d {mustBeMember(opts.r2_d, ["r2", "d"])} = "r2"
    opts.token {mustBeTextScalar}
    
    opts.timeout (1,1) double = 5e3
end

getURL = "https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?token="+opts.token;
data = struct('snps', join(snps, newline), 'pop', join(opts.pop, "+"), ...
    'r2_d', string(opts.r2_d));
headers = {'Content-Type' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', opts.timeout);
try
    tmp = webwrite(getURL, data, options);
    tmp = string(splitlines(tmp));
catch
    disp(tmp)
    res.snps = [];
    res.ld = [];
    return
end


tmp(tmp == "") = [];
tmp = split(tmp, char(9));
res.snps = tmp(2:end, 1);
res.ld = double(tmp(2:end, 2:end));
end % END
