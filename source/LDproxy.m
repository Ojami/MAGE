function res = LDproxy(snps, opts)
% API access to LDproxy tool: https://ldlink.nci.nih.gov/?tab=apiaccess
% snps can be as rsid or chrN:POS
% Oveis Jamialahmadi, Sahlgrenska Akademy, April 2023.

arguments
    snps {mustBeNonzeroLengthText, mustBeText}
    opts.genome_build {mustBeMember(opts.genome_build, ["grch37", "grch38", "grch38_high_coverage"])} = "grch37"
    opts.pop {mustBeMember(opts.pop, {'ALL';'AFR';'YRI';'LWK';'GWD';'MSL';'ESN';'ASW';'ACB';'AMR';'MXL';'PUR';'CLM';'PEL';'EAS';'CHB';'JPT';'CHS';'CDX';'KHV';'EUR';'CEU';'TSI';'FIN';'GBR';'IBS';'SAS';'GIH';'PJL';'BEB';'STU';'ITU'})} = 'GBR' % can be accessed via list_pop in R
    opts.r2_d {mustBeMember(opts.r2_d, ["r2", "d"])} = "r2"
    opts.window {mustBeInRange(opts.window, 0, 1e6)} = 5e4
    opts.token {mustBeTextScalar}
    
    opts.timeout (1,1) double = 5e3
    opts.verbose (1,1) logical = true
end

headers = {'Content-Type' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', opts.timeout);

fis = setdiff(fieldnames(opts), "verbose");
res = cell(numel(snps), 1);
for k = 1:numel(snps)
    url = "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=" + snps(k);
    for j = 1:numel(fis)
        url = url + "&" + fis(j) + "=" + opts.(fis(j));
    end

    tmp = webread(url, options);
    if isstruct(tmp)
        if opts.verbose, disp(tmp.error); end
        continue
    end
    tmp = string(splitlines(tmp));
    tmp(tmp == "") = [];
    tmp = split(tmp, char(9));
    tmp = array2table(tmp(2:end, :), VariableNames=tmp(1, :));
    tmp.Query(:) = snps(k);
    tmp = convertvars(tmp, ["MAF", "Distance", "Dprime", "R2"], @double);
    res{k} = tmp;
end
res(cellfun(@isempty, res)) = [];
if ~isempty(res), res = vertcat(res{:}); end

end % END
