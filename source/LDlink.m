function R = LDlink(snps, varargin)
% Under development.
% Implements methods https://ldlink.nci.nih.gov/?tab=apiaccess

p = inputParser;
p.CaseSensitive = false;
p.StructExpand = true;
addRequired(p, 'snps', @(x)validateattributes(x,{'cellstr','string'},...
    {'nonempty'}));
addParameter(p, 'pop', 'GBR', @(x)validateattributes(x,{'char','string'},...
    {'nonempty', 'scalartext'})); % pop can be a string array (multiple pops); this must be changed in future
addParameter(p, 'r2_d', 'r2', @(x) strcmp(x, {'r2', 'd'}));
addParameter(p, 'APItoken', 'YOURTOKEN', ...
    @(x)validateattributes(x,{'char','string'},...
    {'nonempty', 'scalartext'}));
p.parse(snps, varargin{:});
p = p.Results;

if numel(snps) < 2
    error('rslist must be > 2!')
elseif numel(snps) > 500 % use openGWAS
    fprintf('\nnumel snps > 500, will use OpenGWAS API instead.\n')
    R = call_opengwas(p);
else % use ldlink
    R = call_opengwas(p);
    % R = call_ldlink(p);
end

end %END

%%
function R = call_ldlink(p)
getURL = "https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?&token="+...
    string(p.APItoken);
data = struct('snps', join(p.snps, newline), 'pop', p.pop, 'r2_d',p.r2_d);
headers = {'Content-Type' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', 5000);

r = webwrite(getURL, data, options);
R = struct;
snps = regexp(r, '[\s](\w+)[\s][^\n]*', 'tokens')';
R.snp = string(vertcat(snps{:})); R.snp(1) = [];
r = string(split(r)); r(1) = []; r = reshape(r, numel(snps), numel(snps));
R.ld = double(r(1:end-1, 2:end)); clear r snps
end

%%
function R = call_opengwas(p)
data = struct('rsid',p.snps);
headers = {'Content-Type' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', 5000);
url = "http://gwas-api.mrcieu.ac.uk/ld/matrix";
R = webwrite(url, data, options);
R.matrix = horzcat(R.matrix{:});
R.matrix = double(string(R.matrix)).^2;

underlineIdx = cellfun(@(x)x(end-1), regexp(R.snplist, '_'));
R.snp = cellfun(@(x,y)x(1:y-1), R.snplist, num2cell(underlineIdx), 'UniformOutput', false);
if ~all(ismember(R.snp, p.snps))
    error('call_opengwas:failed to regenerate rsids from API!')
end

[idx1, idx2] = ismember(p.snps, R.snp); idx2(idx2<1) = [];
R.snplist = R.snplist(idx2); R.snp = R.snp(idx2);
R.matrix = R.matrix(:, idx2);
R.matrix = R.matrix(idx2, :);

newmat = spalloc(numel(p.snps), numel(p.snps), 0);
newmat(idx1, idx1) = R.matrix;
R = rmfield(R, {'matrix', 'snplist'});
R.ld = full(newmat); clear newmat
R.snp = string(R.snp);
end
