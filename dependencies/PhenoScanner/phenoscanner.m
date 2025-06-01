function [S, snp] = phenoscanner(key, varargin)

% fetches summary stat from PhenoScanner database:
% http://www.phenoscanner.medschl.cam.ac.uk/
% 
% Note: this function is still under development.
% 
% INPUTS
% key (required):        A scalar string correponding to method.
%                        entryPoint.
% method (optional): Scalar string or char array for entry points to be
%                        used: "gene", "region", "snp" (default) or
%                        "trait".

% input parser ------------------------------------------------------------
[S, snp] = deal([]);
p = inputParser;
p.CaseSensitive = false;
p.StructExpand = true;
addRequired(p, 'key', ...
    @(x) validateattributes(x,{'char','string'}, {'nonempty'}));
addOptional(p, 'method', 'snp',...
    @(x) (ischar(x) || isstring(x)) && ...
    ismember(x, ["snp", "region", "gene", "trait"]) ...
    && isscalar(string(x)));
addParameter(p, 'r2', 0.8,...
    @(x) validateattributes(x,{'numeric'},...
    {'nonempty', 'scalar', '>', 0, '<=', 1}));
addParameter(p, 'p', 1e-5,...
    @(x) validateattributes(x,{'numeric'},...
    {'nonempty', 'scalar', '>', 0, '<=', 1}));
addParameter(p, 'build', 37,...
    @(x) isnumeric(x) && (x == 37 || x == 38) && isscalar(x));
addParameter(p, 'catalogue', "GWAS",...
    @(x) (ischar(x) || isstring(x)) && ...
    ismember(x, ["GWAS", "All", "eQTL", "pQTL", "mQTL", "methQTL"]) ...
    && isscalar(string(x)));
addParameter(p, 'proxies', "None",...
    @(x) (ischar(x) || isstring(x)) && ...
    ismember(x, ["None", "AFR", "AMR", "EAS", "EUR", "SAS"]) ...
    && isscalar(string(x)));

p.parse(key, varargin{:});

p = p.Results;

if ~isstring(p.key)
    p.key = string(p.key);
end

% -------------------------------------------------------------------------
fprintf('\n********************* PhenoScanner *************************\n')

switch string(p.method)
    case "trait"
        p.method = "?query=";
    case "snp"
        p.method = "api/?snpquery=";
    case "region"
        p.method = "api/?regionquery=";   
    case "gene"
        p.method = "api/?genequery=";
end

getURL = "http://www.phenoscanner.medschl.cam.ac.uk/";
headers = {'Content-Type' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', 10000);

getURL = getURL + p.method + p.key + "&catalogue=" + p.catalogue + ...
    "&p=" + p.p + "&proxies=" + p.proxies + "&r2=" + p.r2 +...
    "&build=" + p.build;
fprintf('connecting to PhenoScanner...\n')
try
    S = webread(getURL, options);
catch ME
    fprintf('phenoscanner error: %s\n', ME.message)
    S = [];
    return
end
fprintf('parsing fetched data...\n')

if p.method == "?query="
    S = parseTrait(S);
else
    if isfield(S, 'error')
        fprintf('%s :(\n', S.error)
        S = [];
        return
    else
        [S, snp] = parseAPIscanner(S);
    end
end

end % END

%% subfunctions ===========================================================
function S = parseTrait(S)
S = regexp(S, '(?<=<h2>Top 1000 associations</h2>).*', 'match');
hdr = regexp(S{1}, '(?:th\s{1})(.*?)[>](?<hdr>.*?)(</th>)', 'names');
hdr = string({hdr.hdr})';
S = regexp(S{1}, '(<td>)(?<body>.*?)(</td>)', 'names');
S = string({S.body});
S = reshape(S, numel(hdr), numel(S)/numel(hdr))';
chrpos = S(:, 3); % separate chr and pos
chrpos = split(chrpos, ':');
chr = replace(chrpos(:, 1), 'chr', '');
chrpos = chrpos(:, 2);
S(:, 3) = []; hdr(3) = [];
hdr = [hdr(1:2); "chr"; "pos"; hdr(3:end)];
S = [S(:, 1:2), chr, chrpos, S(:, 3:end)];
S = array2table(S, 'VariableNames', hdr);
S.pos = double(S.pos); S.N = double(S.N); S.Beta = double(S.Beta);
S.P = double(S.P);
end % END

%% ------------------------------------------------------------------------
function [S, meta] = parseAPIscanner(S)
fnames = fieldnames(S);
fnames(ismember(fnames, 'results')) = [];
meta = S.(fnames{1}); meta = string(horzcat(meta{:}));
meta = array2table(meta(:, 2), 'RowNames', meta(:, 1));
meta.Properties.VariableNames = {'info'};

S = S.results;
S = string(horzcat(S{:})');
hdr = S(1, :); S(1, :) = [];
hg19_coordinates = hdr == "hg19_coordinates";
hg38_coordinates = hdr == "hg38_coordinates";
f_coor = find(hg19_coordinates | hg38_coordinates);
hdrpos = hdr(f_coor); hdrpos = replace(hdrpos, "_coordinates", "_pos");
hdr(f_coor) = [];

if isempty(S)
    hdr = [hdr(1:min(f_coor)-1), "chr", hdrpos, hdr(min(f_coor):end)];
    S = array2table(repmat("-", 1, numel(hdr)));
    S.Properties.VariableNames = hdr;
    S.direction = [];
    S.hg19_pos = double(S.hg19_pos); S.hg38_pos = double(S.hg38_pos);
    S.beta = double(S.beta); S.se = double(S.se); S.p = double(S.p);
    S.n = double(S.n); S.n_cases = double(S.n_cases);
    S.n_controls = double(S.n_controls); S.n_studies = double(S.n_studies);
    return
end

chrpos = S(:, f_coor); S(:, f_coor) = [];
chrpos = split(chrpos, ':');
chr = replace(chrpos(:, 1, 1), 'chr', '');
chrpos = chrpos(:, :, 2);
S = [S(:, 1:min(f_coor)-1), chr, chrpos, S(:, min(f_coor):end)];
hdr = [hdr(1:min(f_coor)-1), "chr", hdrpos, hdr(min(f_coor):end)];
S = array2table(S, 'VariableNames', hdr);
S.hg19_pos = double(S.hg19_pos); S.hg38_pos = double(S.hg38_pos);
S.beta = double(S.beta); S.se = double(S.se); S.p = double(S.p);
S.n = double(S.n); S.n_cases = double(S.n_cases);
S.n_controls = double(S.n_controls); S.n_studies = double(S.n_studies);
S.direction = [];
end % END