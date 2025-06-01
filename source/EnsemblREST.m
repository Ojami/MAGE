function outdata = EnsemblREST(indata, rest_method, opts) 
% EnsemblREST will cover useful functions from Ensembl REST API Endpoints:
% https://rest.ensembl.org/
% Note: for POST method, maximum size is 1000 queries
%
% INPUTS:
% indata: a string array of input data, depending on rest_method
% rest_method: a character array based on Ensembl Endpoints.

% [EDITED]: biotools.fr API was added for ID conversions. 6 August 2021.
% [EDITED]: res_vep was modified with single region query. 23 December
%           2021.
% [EDITED]: 'map' endpoint was implemented to convert the co-ordinates of
%           one assembly to another. 'refGenome' is treated as the source
%           assembly. Input format can either be chr:start-end or chr:start
%           24 December 2021.
% [EDITED]: Ensemble has removed 'eqtl' endpoint and replaced by EBI API
%           (not official). Input (only single query) should be in
%           chr$pos$ref$alt format, $ can be :|-|_. 24 December 2021.
% [EDITED]: a bug was fixed. 07 January 2022.
% [EDITED]: a bug was fixed with 'lookup' method. 31 January 2022.
% [EDITED]: 'nearest' method was added to find the nearest gene to a
%           variant. Query variant should be in format of chr$pos, where $
%           can be :|-|_. 12 June 2022.
% [EDITED]: 'xref' option was implemented. 20 June 2023.
% @31JAN2025: 'convert' was added, which is equivalent to xref_name
% 
% TODO: 
% https://rest.ensembl.org/overlap/translation/ENSP00000400821?content-type=application/xml
% https://rest.ensembl.org/map/translation/ENSP00000400821/1..2000?content-type=application/xml


arguments
    indata {mustBeNonzeroLengthText, mustBeText}
    rest_method {mustBeMember(rest_method, {'vep', 'lookup', 'eqtl', 'biotoolsfr', 'map', 'nearest', 'xref', 'convert'})}
    opts.refGenome {mustBeMember(opts.refGenome, {'37', '38'})} = '38';
    opts.forcePost (1,1) logical = false
    opts.getsymbol (1,1) logical = true % get gene symbol for eqtl gene IDs
    opts.geneSymbol (1,1) logical = false; % for 'lookup' method: if true, queries are gene symbols instead of Ensembl IDs
    opts.parseFlag (1,1) logical = true; % extract useful fields from raw structures (only works with 'lookup')
    opts.verbose (1,1) logical = false
end

if numel(indata) >= 1000 && rest_method ~= "convert"
    error("EnsemblREST:indata size must be < 1000!")
end

% Check indata size:
if numel(indata) > 1
    postFlag = true;
else
    postFlag = false;
end

if opts.forcePost
    postFlag = true;
end

if iscell(indata)
    indata = string(indata);
end

switch rest_method
    case 'lookup'
        outdata = rest_lookup(indata, postFlag, opts.refGenome, opts.geneSymbol, opts.parseFlag);
    case 'vep'
        outdata = rest_vep(indata, postFlag, opts.refGenome);
    case 'eqtl'
        outdata = rest_eqtl(indata, opts);
    case 'biotoolsfr'
        outdata = rest_biotoolsfr(indata);
    case 'map'
        outdata = rest_map(indata, opts);
    case 'nearest'
        outdata = rest_nearest(indata, opts);
    case 'xref'
        outdata = rest_xref(indata, opts);
    case 'convert'
        outdata = rest_convert(indata, opts);
    otherwise
        error('Method not supported!')
end

end % END

%% lookup using external reference ========================================
function outdata = rest_convert(indata, opts)
% should be gene symbols

if strcmp(opts.refGenome, '37')
    getURL = "https://grch37.rest.ensembl.org/xrefs/name/human/";
else
    getURL = "https://rest.ensembl.org/xrefs/name/human/";
end
headers = {'Content-Type' 'application/xml'};
options = weboptions('HeaderFields', headers, 'Timeout', 5000);

indata = string(indata);
outdata = cell(numel(indata), 1);
for k = 1:numel(indata)
    tmp = webread(getURL + indata(k), options);
    if isempty(tmp), continue; end
    
    if iscell(tmp)
        tmp = vertcat(tmp{:});
    end

    try
        tmp = struct2table(tmp);
    catch
        tmp = struct2table(tmp, AsArray=true);
    end

    outdata{k} = aggregateCells(tmp);
    clear tmp
end

outdata(cellfun(@isempty, outdata)) = [];
outdata = vertcat(outdata{:});

end

%% cross ref using external symbol ========================================
function outdata = rest_xref(indata, opts)
% useful when a gene name and id are deprecated in one reference genome,
% and gene name can be used to find it in the other reference genome.
% example: HIST1H4B (ENSG00000124529) is deprecated in GRCh38. To find the
% ID associated in GRCh38, one can use this option.
if numel(indata) > 1
    error('EnsemblREST::xref can only work with single query!')
end

if strcmp(opts.refGenome, '37')
    getURL = "https://grch37.rest.ensembl.org/xrefs/symbol/homo_sapiens/" + indata + "?external_db=HGNC";
else
    getURL = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/" + indata + "?external_db=HGNC";
end
headers = {'Content-Type' 'application/xml'};
options = weboptions('HeaderFields', headers, 'Timeout', 5000);
outdata = webread(getURL, options);
if isempty(outdata), return; end

if iscell(outdata)
    outdata = vertcat(outdata{:});
end
outdata = struct2table(outdata);
outdata = convertvars(outdata, 1:width(outdata), @string);
end

%% nearest (coding gene) ==================================================
function outdata = rest_nearest(indata, opts)
if numel(indata) > 1
    error('EnsemblREST::map can only work with single query!')
end

% indata format should be chr:start-end
indata = split(indata, [":", "-", "_"]);

if strcmp(opts.refGenome, '37')
    addurl = "grch37.";
else
    addurl = "";
end
url = @(x) "https://" + addurl + "rest.ensembl.org/overlap/region/human/" + ...
    x + "?feature=gene";
headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', 5000);

bpregion = 1e6;
outdata = [];
while isempty(outdata)
    rstart = double(indata(2)) - bpregion;
    if rstart<= 0, rstart = 1; end
    rstop = double(indata(2)) + bpregion;
    thisurl = url(indata(1) + ":" + rstart + "-" + rstop);
    outdata = webread(thisurl, options);
    bpregion = bpregion*1.2;
    if bpregion > 5e6, break; end
end

try
    outdata = struct2table(vertcat(outdata{:}));
catch % non-matching fields
    cols = ["biotype", "external_name", "id", "gene_id", "start", "end", "description"];
    res = array2table(strings(numel(outdata), 7), ...
        VariableNames=cols);
    res = convertvars(res, ["start", "end"], @double);
    
    for i = 1:numel(outdata)
        for j = 1:numel(cols)
            if isfield(outdata{i}, cols(j))
                res.(cols(j))(i) = outdata{i}.(cols(j));
            end
        end
    end
    outdata = res; clear res
end

% find nearest protein coding gene to this variant
outdata(~strcmp(outdata.biotype, "protein_coding"), :) = [];
outdata.d1 = abs(outdata.start - double(indata(2)));
outdata.d2 = abs(outdata.end - double(indata(2)));
outdata.d = min(outdata.d1, outdata.d2);
[~, idx] = min(outdata.d);
outdata = outdata(idx, :);
if ~isempty(outdata)
    outdata = outdata(:, ["external_name", "id", "gene_id", "start", "end", "description"]);
    outdata = renamevars(outdata, "external_name", "name");
    idx = cellfun(@isempty, outdata.description);
    outdata.description(idx) = {''};
    outdata = convertvars(outdata, ["id", "gene_id", "name", "description"], @string);
end

end

%% map ====================================================================
function outdata = rest_map(indata, opts)
if numel(indata) > 1
    error('EnsemblREST::map can only work with single query!')
end
% indata format should be chr:start-end
if ~any(contains(indata, ["-", "_"])) % exact position (e.g. an SNP)
    endpos = extractAfter(indata, ":");
    indata = indata + "-" + endpos; 
end
indata = replace(indata, ["-", "_"], "..");
% indata = indata + ":1"; %@19NOV2024:removed strandness
opts.refGenome = "GRCh" + string(opts.refGenome);
assemblies = ["GRCh37", "GRCh38"];
indata = opts.refGenome + "/" + indata + "/" + setdiff(assemblies, opts.refGenome);

getURL = "https://rest.ensembl.org/map/human/" + indata;

headers = {'Content-Type' 'application/xml'};
options = weboptions('HeaderFields', headers, 'Timeout', 5000);
outdata = webread(getURL, options);
try
    outdata = outdata.mappings.mapped;
catch % error
end
end

%% eQTL ===================================================================
function outdata = rest_eqtl(indata, opts)

% Ensembl has removed this endpoint and replaced it by a new API
% if strcmp(opts.refGenome, '37')
%     getURL = 'https://grch37.rest.ensembl.org/eqtl/variant_name/homo_sapiens';
% else
%     getURL = 'https://rest.ensembl.org/eqtl/variant_name/homo_sapiens';
% end
% beta = ismember(cellfun(@(x)x.statistic, outdata, 'uni', false), 'beta');
% 
% outdata2 = deal(cellfun(@(x)string({x.value, x.tissue, x.gene}) , outdata(beta), 'uni', false));
% outdata2 = splitvars(table(vertcat(outdata2{:})));
% outdata2.Properties.VariableNames = {'beta', 'tissue', 'gene'};
% 
% outdata1 = deal(cellfun(@(x)string({x.value, x.tissue, x.gene}) , outdata(~beta), 'uni', false));
% outdata1 = splitvars(table(vertcat(outdata1{:})));
% outdata1.Properties.VariableNames = {'p', 'tissue', 'gene'};
% 
% outdata = join(outdata1, outdata2);
% outdata.beta = double(outdata.beta); outdata.p = double(outdata.p);
% outdata = sortrows(outdata, 'p', 'ascend');

if numel(indata) > 1
    error('EnsemblREST::eqtl can only work with single query!')
end
indata = split(indata, ["_", "-", ":"]);
if numel(indata) ~= 4 && numel(indata) ~= 1
    error('EnsemblREST::eqtl query must be in chr$pos$ref$alt [$ is :/-/_] format!')
end

if numel(indata) ~= 1 % rsid
    % variant ID should be in CHR_BP_REF_ALT  format.
    if strcmp(opts.refGenome, '37') % new API only supports GRCh38
        g38 = rest_map(indata(1)+":"+indata(2), opts); % map to GRCh38
        indata(2) = string(g38.start);
    end
    
    indata = "chr" + join(indata, "_");
end

options = weboptions('Timeout', 5000);
studies = webread("https://www.ebi.ac.uk/eqtl/api/studies", options);
studies = string({studies.x_embedded.studies.study_accession}.');
tab = [];
% studies = "GTEx_V8";
for i = 1:numel(studies)
    if opts.verbose
        fprintf('%d (of %d)-study: %s\n', i, numel(studies), studies(i))
    end
    contMe = true; start = 0;
    
    while contMe
    %     getURL = "https://www.ebi.ac.uk/eqtl/api/associations/"+indata+"?study=GTEx_V8&start="+start+"&size=1000";
    %     getURL = "https://www.ebi.ac.uk/eqtl/api/associations/"+indata+"?start="+start+"&size=1000";
        getURL = "https://www.ebi.ac.uk/eqtl/api/associations?variant_id="+indata+"&start="+start+"&size=1000&study="+studies(i);
        headers = {'Content-Type' 'application/xml'};
        options = weboptions('HeaderFields', headers, 'Timeout', 5000);
        try
            outdata = webread(getURL, options);
        catch
            contMe = false;
        end
        try
            outdata = struct2cell(outdata.x_embedded.associations);
            if numel(outdata) < 1000
                contMe = false;
            end
            res = struct2table(vertcat(outdata{:}), 'AsArray', true);
            res(res.pvalue > 1, :) = [];
            if ~isempty(res) && ~isnumeric(res.r2) % unify r2 column if is of type cell
                emptyIdx = cellfun(@isempty, res.r2);
                res.r2(emptyIdx) = {nan};
                res.r2 = double(string(res.r2));
            end
            tab = [tab; res];
            if ~isnumeric(tab.r2)
                tab.r2 = double(string(tab.r2));
            end
        catch
            continue
        end
        start = start + 1e3;
    end
end

try
    outdata = sortrows(tab, 'pvalue', 'ascend');
    varnames = string(outdata.Properties.VariableNames);
    varnames = varnames(ismember(varnames, ["type", "ref", "rsid", ...
        "chromosome", "alt", "variant", "study_id", "qtl_group",...
        "molecular_trait_id", "gene_id", "tissue"]));
    for i = 1:numel(varnames)
        outdata.(varnames(i)) = string(outdata.(varnames(i)));
    end
catch
    outdata = [];
    return
end

if opts.getsymbol && ~isempty(outdata)% get gene IDs
    gdata = ({});
    genes = string(unique(outdata.gene_id));
    grange = unique([1:999:numel(genes), numel(genes) + 1]);
    for i = 1:numel(grange)-1
        gdata{i, 1} = EnsemblREST(genes(grange(i):grange(i+1)-1), 'lookup', 'refGenome', '38');
    end
    gdata = vertcat(gdata{:});
    [idx1, idx2] = ismember(outdata.gene_id, gdata.id);
    outdata.symbol = strings(height(outdata), 1);
    outdata.symbol(idx1) = gdata.name(idx2(idx1));   
end
end

%% Lookup =================================================================
function outdata = rest_lookup(indata, postFlag, refGenome, geneSymbol, parseFlag)
if geneSymbol
    term = 'symbol/homo_sapiens';
else
    term = 'id';
end

if postFlag
    % POST lookup/id
    % Find the species and database for several identifiers. IDs that are not 
    % found are returned with no data.
    if strcmp(refGenome, '37')
        getURL = ['https://grch37.rest.ensembl.org/lookup/', term];
    else
        getURL = ['https://rest.ensembl.org/lookup/', term];
    end
    headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
    options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
        'Timeout', 5000);
    if geneSymbol
        data = struct('symbols',indata);
    else
        data = struct('ids',indata);
    end
    outdata = webwrite(getURL, data, options);
else
    % GET lookup/id/:id
    % Find the species and database for a single identifier e.g. gene,
    % transcript, protein
    if strcmp(refGenome, '37')
        getURL = "https://grch37.rest.ensembl.org/lookup/" + term + "/" + indata + "?expand=1";
    else
        getURL = "https://rest.ensembl.org/lookup/" + term + "/" + indata + "?expand=1";
    end
    headers = {'Content-Type' 'application/xml'};
    options = weboptions('HeaderFields', headers, 'Timeout', 5000);
    outdata = webread(getURL, options);
end

if parseFlag
    if isstruct(outdata)
        if postFlag
            outdata = struct2cell(outdata);
            outdata(cellfun(@isempty, outdata)) = [];
            try
                try
                    outdata = deal(cellfun(@(x)string({x.id, x.biotype,...
                        x.start, x.end, x.seq_region_name, x.display_name}) , outdata, 'uni', false));
                catch % if there is no 'display_name' field for some genes
                    tmp = deal(cellfun(@(x)string({x.id x.biotype,...
                        x.start, x.end, x.seq_region_name}) , outdata, 'uni', false));
                    for i = 1:numel(outdata)
                        if ~isfield(outdata{i}, 'display_name')
                            tmp{i}(6) = "";
                        else
                            tmp{i}(6) = string(outdata{i}.display_name);
                        end
                    end
                    outdata = tmp;
                end
                outdata = splitvars(table(vertcat(outdata{:})));
                outdata.Properties.VariableNames = {'id', 'biotype', 'start', 'end', 'chr', 'name'};
                outdata.start = double(outdata.start);
                outdata.end = double(outdata.end);
            catch
                % just return raw data
            end
        else
            outdata = struct2table(outdata, 'AsArray', true);
            keepcols = {'id', 'display_name', 'biotype', 'start', 'end', 'seq_region_name'};
            [keepidx, idx] = ismember(outdata.Properties.VariableNames, keepcols);
            keepcols = keepcols(idx(keepidx));
            outdata = outdata(:, keepidx);
            keepcols = keepcols(ismember(keepcols, {'id', 'display_name', 'biotype', 'seq_region_name'}));
            outdata = convertvars(outdata, keepcols, "string");
            if ~any(ismember(keepcols, 'display_name'))
                outdata.name = repmat("", height(outdata), 1);
            else
                outdata.Properties.VariableNames = replace(outdata.Properties.VariableNames, 'display_name', 'name');
            end

            if ~any(ismember(keepcols, 'seq_region_name'))
                outdata.chr = repmat("", height(outdata), 1);
            else
                outdata.Properties.VariableNames = replace(outdata.Properties.VariableNames, 'seq_region_name', 'chr');
            end

        end
    end
end
end
%==========================================================================
%% VEP ====================================================================
function outdata = rest_vep(indata, postFlag, refGenome)
dupRegion = false;
if numel(split(indata(1))) > 1 % VCF
    method = "region";
    if ~postFlag % duplicate indata and use POST method, because GET method needs adjustments to current implementation
        postFlag = true;
        indata = [indata; indata];
        dupRegion = true;
    end
else
    method = "id";
end

if postFlag
    % POST vep/:species/id
    % Fetch variant consequences for multiple ids
    if strcmp(refGenome, '37')
        getURL = 'https://grch37.rest.ensembl.org/vep/human/'+method;
    else
        getURL = 'https://rest.ensembl.org/vep/human/'+method;
    end
    headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
    options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
        'Timeout', 5000);
    
    if method == "id"
        data = struct('ids',indata,'LoF',1,'CADD',1,'dbNSFP', ...
            'clinvar_clnsig,clinvar_id,REVEL_score,PrimateAI_pred','domains',1,...
            'protein',1,'hgvs',1,'SpliceAI',1, 'failed', 1, 'merged', 1, ...
            'AlphaMissense', 1);
    elseif method == "region"
        data = struct('variants',indata,'LoF',1,'CADD',1,'dbNSFP',...
            'clinvar_clnsig,clinvar_id,REVEL_score,PrimateAI_pred','domains',1,...
            'protein',1,'hgvs',1,'SpliceAI',1, 'failed', 1, 'merged', 1, ...
            'AlphaMissense', 1);
    end
    try
        outdata = webwrite(getURL, data, options);
    catch ME
        %@22MARCH2025: LoF doesn't work anymore for some reasons with
        %grch37
        data = rmfield(data, "LoF");
        try
            outdata = webwrite(getURL, data, options);
        catch ME
            error(ME.message)
        end
    end
else
    % GET vep/:species/id/:id
    % Fetch variant consequences based on a variant identifier
    if strcmp(refGenome, '37')
        getURL = "https://grch37.rest.ensembl.org/vep/human/"+ method +...
            "/" + indata + "?LoF=1";
    else
        getURL = "https://rest.ensembl.org/vep/human/" + method + "/"...
            + indata + "?LoF=1";
    end
    headers = {'Content-Type' 'application/xml'};
    options = weboptions('HeaderFields', headers, 'Timeout', 5000);
    outdata = webread(getURL, options);
end

try
    if dupRegion; outdata = outdata{1}; end
catch % empty outdata
end

end
%==========================================================================
%% biotoolsfr =============================================================
function outdata = rest_biotoolsfr(indata)
% see: https://www.biotools.fr/human/ensembl_symbol_converter

if numel(indata) > 1 % POST
    getURL = "https://biotools.fr/human/ensembl_symbol_converter/";
    options = weboptions('RequestMethod', 'post', 'Timeout', 5000);
    try
        outdata = webwrite(getURL, 'ids', jsonencode(indata), 'api', 1, options);
    catch ME
        error(ME.message)
    end
    
else % GET
    getURL = "https://biotools.fr/human/ensembl_symbol_converter/?api=1&id=";
    headers = {'Content-Type' 'application/xml'};
    options = weboptions('HeaderFields', headers, 'Timeout', 5000);
    outdata = webread(getURL + indata, options);
end

ids = string(fieldnames(outdata));
names = struct2cell(outdata);
emptyidx = cellfun(@isempty, names);
names(emptyidx) = []; ids(emptyidx) = [];
if ~isempty(names)
    outdata = table(ids, string(names), 'VariableNames', {'id', 'name'});
end

end

%%
function df = aggregateCells(df, sep)
% loops over columns and if columns are cell, merges them using the desired
% delimiter

arguments
    df {mustBeA(df, "table")}
    sep {mustBeTextScalar} = ","
end

cols = colnames(df);
df = df{:, :};
idx = ismissing(df) | cellfun(@isempty, df);
[df{idx}] = deal("");
df = cellfun(@(x)join(string(x), sep), df);
df = array2table(df, VariableNames=cols);
df = standardizeMissing(df, "");

end