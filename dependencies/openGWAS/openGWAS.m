function [res, ginfo] = openGWAS(method, opts)
% implements OpenGWAS API methods: http://gwas-api.mrcieu.ac.uk/docs/

arguments
    method {mustBeTextScalar, mustBeMember(method, ["batches", ...
        "gwasinfo", "associations", "tophits", "phewas", "variants", ...
        "ld-clump", "ld-matrix", "ld-reflookup"])} = "tophits"

    % search query for 'tophits', 'associations' and 'gwasinfo'. Should be
    % used if 'id' is unsued (overrides 'id' option)
    opts.query {mustBeText, mustBeVector}
    opts.id {mustBeText, mustBeVector} % List of GWAS ID(s)
    opts.all (1,1) logical = false % to select all matched strings in 'query' or 'id' fields (non-specific for 'query')
    opts.pval double {mustBeVector, mustBeInRange(opts.pval, 1e-400, 1)} = 5e-8

    % tophits specific
    opts.preclumped (1,1) logical = true % Whether to use pre-clumped tophits
    opts.clump (1,1) logical = true % Whether to clump (1) or not (0)
    opts.bychr (1,1) logical = true % Whether to extract by chromosome (1) or all at once (0). There is a limit on query results so bychr might be required for some well-powered datasets
    opts.r2 (1,1) double {mustBeInRange(opts.r2, 1e-400, 1)} % Clumping parameter ('associations' method: minimum LD r2 for a proxy)
    opts.kb (1,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonzero} = 10000 
    opts.pop {mustBeTextScalar, mustBeMember(opts.pop, ["EUR", "SAS", "EAS", "AFR", "AMR", "legacy"])} = "EUR"
    
    % associations specific
    opts.proxies (1,1) logical = true
    opts.align_alleles (1,1) logical = true
    opts.palindromes (1,1) logical = true
    opts.maf_threshold (1,1) double {mustBeInRange(opts.maf_threshold, 1e-400, 0.5)} = 0.01


    % phewas specific
    opts.variant {mustBeText, mustBeVector} % Comma-separated list of rs IDs, chr:pos or chr:pos range (hg19/b37). e.g rs1205,7:105561135,7:105561135-105563135
    opts.index_list {mustBeText, mustBeVector} % List of study indexes (batch). If empty then searches across all indexes

    % variants specific
    opts.rsid {mustBeText, mustBeVector} % List of variant rs IDs
    opts.chrpos {mustBeText, mustBeVector} % List of variant chr:pos format on build 37 (e.g. 7:105561135)
    opts.radius (1,1) double {mustBeGreaterThanOrEqual(opts.radius, 0)} = 0 % Range to search either side of target locus (for chrpos only)
    opts.gene {mustBeTextScalar} % A gene identifier, either Ensembl or Entrez, e.g. ENSG00000123374 or 1017
    
    % ld specific
    opts.pthresh (1,1) double {mustBeInRange(opts.pthresh, 1e-400, 1)} = 5e-8

    opts.timeout (1,1) double = 5000 % in sec for weboptions 'post'
    opts.verbose (1,1) logical = true
end

url = "http://gwas-api.mrcieu.ac.uk/";
headers = {'Content-Type' 'application/json'};
options.r = weboptions('HeaderFields', headers, 'Timeout', 5000);
headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
options.p = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
    'Timeout', opts.timeout);

% check status ------------------------------------------------------------
if opts.verbose
    raw = webread(url+"status", options.r);
    fprintf('OpenGWAS project:\nAPI version: %s\n', raw.x_API_Version)
    % fprintf('Access: %s\n', raw.Access)
end

% 'X-Api-Token'=access_token,
% 'X-Api-Source'=ifelse(is.null(options()$mrbase.environment), 'R/TwoSampleMR', 'mr-base-shiny')

% if 'query' is used (overrides 'id' option), gwas id(s) are extracted from
% "gwasinfo" and "batches"
if isfield(opts, 'query') && ismember(method, ["tophits", "associations", "gwasinfo"])
    opts = checkMRCfield(opts, "query");
else
    opts.query = [];
end

if ~isempty(opts.query) || ~isfield(opts, 'id')
    opts.id = []; % 'query' overrides 'id'
elseif isfield(opts, 'id')
    opts = checkMRCfield(opts, "id");
end

% check other string fields
opts = checkMRCfield(opts, ["index_list", "variant", "rsid", "chrpos", "gene"]);

if (~isempty(opts.query) && ...
        ismember(method, ["tophits", "associations", "gwasinfo"])) ...
        || nargout > 1 % ginfo
    if opts.verbose, fprintf('getting gwasinfo...'); end
    ginfo = webread(url + "gwasinfo", options.r);
    batches = webread(url + "batches", options.r);
    ginfo = getgwastab(ginfo, batches); % make a table for ginfo
    clear batches
    if opts.verbose, fprintf("\b\b done.\n"); end
else
    ginfo = [];
end

% initialize r2 field
if ~isfield(opts, 'r2')
    if method == "associations"
        opts.r2 = 0.8;
    else % tophits
        opts.r2 = 0.001;
    end
end

if method == "gwasinfo" % =================================================
    if isempty(ginfo)
        if isempty(opts.id) % fetch all gwas ids
            ginfo = webread(url + method, options.r);
            ginfo = getgwastab(ginfo, []);
        else
            if numel(opts.id) == 1 % GET
                ginfo = webread(url + method + "/" + opts.id, options.r);
            else % POST
                data = getmrcstruct(opts, "id");
                ginfo = webwrite(url + method, data, options.p);
            end
            ginfo = struct2table(ginfo);
            try
                ginfo = ginfo(:, {'trait','year','sample_size','population','consortium','id'});
            catch
                ginfo = ginfo(:, {'trait','year','population','consortium','id'});
                ginfo.sample_size = nan(height(ginfo), 1);
                ginfo = movevars(ginfo, 'sample_size', 'after', 'year');
            end
        end
        ginfo = convertvars(ginfo, 1:width(ginfo), "string");
        ginfo.sample_size = double(ginfo.sample_size);
    end
    
    if ~isempty(opts.query) || ~isempty(opts.id)
        [~, ginfo] = findIdFromQuery(ginfo, opts, nargout);
    end
    res = ginfo;
    return
  
elseif method == "batches" % ==============================================
    if isempty(ginfo)
        ginfo = webread(url + method, options.r);
        ginfo = getgwastab([], ginfo);
    end
    res = ginfo;
    return

elseif method == "associations" % =========================================
    [~, ginfo, opts] = findIdFromQuery(ginfo, opts, nargout);
    
    if isempty(opts.id) || isempty(opts.variant)
        error('openGWAS:associations', 'both id and variant should be provided!')
    end
    data = getmrcstruct(opts, ["id", "variant", "proxies", "r2", ...
        "align_alleles", "palindromes", "maf_threshold"]);
    res = webwrite(url + method, data, options.p);

elseif method == "tophits" % ==============================================

    % if both 'query' and 'id' are empty, gwas summary for all studies will
    % be read (not recommended)
    [~, ginfo, opts] = findIdFromQuery(ginfo, opts, nargout);

    data = getmrcstruct(opts, ["id", "pval", "preclumped", "clump", ...
        "r2", "kb", "pop", "bychr"]);
    res = webwrite(url + method, data, options.p);
    
elseif method == "phewas" % ===============================================
    data = getmrcstruct(opts, ["variant", "pval", "index_list"]);
    res = webwrite(url + method, data, options.p);

elseif method == "variants" % =============================================
    if isempty(opts.gene)
        data = getmrcstruct(opts, ["rsid", "chrpos", "radius"]);
        res = webwrite(url + method + "/afl2", data, options.p);
    else
        options.r.HeaderFields = []; % donno why
        try
            res = webread(url + method + "/gene/" + opts.gene + "?radius=" + opts.radius, options.r);
            if ~res
                ens = EnsemblREST(opts.gene, "lookup", "geneSymbol", true, 'refGenome', '37');
                res = webread(url + method + "/gene/" + ens.id + "?radius=" + opts.radius, options.r);
            end
        catch % gene symbol instead of gene ID?
            ens = EnsemblREST(opts.gene, "lookup", "geneSymbol", true, 'refGenome', '37');
            res = webread(url + method + "/gene/" + ens.id + "?radius=" + opts.radius, options.r);
        end
    end

elseif method == "ld-clump" % =============================================
    if isempty(opts.rsid)
        error('openGWAS:ldClump', 'rsid cannot be empty!')
    end

    if numel(opts.rsid) ~= numel(opts.pval)
        error('openGWAS:ldClump', 'rsid and pval must be of the same length!')
    end
    data = getmrcstruct(opts, ["rsid", "pval", "pthresh", "r2", "kb", "pop"]);
    res = webwrite(url + replace(method, "-", "/"), data, options.p);

elseif method == "ld-matrix" % ============================================
    data = getmrcstruct(opts, ["rsid", "pop"]);
    res = webwrite(url + replace(method, "-", "/"), data, options.p);
    res.matrix = horzcat(res.matrix{:});
    res.matrix = double(string(res.matrix)).^2;
    
    underlineIdx = cellfun(@(x)x(end-1), regexp(res.snplist, '_'));
    res.snp = cellfun(@(x,y)x(1:y-1), res.snplist, num2cell(underlineIdx), 'UniformOutput', false);
    if ~all(ismember(res.snp, data.rsid))
        error('openGWAS:failedLdMatrix', 'failed to regenerate rsids from API!')
    end
    
    [idx1, idx2] = ismember(data.rsid, res.snp); idx2(idx2<1) = [];
    res.snplist = res.snplist(idx2); res.snp = res.snp(idx2);
    res.matrix = res.matrix(:, idx2);
    res.matrix = res.matrix(idx2, :);
    
    newmat = spalloc(numel(data.rsid), numel(data.rsid), 0);
    newmat(idx1, idx1) = res.matrix;
    res = rmfield(res, {'matrix', 'snplist'});
    res.ld = full(newmat); clear newmat
    res.snp = string(res.snp);
    return

elseif method == "ld-reflookup" % =========================================
    data = getmrcstruct(opts, ["rsid", "pop"]);
    res = webwrite(url + replace(method, "-", "/"), data, options.p);

end


% format results into a table
if iscell(res)
    try
        res = vertcat(res{:});
    catch % structs with different field names
        fi = cellfun(@fieldnames, res, 'uni', false);
        fi = string(unique(vertcat(fi{:})));
        r = res; res = struct;
        for i = 1:numel(r)
            for j = 1:numel(fi)
                if isfield(r{i}, fi(j))
                    res(i).(fi(j)) = r{i}.(fi(j));
                else
                    res(i).(fi(j)) = missing;
                end
            end
        end
    end
end

if ~isempty(res)
    if isstruct(res)
        if size(res, 1) > 1
            res = struct2table(res);
        else
            res = struct2table(res, 'AsArray', true);
        end
    else
        res = string(res);
    end
end

end % END

%% subfunctions ===========================================================
function gwasinfo = getgwastab(gwasinfo, batches)

if ~isempty(gwasinfo)
    fnames = string(fieldnames(gwasinfo));
end

if ~isempty(batches)
    batches = [string(cellfun(@(x)x.id, batches, 'UniformOutput', false)),...
    string(cellfun(@(x)x.description, batches, 'UniformOutput', false))];

    if ~isempty(gwasinfo)
        batchesNew = strings(numel(fnames), 1);
        for i = 1:size(batches, 1)
            idx = contains(fnames, batches(i, 1));
            batchesNew(idx, 1) = batches(i, 1);
            batchesNew(idx, 2) = batches(i, 2);
        end
    else
        batchesNew = batches;
    end
end

if isempty(gwasinfo)
    gwasinfo = batchesNew;
    gwasinfo = array2table(gwasinfo, 'VariableNames', {'batch_id',...
        'batch_name'});
    gwasinfo.batch_id = replace(gwasinfo.batch_id, '_', '-');
elseif isempty(batches)
    gwasinfo = [gimmeFileds(gwasinfo, 'trait'), ...
        gimmeFileds(gwasinfo, 'year'), gimmeFileds(gwasinfo, 'sample_size'),...
        gimmeFileds(gwasinfo, 'population'), ...
        gimmeFileds(gwasinfo, 'consortium'), gimmeFileds(gwasinfo, 'id')];
    gwasinfo = array2table(gwasinfo, 'VariableNames', ...
        {'trait','year','sample_size','population','consortium','id'});
    gwasinfo.year = double(gwasinfo.year);
    gwasinfo.sample_size = double(gwasinfo.sample_size);
else
    gwasinfo = [gimmeFileds(gwasinfo, 'trait'), ...
        gimmeFileds(gwasinfo, 'year'), gimmeFileds(gwasinfo, 'sample_size'),...
        gimmeFileds(gwasinfo, 'population'), ...
        gimmeFileds(gwasinfo, 'consortium'), gimmeFileds(gwasinfo, 'id'), batchesNew];
    gwasinfo = array2table(gwasinfo, 'VariableNames', ...
        {'trait','year','sample_size','population','consortium','id','batch_id',...
        'batch_name'});
    gwasinfo.batch_id = replace(gwasinfo.batch_id, '_', '-');
    gwasinfo.year = double(gwasinfo.year);
    gwasinfo.sample_size = double(gwasinfo.sample_size);
end

end

%% ------------------------------------------------------------------------
function outstruct = gimmeFileds(instruct, infield)
try
    outstruct = string(struct2cell(structfun(@(x) x.(infield), instruct, ...
        'UniformOutput', false)));
catch
    outstruct = strings(numel(instruct), 1);
    fnames = fieldnames(instruct);
    for i = 1:numel(fnames)
        try
            outstruct(i, 1) = instruct.(fnames{i}).(infield);
        catch
            outstruct(i, 1) = "-";
        end
    end
end
end

%% ------------------------------------------------------------------------
function data = getmrcstruct(opts, fi)

data = struct;
fi = string(fi);

for i = 1:numel(fi)
    if isfield(opts, fi(i)), data.(fi(i)) = opts.(fi(i)); end
end
end

%% ------------------------------------------------------------------------
function [res, ginfo, opts] = findIdFromQuery(ginfo, opts, nout)

res = [];
if ~isempty(opts.query)
    fidx = contains(lower(ginfo.trait), opts.query.lower);
    if ~any(fidx)
        fidx = ismember(ginfo.id, opts.query);
    end
elseif isfield(opts, 'id') && (~isempty(ginfo) || nout > 1)
    fidx = ismember(ginfo.id, opts.id);
end

if ~isempty(ginfo)
    if ~any(fidx)
        fprintf('nothing found!\n')
        ginfo = [];
        return
    end

    ginfo = ginfo(fidx, :);
    ginfo = sortrows(ginfo, 'id');
    ginfo.index = (1:size(ginfo,1))';
    ginfo = movevars(ginfo, 'index', 'before', 'trait');
    ginfo = movevars(ginfo, 'id', 'after', 'index');
    ginfo = movevars(ginfo, 'year', 'after', 'population');

    if any(ismember(ginfo.Properties.VariableNames, "batch_id"))
        ginfo = movevars(ginfo, 'batch_id', 'after', 'trait');
    end
    
    if opts.verbose, disp(ginfo); end
    
    if height(ginfo) > 1 % multiple matches
        if ~opts.all
            format("compact")
            if ~opts.verbose, disp(ginfo); end % show the table only once 
            sel = input("index(indices) [enter to select all]? ");
            if isempty(sel) || any(sel > height(ginfo)) || any(sel < 0)
                sel = 1:height(ginfo);
            end
            format("default")
        else
            sel = 1:height(ginfo);
        end

        ginfo = ginfo(sel, :);
    end

    opts.id = ginfo.id;
end

end

%% ------------------------------------------------------------------------
function opts = checkMRCfield(opts, fi)
fi = string(fi);

for i = 1:numel(fi)
    if isfield(opts, fi(i))
        opts.(fi(i)) = string(opts.(fi(i)));
        rmidx = opts.(fi(i)) == "" | ismissing(opts.(fi(i)));
        opts.(fi(i))(rmidx) = [];
    else
        opts.(fi(i)) = [];
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [tophits, foundkeys] = openGWAS(keyword, method, p)
% % http://gwas-api.mrcieu.ac.uk/docs/
% 
% [tophits, foundkeys] = deal([]);
% 
% if nargin < 3
%     % only allowed to be empty with "tophits"; it must be a struct when
%     % method = "phewas" with at least one field named 'variant'. It can
%     % also be a string array or struct (with field: 'id') in case method is
%     % set to "gwasinfo".
%     p = [];
% end
% if nargin < 2
%     method = "tophits";
% else
%     method = string(method);
%     if ~isscalar(method)
%         fprintf('ERROR: method must be a scalar string!\n')
%         return
%     elseif ~any(ismember(method, ["tophits", "phewas", "gwasinfo"]))
%         fprintf('ERROR: method can only be "tophits", "phewas" or "gwasinfo"!\n')
%         return
%     elseif isempty(p) && method == "phewas"
%         fprintf('ERROR: p struct must not be empty!\n')
%         fprintf(['it must containt field variant with rs IDs, chr:pos',...
%             ' or chr:pos range',...
%             ' (hg19/b37).\n e.g rs1205,7:105561135,7:105561135-105563135\n'])
%         return
%     elseif ~isfield(p, 'variant') && method == "phewas"
%         fprintf('ERROR: field variant is missing from p!\n') 
%         fprintf(['it must containt rs IDs, chr:pos or chr:pos range',...
%             ' (hg19/b37).\n e.g rs1205,7:105561135,7:105561135-105563135\n'])
%         return
%     elseif isempty(p) && method == "gwasinfo"
%         fprintf('ERROR: p must not be empty!\n')
%         fprintf(['it can be a struct (with field ''id'')\n',...
%             'or a string array containing list of GWAS IDs',...
%             '.\n e.g ["prot-a-235", "ebi-a-GCST006804"]\n'])
%         return
%     end
%     
%     % set default field values
%     if method == "phewas"
%         if ~isfield(p, 'pval') && method == "phewas"
%             p.pval = 1e-3;
%         end
%         if ~isfield(p, 'index_list') && method == "phewas"
%             p.index_list = [];
%         end
%     end
% end
% if nargin < 1 
%     keyword = [];
% end
% fprintf('\n********************* openGWAS *************************\n')
% % implements OpenGWAS API methods: http://gwas-api.mrcieu.ac.uk/docs/
% getURL = "http://gwas-api.mrcieu.ac.uk/";
% headers = {'Content-Type' 'application/json'};
% options = weboptions('HeaderFields', headers, 'Timeout', 5000);
% home = fileparts(which('openGWAS.m'));
% try
%     gwasinfo = load(fullfile(home, 'gwasinfo.mat'));
%     gwasinfo = gwasinfo.gwasinfo;
% catch
%     fprintf('getting gwasinfo, this takes few minutes...\n')
%     gwasinfo = webread(getURL+"gwasinfo", options);
%     batches = webread(getURL+"batches", options);
%     gwasinfo = getgwastab(gwasinfo, batches); % make a table for gwasinfo
%     save(fullfile(home, 'gwasinfo.mat'), 'gwasinfo')
% end
% 
% % check status ------------------------------------------------------------
% raw = webread(getURL+"status", options);
% fprintf('OpenGWAS project:\nAPI version: %s\n', raw.APIVersion)
% fprintf('Access: %s\n', raw.Access)
% 
% headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
% % 'X-Api-Token'=access_token,
% % 'X-Api-Source'=ifelse(is.null(options()$mrbase.environment), 'R/TwoSampleMR', 'mr-base-shiny')
% options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
%     'Timeout', 5000);
% 
% if method == "tophits"
%     % search within traits
%     searchFlag = true;
%     while searchFlag
%         if isempty(keyword)
%             keyword = input('search trait keyword (case insensitive): ', 's');
%         end
%         foundkeysIDX = contains(lower(gwasinfo.trait), keyword.lower);
%         if ~any(foundkeysIDX)
%             foundkeysIDX = ismember(gwasinfo.id, keyword);
%         end
%         
%         if ~any(foundkeysIDX)
%             contAns = input('nothing found :( wanna try again?[Y/n] ', 's');
%             if strcmpi(contAns, 'n')
%                 searchFlag = false;
%             end
%             keyword = [];
%         else
%             foundkeys = (gwasinfo(foundkeysIDX, :));
%             foundkeys = sortrows(foundkeys, 'id');
%             foundkeys.num = (1:size(foundkeys,1))';
%             foundkeys = movevars(foundkeys, 'num', 'before', 'trait');
%             foundkeys = movevars(foundkeys, 'year', 'after', 'population');
%             foundkeys = movevars(foundkeys, 'batch_id', 'after', 'trait');
% 
%             disp(foundkeys);
% %             if size(foundkeys, 1) == 1
% %                 selKey = 1;
% %             else
% %                 selKey = [];  
% %                 while isempty(selKey) || numel(selKey) > 1 || ...
% %                     selKey > numel(foundkeys) || selKey < 1
% %                     selKey = input('? ');
% %                 end
% %             end
% %             foundkeys = foundkeys(selKey, :);
%             searchFlag = false;
%         end
%     end
%     fprintf('***************************************************\n')
%     % get tophits -------------------------------------------------------------
%     data = struct('id', foundkeys.id, 'pval', 5e-8,...
%         'preclumped', 0, 'clump', 1, ...
%         'r2', 0.001, 'kb', 10000, 'pop', 'EUR');
%     tophits = webwrite(getURL+"tophits", data, options);
%     if iscell(tophits)
%         tophits = [tophits{:}];
%     end
%     if ~isempty(tophits)
%         if size(tophits, 1) > 1
%             tophits = struct2table(tophits);
%         else
%             tophits = struct2table(tophits, 'AsArray', true);
%         end
%     end
%     
% elseif method == "phewas" % ===============================================
%     data = struct('variant', p.variant, 'pval', p.pval,...
%         'index_list', p.index_list);
%     try
%         tophits = webwrite(getURL+"phewas", data, options);
%     catch ME
%         fprintf('openGWAS error: %s\n', ME.message)
%         return
%     end
%     if size(tophits, 1) > 1
%         tophits = struct2table(horzcat(tophits{:}));
%     else
%         tophits = struct2table(horzcat(tophits{:}), 'AsArray', true);
%     end
%     
% elseif method == "gwasinfo" % =============================================
%     if isstruct(p)
%         data = struct('id', p.id);
%     else
%         data = struct('id', p);
%     end
%     
%     try
%         tophits = webwrite(getURL+"gwasinfo", data, options);
%     catch ME
%         fprintf('openGWAS error: %s\n', ME.message)
%         return
%     end
% %     if size(tophits, 1) > 1
% %         tophits = struct2table(horzcat(tophits{:}));
% %     else
% %         tophits = struct2table(horzcat(tophits{:}), 'AsArray', true);
% %     end
% end
% 
% end % END
% 
% %% subfunctions ===========================================================
% function gwasinfo = getgwastab(gwasinfo, batches)
% fnames = string(fieldnames(gwasinfo));
% batches = [matlab.lang.makeValidName(string(cellfun(@(x)x.id, batches,...
%     'UniformOutput', false))),...
% string(cellfun(@(x)x.description, batches, 'UniformOutput', false))];
% batchesNew = strings(numel(fnames), 1);
% for i = 1:size(batches, 1)
%     idx = contains(fnames, batches(i, 1));
%     batchesNew(idx, 1) = batches(i, 1);
%     batchesNew(idx, 2) = batches(i, 2);
% end
% gwasinfo = [gimmeFileds(gwasinfo, 'trait'), ...
%     gimmeFileds(gwasinfo, 'year'), gimmeFileds(gwasinfo, 'sample_size'),...
%     gimmeFileds(gwasinfo, 'population'), ...
%     gimmeFileds(gwasinfo, 'consortium'), fnames, batchesNew];
% gwasinfo = array2table(gwasinfo, 'VariableNames', ...
%     {'trait','year','sample_size','population','consortium','id','batch_id',...
%     'batch_name'});
% gwasinfo.id = replace(gwasinfo.id, gwasinfo.batch_id+"_",...
%     replace(gwasinfo.batch_id, '_', '-')+"-");
% gwasinfo.batch_id = replace(gwasinfo.batch_id, '_', '-');
% gwasinfo.year = double(gwasinfo.year);
% gwasinfo.sample_size = double(gwasinfo.sample_size);
% end
% 
% %% 
% function outstruct = gimmeFileds(instruct, infield)
% try
%     outstruct = string(struct2cell(structfun(@(x) x.(infield), instruct, ...
%         'UniformOutput', false)));
% catch
%     outstruct = strings(numel(instruct), 1);
%     fnames = fieldnames(instruct);
%     for i = 1:numel(fnames)
%         try
%             outstruct(i, 1) = instruct.(fnames{i}).(infield);
%         catch
%             outstruct(i, 1) = "-";
%         end
%     end
% end
% end