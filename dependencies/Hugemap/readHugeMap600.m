function hits = readHugeMap600(query, opts)
% 24 JUN 2025

arguments
    query {mustBeTextScalar}
    opts.index {mustBeMember(opts.index, ["gene-associations-600trait", "gene-finder-600trait"])} = "gene-associations-600trait"
    opts.ancestry = "Mixed"
    opts.cohort = "UKB_450k_AoU_250k_MGB_53k_META_overlapcorrected"
end

% set-up --------------------------------------------------------
base = "https://bioindex.hugeamp.org/api/bio/";


% build raw q= string (keep commas literal!)
qraw = opts.ancestry + "," + opts.cohort + "," + query;

url  = base + "query/" + opts.index + "?q=" + qraw;   % no urlencode needed
wopts = weboptions(Timeout = 60);

% pull first page ----------------------------------------------
out  = webread(url, wopts);
hits = out.data;
cont = string(out.continuation);    % token for next page

% follow continuation tokens -----------------------------------
while ~ismissing(cont)
    page = webread(base + "cont?token=" + cont, wopts);
    hits = [hits ; page.data];      %#ok<AGROW>
    cont = string(page.continuation);
end

if iscell(hits)
    hits = parseHits(hits);
else
    hits = struct2table(hits);
end

isCellOfText = @(v)  iscell(v) && ...
    all(cellfun(@(c) ischar(c) || isstring(c), v));
hits = convertvars(hits, isCellOfText, "string");

end % END

%% subfunctions ===========================================================
function hits = parseHits(hits)

% flatten the cell so every element is a scalar struct 
hits = cellfun(@(s) num2cell(s(:).'), hits, uni=false);
hits = [hits{:}];

%   num2cell(s(:).') turns a K×1 struct array into a 1×K cell of scalar structs
%   The horizontal concat [ …{:} ] gives a 1×M cell, M = total structs

% build the master field list
allFields = cellfun(@(x)string(fieldnames(x)), hits, uni=false);
allFields = unique(vertcat(allFields{:}));

% prototype with sensible empties 
proto = struct();
for f = allFields'
    idx     = find(cellfun(@(s) isfield(s,f), hits), 1);
    example = hits{idx}.(f);

    if  ischar(example) || isstring(example)
        proto.(f) = "";
    elseif isnumeric(example) || islogical(example)
        proto.(f) = NaN;
    else
        proto.(f) = [];
    end
end

% pad each struct & align field order 
pad = @(s) orderfields(mergeInto(proto,s), proto);
hits = cellfun(pad, hits, uni=false);

% smash together & convert 
hits = vertcat(hits{:});    
hits = struct2table(hits);

end % END

%% ------------------------------------------------------------------------
function out = mergeInto(base, extra)

out = base;
fn  = fieldnames(extra);
for k = 1:numel(fn)
    out.(fn{k}) = extra.(fn{k});
end

end % END