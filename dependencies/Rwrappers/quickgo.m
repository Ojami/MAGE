function res = quickgo(goId, opts)
% connects and fetches data from QuickGO API.
% Under development
% Jan 2023, Oveis Jamialahmadi, University of Gothenburg

arguments
    goId {mustBeText} % list of GO terms should be in format of GO:XXXX
    opts.goUsage {mustBeMember(opts.goUsage, ["exact", "descendants", "slim"])} = "exact"
    opts.taxonId (1,1) double = 9606 % homo sapiens
    opts.goUsageRelationships {mustBeMember(opts.goUsageRelationships, ["is_a", "part_of", "occurs_in", "regulates"])} = ["is_a"]
    opts.geneProductType {mustBeMember(opts.geneProductType, ["miRNA", "complex", "protein"])} = "protein"
    opts.includeFields {mustBeMember(opts.includeFields, ["goName", "taxonName" , "name", "synonyms"])} = ["goName" , "name", "synonyms"]
    opts.limit (1,1) double {mustBeInRange(opts.limit, 1, 200)} = 200 % results per page
    opts.page (1,1) double {mustBeInRange(opts.page, 1, 10000)} = 1 % starts with 1 page and fetches all the rest (do not change)
end

tmp = innerGoSearch(goId, opts);

pageinfo = tmp.pageInfo;
res{pageinfo.current} = tmp;

% check how many pages left
if (pageinfo.current < pageinfo.total)
    for k = (pageinfo.current+1):pageinfo.total
        ropts = opts;
        ropts.page = k;
        res{k} = innerGoSearch(goId, ropts);
    end
end

res = cellfun(@(x)x.results, res, uni=false);
res = vertcat(res{:});
res = struct2table(res);
cols = colnames(res);
for k = 1:numel(cols)
    tmp = res.(cols(k));
    if isnumeric(tmp)
        continue
    end
    idx = cellfun(@isempty, tmp);
    if all(idx)
        tmp(:) = {missing};
    else
        tmp(idx) = {''};
    end

    try
        res.(cols(k)) = string(tmp);
    catch
        res.(cols(k)) = tmp;
    end
end
res = rmmissing(res, 2, "MinNumMissing", height(res));

end % END

%% subfunctions ============================================================
function tmp = innerGoSearch(goId, opts)

url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?";

fis = string(fieldnames(opts));
terms = strings(numel(fis), 1);
for k = 1:numel(fis)
    terms(k) = fis(k) + "=" + join(string(opts.(fis(k))), ",");
end
terms = join(terms, "&");

headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', 5e4);
tmp = webread(url + "goId=" + goId.join(",") + "&" + terms, options);

end