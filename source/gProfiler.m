function out = gProfiler(query, opts)
%@31JAN2024 An API for g:Profiler convert method.

arguments
    query {mustBeText, mustBeVector}
    opts.target {mustBeMember(opts.target, ["ENTREZGENE_ACC", "ENSG", ...
        "HGNC"])} = "ENSG"
    opts.organism {mustBeTextScalar} = "hsapiens"
    opts.Timeout (1,1) double = 5e3
end

headers = {'Content-Type' 'application/json'};
options = weboptions('HeaderFields', headers, 'Timeout', opts.Timeout);
getURL = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/";

data = struct("query", string(query), ...
    "target", opts.target, ...
    "organism", opts.organism, ...
    "numeric_ns", "ENTREZGENE_ACC");
out = webwrite(getURL, data, options);
out = jsondecode(out);
out = struct2table(out.result, AsArray=true);
out = convertvars(out, 1:width(out), @string);

end %END