function R = opentargets_genetics(query, p)
% implements GraphQL queries to access genetic.opentargets server.
% https://api.genetics.opentargets.org/graphql/browser
% Oveis Jamialahmadi, University of Gothenburg, Dec 2021.
% 
% @23/12/2021: a provisional modification to 'parseType' subfunction was
%              made to avoid error when fields of each data type have
%              specific arguments as well. Those fields are removed
%              currently.
% @29APR2025: opentargets genetics is now retired.

arguments
    query {mustBeA(query, ["string", "struct"])} % string when there is 1 argument, otherwise should a struct with fieladnames as arguments.
    p.method {mustBeMember(p.method, ["search", "genes", "geneInfo",...
        "variantInfo", "studiesForGene", "manhattan", "tagVariantsAndStudiesForIndexVariant",...
        "indexVariantsAndStudiesForTagVariant", "pheWAS", "genesForVariant",...
        "colocalisationsForGene", "gwasColocalisationForRegion",...
        "studiesAndLeadVariantsForGene", "gecko"])} = "search"
end

getURL = "https://api.genetics.opentargets.org/graphql";
headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', 20000);

p = parseSchema(p); 
data = getgraphql(query, p);
try
    R = webwrite(getURL, data, options);
catch ME
    R = ME.message;
end

end % END

%% subfunctions ===========================================================
function data = getgraphql(query, p)
args = p.qr.args.str(p.qr.args.need);
if isstruct(query)
    fnames = string(fieldnames(query));
    if any(~ismember(args, fnames))
        error("opentargets:query fields must be consistent with the method arguments:%s", join(args, ","))
    end
    fnames = setdiff(fnames, args);
    if ~isempty(fnames); query = rmfield(query, fnames); end
    instruct = query;
else
    if numel(args) > 1
        error("opentargets:query must be a struct since there are multiple arguments:\n%s", join(args, ","))
    end
    instruct.(args) = query;
end

data = struct('query', p.cmd, 'variables', instruct);

end

%% ------------------------------------------------------------------------
function p = parseSchema(p)
file = "genetics.opentargets.schema.mat";
pth = fullfile(fileparts(which('opentargets.m')), file);
if ~exist(pth, 'file')
    schema = string(webread("https://api.genetics.opentargets.org/graphql/schema"));
    save(pth, 'schema');
else
    schema = load(pth).schema;
end

raw = extract(schema, p.method + "(" + wildcardPattern + ")" + wildcardPattern + lineBoundary('end')).strip;
if isempty(raw)
    error('opentargets::this method cannot be found in schema!')
end
raw = splitlines(raw).strip;
raw(1) = []; % contains method header
[raw, qr.info] = getInfo(raw);

if numel(raw) > 1; raw = raw.join; end
% method's args 
[qr.args.str, qr.args.type, qr.args.need] = sepStrType(raw);

% method's type
qr.type.name = extractAfter(raw, '):').strip;
qr.type.name = erase(qr.type.name, '!' + textBoundary("end"));
qr.type.name = replace(qr.type.name, ["[", "!", "]"], '');
qr.type.fields = parseType(schema, qr.type.name);

% sub-types (and sub-fields) for the main method's type
qr.type.tree = struct;
qr.type.tree = str2struct(qr.type.fields.str', qr.type.fields.type.');
qr.type.tree = addSubfields(schema, qr.type.tree);
qr.method = p.method;
p.cmd = join(constructQuery(qr), newline);
p.qr = qr;
end

%% ------------------------------------------------------------------------
function [str, type, need] = sepStrType(raw)
if numel(raw) > 1; raw = raw.join; end
str = split(extractBefore(raw, ')'), ',').strip; % arguments
need = endsWith(str, '!'); % mandatory arguments in the query
type = extractBetween(str, ': ', textBoundary('end'));
type = erase(type, '!' + textBoundary("end"));
str = erase(str, ':' + wildcardPattern + textBoundary);
end

%% ------------------------------------------------------------------------
function [raw, info] = getInfo(raw)
raw(raw == "") = [];
info_idx = contains(raw, '"' + wildcardPattern + '"');
info = raw(info_idx);  % description of the method
info = replace(info, '"', '');
raw(info_idx) = [];
end

%% ------------------------------------------------------------------------
function tab = parseType(schema, name)
info = extractBetween(schema, lineBoundary("start") + "type " + name + " {", lineBoundary("start") + "}").strip; % type fields
if isempty(info)
    tab = [];
    return
end

info = splitlines(info).strip;
% remove those fields which have themselves arguments (e.g.
% topOverlappedStudies field for Manhattan type)
checkPar = contains(info, "(") & count(info, '"') < 2;
if any(checkPar)
    checkinfo = info;
    remidx = find(contains(checkinfo, ["(", ")"]));
    remidx = reshape(remidx, 2, []);
    for i = 1:size(remidx, 2)
        if count(checkinfo(remidx(1, i) - 1), '"') >= 2 % description line
            remidx(1, i) = remidx(1, i) - 1;
        end
        checkinfo(remidx(1, i):remidx(2, i)) = missing;
    end
    checkinfo(ismissing(checkinfo) | checkinfo == "") = [];
    info = checkinfo;
   
end

[str, info] = getInfo(info);
type = extractBetween(str, ': ', textBoundary('end'));
type = replace(type, ["[", "]", "!"], '');
str = erase(str, ':' + wildcardPattern + textBoundary);
if isempty(info); info = strings(numel(str), 1); end
tab = struct('str', str, 'type', type, 'info', info);
end

%% ------------------------------------------------------------------------
function qr = addSubfields(schema, qr)
types = string(fieldnames(qr));

for i = 1:numel(types)
    tt = parseType(schema, qr.(types(i)));
    if ~isempty(tt)
        qr.(types(i)) = str2struct(tt.str, tt.type);
        qr.(types(i)) = addSubfields(schema, qr.(types(i)));
    end
end
end

%% ------------------------------------------------------------------------
function str = str2struct(fields, vals)
for i = 1:numel(fields)
    str.(fields(i)) = vals(i);
end
end

%% ------------------------------------------------------------------------
function cmd = constructQuery(qr)
% constructs the query for open targets graphql API.

% query header: ignore optional arguments
cmd(1) = "query target("; cmd(2) = blanks(2) + string(qr.method) + "(";
qr.args.str(~qr.args.need) = []; qr.args.type(~qr.args.need) = [];
for i = 1:numel(qr.args.str)
    cmd(1) = cmd(1) + "$" +  qr.args.str(i) + ": " + qr.args.type(i);
    cmd(1) = cmd(1) + "!"; 
    cmd(2) = cmd(2) + qr.args.str(i) + ": $" + qr.args.str(i); 
    if i < numel(qr.args.str)
        cmd(:) = cmd(:) + ", ";
    else
        cmd(:) = cmd(:) + ") {";
    end
end

cmd = [cmd.'; constructQueryBody(qr); blanks(2) + "}"; "}"];

end

%% ------------------------------------------------------------------------
function body = constructQueryBody(qr, indent)

if nargin < 2
    str = qr.type.tree;
    body = string;
    indent = 4;
else
    str = qr;
end
cnt = 1;

fnames = string(fieldnames(str));
for i = 1:numel(fnames)
    body(cnt, 1) = blanks(indent) + fnames(i);
    if isstruct(str.(fnames(i)))
        body(cnt, 1) = body(cnt, 1) + " {";
        bodyInner = constructQueryBody(str.(fnames(i)), indent + 2);
        body = [body; bodyInner];
        body(end) = body(end) + "}";
        cnt = numel(body);
    end
    cnt = cnt + 1;
end

end