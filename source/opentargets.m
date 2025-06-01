function R = opentargets(query, p)
% OPENTARGETS  Lightweight GraphQL client for the Open Targets platform.
%   R = OPENTARGETS(QUERY) submits the GraphQL QUERY (string or struct of
%   variables) to https://api.platform.opentargets.org/api/v4/graphql and
%   returns the decoded JSON response.
%
%   R = OPENTARGETS(QUERY, P) allows additional parameters:
%     P.method : name of the root Query method (e.g. "target", "disease").
%                If omitted the client tries to guess it (see below).
%     P.help   : true/false (default false).  When true OPENTARGETS prints a
%                list of available root‑level methods (name + arg list) and
%                returns without issuing a web request.
%
% The function keeps the original structure of the user‑supplied version
% wherever no change was required; only small patches were added to support
% the new features (schema caching, help flag, automatic method guess).
%
% Oveis Jamialahmadi · University of Gothenburg · first release  Apr 2025
% This patch: 29 Apr 2025

arguments
    query {mustBeA(query,["string","struct"])}
    p.method {mustBeMember(p.method, ["target", "targets", "disease", ...
        "diseases", "drug", "drugs", "search", "facets", "mapIds", ...
        "geneOntologyTerms", "variant", "study", "studies", ...
        "credibleSet", "credibleSets"])} = "search"
    p.help   logical = false
end

url     = "https://api.platform.opentargets.org/api/v4/graphql";
headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
options = weboptions('RequestMethod','post','HeaderFields',headers,'Timeout',20000);

% -----------------------------------------------------------------------
% 0.  Obtain (cached) SDL schema ----------------------------------------
% -----------------------------------------------------------------------
schema = getSchemaCached(url);

% -----------------------------------------------------------------------
% 0a. Help flag: list methods & exit -------------------------------------
% -----------------------------------------------------------------------
if p.help
    R = listAvailableMethods(schema);
    return
end

% -----------------------------------------------------------------------
% 1.  Deduce method name if not supplied ---------------------------------
% -----------------------------------------------------------------------
if strlength(p.method)==0
    p.method = guessMethod(query,schema);
end

% -----------------------------------------------------------------------
% 2.  Parse schema for that method (no behaviour change otherwise) --------
% -----------------------------------------------------------------------

p = parseSchema(p,schema);   % *** only signature was changed ***

% -----------------------------------------------------------------------
% 3.  Build request JSON & POST ------------------------------------------
% -----------------------------------------------------------------------
body = getgraphql(query,p);
try
    R = webwrite(url, body, options);
catch ME
    R = ME.message;
end
end

%% ======================================================================
function schema = getSchemaCached(baseUrl)
% Download https://.../schema once per MATLAB session and cache on disk.
persistent SCHEMA
cacheFile = fullfile(fileparts(mfilename('fullpath')),'opentargets.schema.mat');
if isempty(SCHEMA)
    try
        SCHEMA = load(cacheFile).schema;
    catch
        % download fresh copy
        SCHEMA = string(webread(baseUrl+"/schema"));
        try save(cacheFile,'SCHEMA','-v6'); catch, end
    end
end
schema = SCHEMA;
end

%% ----------------------------------------------------------------------
function methods = listAvailableMethods(schema)
%LISTAVAILABLEMETHODS  Print only valid root-level Query methods and return them
%
%   methods = listAvailableMethods(schema) prints each root-level method
%   and returns a struct array with fields:
%       methods(i).name : string
%       methods(i).args : cell array of strings 'arg:type'
%
%   `schema` is the GraphQL SDL text as a MATLAB string or char array.

    % 1) Normalize schema to a character vector
    if isstring(schema) && isscalar(schema)
        txt = char(schema);
    elseif ischar(schema)
        txt = schema;
    else
        error('opentargets:schema','Schema must be a string or char array.');
    end

    % 2) Locate the start of the Query block
    startTok = 'type Query {';
    sPos = strfind(txt, startTok);
    if isempty(sPos)
        error('opentargets:schema', 'No ''type Query {'' in schema.');
    end
    idx = sPos(1) + numel(startTok);

    % 3) Balance braces to find the matching closing '}'
    bc = 1;             % we've just passed the first '{'
    len = numel(txt);
    endPos = 0;
    for p = idx+1:len
        switch txt(p)
            case '{'
                bc = bc + 1;
            case '}'
                bc = bc - 1;
                if bc == 0
                    endPos = p;
                    break;
                end
        end
    end
    if bc ~= 0
        error('opentargets:schema','Unbalanced braces in Query block.');
    end

    % 4) Extract only the inner text of the Query block
    blk = txt(idx+1:endPos-1);

    % 5) Use one regex to pull out all "name(arg1:Type,...)
    toks = regexp(blk, ...
        '^\s*([A-Za-z]\w*)\s*\(\s*([^)]*?)\s*\)', ...
        'tokens', 'lineanchors');

    % Initialize output
    methods = struct('name', {}, 'args', {});
    if isempty(toks)
        fprintf('No root methods found.\n');
        return;
    end

    % 6) Print and collect
    fprintf('Available root methods with arguments:\n\n');
    for k = 1:numel(toks)
        name    = toks{k}{1};
        argsStr = strtrim(toks{k}{2});
        if ~isempty(argsStr)
            argPairs = strtrim(strsplit(argsStr, ','));
            argList  = strjoin(argPairs, ', ');
        else
            argPairs = {};
            argList  = '';
        end
        fprintf('  %-20s %s\n', name, argList);
        methods(k).name = string(name);
        methods(k).args = argPairs;
    end
    fprintf('\n');
end


%% ----------------------------------------------------------------------
function m = guessMethod(query,schema)
if isstruct(query) % use first field‑name
    fn = fieldnames(query);
    m  = string(fn{1});
    return
end
% scalar query → look for single‑argument root methods
blk = extractBetween(schema,"type Query {","}");
lines = splitlines(blk); lines = strtrim(lines);
cand = string.empty;
for L = lines.'
    tok = regexp(L,'^(\w+)\([^!:]+!:','tokens','once');
    if ~isempty(tok); cand(end+1) = string(tok{1}); end %#ok<AGROW>
end
if isscalar(cand)
    m = cand;
else
    error("opentargets:method","Cannot infer unique method; candidates: %s",join(cand,","));
end
end

%% ----------------------------------------------------------------------
function data = getgraphql(query,p)
args = p.qr.args.str(p.qr.args.need);
if isstruct(query)
    fn = string(fieldnames(query));
    if any(~ismember(args,fn))
        error("opentargets:args","Query struct must contain: %s",join(args,","));
    end
    query = rmfield(query,setdiff(fn,args));
    instruct = query;
else
    if numel(args)~=1
        warning("opentargets:args","Method %s needs %d args – supply struct",p.method,numel(args));
    end
    % instruct.(args) = query;
    instruct = query;
end

data = struct('query',p.cmd,'variables',instruct);
end

%% ----------------------------------------------------------------------
function p = parseSchema(p,schema)
% ==== (original code unchanged except input arg) =========================
% file = "opentargets.schema.mat";
% no longer load/save here – schema already provided

% raw = extract(schema, p.method + "(" + wildcardPattern + ")" + wildcardPattern + lineBoundary('end')).strip;
patt = "(?s)" + p.method + "\s*\([\s\S]*?\)\s*:\s*\w+[!\]]*";
raw = regexp(schema, patt, "match", "once");

if isempty(raw)
    error('opentargets:method','Method not found in schema.');
end
raw = splitlines(raw).strip;
raw(1) = [];                 % drop header line
[raw, qr.info] = getInfo(raw);
if numel(raw)>1, raw = join(raw); end

% method's args 
[qr.args.str, qr.args.type, qr.args.need] = sepStrType(raw);

% method's type
qr.type.name = extractAfter(raw," ):"|"):").strip;
qr.type.name = regexprep(qr.type.name,'[!\[\]]','');
qr.type.fields = parseType(schema, qr.type.name);

% sub-types (and sub-fields) for the main method's type
qr.type.tree = struct;
qr.type.tree = str2struct(qr.type.fields.str', qr.type.fields.type.');
qr.type.tree = addSubfields(schema, qr.type.tree);
qr.method = p.method;
p.cmd = join(constructQuery(qr), newline);
p.qr = qr;

end

% -------------------------------------------------------------------------
function [str,type,need] = sepStrType(raw)
if numel(raw)>1, raw = join(raw); end
str  = split(extractBefore(raw,')'),',').strip;
need = endsWith(str,'!');
type = extractBetween(str,': ',textBoundary('end'));
type = erase(type,'!'+textBoundary('end'));
str  = regexprep(str,':.*$','');
end

function [raw, info] = getInfo(raw)
    % Remove blank lines
    raw(raw=="") = [];

    % Identify lines that are *only* descriptions (start & end with ")
    isDesc = startsWith(strtrim(raw), '"') & endsWith(strtrim(raw), '"');

    % Pull those out as info (strip the quotes), drop them from raw
    info    = strtrim(raw(isDesc));
    info    = regexprep(info, '^"|"$', '');   % remove leading/trailing "
    raw(isDesc) = [];

    % Done
end

% -------------------------------------------------------------------------
function tab = parseType(schema,name)
% Returns struct array with .str (field names) and .type (return types)

    % ── 1. pull the raw block ------------------------------------------------
    blk = extractBetween(schema, ...
            lineBoundary("start") + "type " + name + " {", ...
            lineBoundary("start") + "}");
    blk = char(blk);          % char vector, easier for regexp
    blk(blk == 13) = [];      % strip CRs

    if isempty(blk), tab = []; return; end

    % ── 2. regex:  field  (optionalArgs)  :  ReturnType -----------------------
    % ^\s*   – start of line, allow indent
    % ([A-Za-z_]\w*) – field name
    % \s*(?:\([\s\S]*?\))? – non-greedy “( … )” block (may span newlines)
    % \s*:\s* – colon
    % ([A-Za-z_\[\]!]+) – return type token (allow [], !)
    pat = '(?m)^\s*([A-Za-z_]\w*)\s*(?:\([\s\S]*?\))?\s*:\s*([A-Za-z_\[\]!]+)';

    toks = regexp(blk, pat, 'tokens');   % one token-pair per field

    if isempty(toks), tab = []; return; end

    str   = string(cellfun(@(c) c{1}, toks, 'uni', 0));
    ptype = string(cellfun(@(c) c{2}, toks, 'uni', 0));
    ptype = regexprep(ptype, '[!\[\]]', '');   % clean “!”, “[” “]”

    tab = struct('str', str, 'type', ptype, 'info', strings(numel(str),1));
end

%% ------------------------------------------------------------------------
function qr = addSubfields(schema, qr, visited)
% expand only OBJECT types exactly once per path

    if nargin < 3, visited = string.empty; end
    fn = string(fieldnames(qr));

    for k = 1:numel(fn)
        tname = qr.(fn(k));                     % field’s *type*

        % a) already expanded or not a plain string  → skip
        if isstruct(tname) || any(tname == visited)
            continue
        end

        % b) what kind of type is it?
        switch typeKind(schema, tname)
            case "object"                       % dive only for objects
                visited(end+1) = tname;         %#ok<AGROW>
                sub = parseType(schema, tname);
                if ~isempty(sub)
                    qr.(fn(k)) = str2struct(sub.str, sub.type);
                    qr.(fn(k)) = addSubfields(schema, qr.(fn(k)), visited);
                end
            otherwise                           % scalar / enum / union
                % treat as leaf
        end
    end
end


% function qr = addSubfields(schema,qr)
%     fn = string(fieldnames(qr));
%     for i = 1:numel(fn)
%         tt = parseType(schema, qr.(fn(i)));
%         if ~isempty(tt)
%             qr.(fn(i)) = str2struct(tt.str, tt.type);
%             qr.(fn(i)) = addSubfields(schema, qr.(fn(i)));
%         end
%     end
% end

%% ------------------------------------------------------------------------
function kind = typeKind(schema, tname)
% Returns  "object" | "enum" | "union" | "scalar"

    % ── 0. Core name & built-in scalars ────────────────────────────────
    tname = strip(tname);
    if any(tname == ["String","Int","Float","Boolean","ID"])
        kind = "scalar";   return
    end

    % ── 1. Normalise & split once ─────────────────────────────────────
    txt       = char(schema);      txt(txt==13) = [];   % kill CR
    lines     = splitlines(string(txt));                % string array
    tpattern  = "type "  + tname + " ";
    epattern  = "enum "  + tname + " ";
    upattern  = "union " + tname + " ";

    % ── 2. Starts-with check on every trimmed line ────────────────────
    if any(startsWith(strtrim(lines),  tpattern))
        kind = "object";
    elseif any(startsWith(strtrim(lines), epattern))
        kind = "enum";
    elseif any(startsWith(strtrim(lines), upattern))
        kind = "union";
    else
        kind = "scalar";          % custom scalar or undefined
    end
end


%% ------------------------------------------------------------------------
function s = str2struct(fields,vals)
    for i = 1:numel(fields)
        s.(fields(i)) = vals(i);
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
function body = constructQueryBody(qr,indent)
    if nargin<2, indent = 4; str = qr.type.tree; else, str = qr; end
    body = strings(0);
    fn = string(fieldnames(str));
    for i = 1:numel(fn)
        line = blanks(indent)+fn(i);
        if isstruct(str.(fn(i)))
            line = line+" {";
            body = [body; line]; %#ok<AGROW>
            inner = constructQueryBody(str.(fn(i)),indent+2);
            body  = [body; inner; blanks(indent) + "}"]; %#ok<AGROW>
        elseif str.(fn(i)) ~= "EntityUnionType"
            body = [body; line]; %#ok<AGROW>
        end
    end
end
