function KEGGview(gene_list, pathids, opts)
% maps data to KEGG pathways using R library "pathview".
% Also depends on "clusterProfiler" package. 
% 
% Oveis Jamialahmad, GU, Oct 2020.
% 
% @04/10/2022: pathids can be left empty, and the function used KEGG API to
%              map input genes to KEGG maps.
% @03MAY2024: 'minN' was added to only draw maps with a minimum number of
%              overlapping genes between gene_list and target pathway 
%              genes.

arguments
    gene_list {mustBeA(gene_list, 'table')} % gene_list must be a table, with id and value columns
    pathids {mustBeText, mustBeVector} = "" % KEGG pathway ids for input gene list
    opts.cpd_data {mustBeA(opts.cpd_data, 'table')} % similar to gene_list but KEGG cpd IDs should be in the first column.
    opts.genome {mustBeMember(opts.genome, ["37", "38"])} = "38" % for 'findPathIds' subfunction
    opts.minN_gene(1,1) double = 3 % minimum number of overlapping genes in the target pathway for visualization
    opts.ratio_gene (1,1) double % overrides 'minN_gene': similar to minN but minimum ratio of overlapping genes
    opts.minN_cpd(1,1) double = 3 % minimum number of overlapping compounds in the target pathway for visualization
    opts.ratio_cpd (1,1) double % overrides 'minN_cpd': similar to minN but minimum ratio of overlapping compounds
    opts.ignore {mustBeText, mustBeVector} % pathways to ignore
    opts.outdir {mustBeTextScalar} % output dir
    opts.forceInnerAPI (1,1) logical = false % to get overlapping genes with each pathway. This is useful, when user wants to have different colorbars per each pathway.
    opts.join {mustBeMember(opts.join, ["union", "intersect"])} = "intersect" % if pathids are provided and 'forceInnerAPI' is true, how to keep 'pathids' and overlapping pathways?
end

if ~isfield(opts, 'outdir')
    opts.outdir = pwd;
end

if ~isfolder(opts.outdir), mkdir(opts.outdir); end
 
gene_list.Properties.VariableNames = {'id', 'value'};
gene_list.id = string(gene_list.id);
if all(startsWith(gene_list.id, "ENSG"))
    fromType = "ENSEMBL";
elseif all(~isnan(double(string(gene_list.id))))
    fromType = "ENTREZID";
else
    fromType = "SYMBOL";
end

% check cpd_data
if isfield(opts, 'cpd_data')
    cpd_list = rmmissing(opts.cpd_data);
    cpd_list.Properties.VariableNames = {'id', 'value'};
    cpd_list.id = string(cpd_list.id);
    cpdValRange = [min(cpd_list.value), max(cpd_list.value)];
    cpdValRange = [floor(cpdValRange(1)*10)/10, ceil(cpdValRange(2)*10)/10];
    opts.cpd_data = cpd_list;
    writetable(cpd_list, fullfile(opts.outdir, 'KEGGview_CPDLIST.txt'), 'Delimiter','\t')
else
    cpdValRange = [-1, 1];
end

writetable(gene_list, fullfile(opts.outdir, 'KEGGview_TMPLIST.txt'), 'Delimiter','\t')
genevalRange = [min(gene_list.value), max(gene_list.value)];
genevalRange = [floor(genevalRange(1)*10)/10, ceil(genevalRange(2)*10)/10];

% create a new table
if fromType ~= "ENTREZID"
    Rdata = strings;
    Rdata(1, 1) = "setwd('" + replace(opts.outdir, '\', '/') + "')";
    Rdata(2, 1) = "df<-read.table(file ='KEGGview_TMPLIST.txt', header = T, sep = '\t')";
    Rdata(3, 1) = 'packages <- c("pathview", "clusterProfiler")';
    Rdata(4, 1) = "if(length(setdiff(packages, rownames(installed.packages()))) != 0){BiocManager::install(setdiff(packages, rownames(installed.packages())))}";
    Rdata(5, 1) = "require(pathview)";
    Rdata(6, 1) = "require(clusterProfiler)";
    Rdata(7, 1) = "library(""org.Hs.eg.db"", character.only = TRUE)";
    Rdata(8, 1) = "ids<-bitr(df$id, " + "fromType = '" + fromType + ...
        "', toType = ""ENTREZID"", OrgDb=org.Hs.eg.db)";
    Rdata(9, 1) = "dedup_ids = ids[!duplicated(ids$" + fromType + "),]";
    Rdata(10, 1) = "df = df[df$id %in% dedup_ids$" + fromType + ",]";
    Rdata(11, 1) = "df$id = dedup_ids$ENTREZID";
    Rdata(12, 1) = "write.table(df, 'KEGGview_TMPLIST.txt', sep = '\t', row.names = F, col.names = T, quote = F)";
    MATLAB2Rconnector('KEGGview.r', "code", Rdata);
    gene_list = readtable(fullfile(opts.outdir, "KEGGview_TMPLIST.txt"), "TextType", "string", "Delimiter", "\t");
    gene_list.id = string(gene_list.id);
    delete(fullfile(opts.outdir, 'KEGGview_TMPLIST.txt'));
    writetable(gene_list, fullfile(opts.outdir, 'KEGGview_TMPLIST.txt'), 'Delimiter','\t')
    fromType = "ENTREZID";

end

opts.pathids = pathids; % to fill in, if provided and 'forceInnerAPI' is true

if any(pathids == "") || opts.forceInnerAPI % find path ids internally using KEGG API
    [pathids, opts] = findPathIds(gene_list.id, fromType, opts);
end

if all(opts.pathids ~= "")
    if opts.forceInnerAPI
        pathids = feval(opts.join, opts.pathids, pathids);
    else
        pathids = opts.pathids;
    end
end

if isfield(opts, "ignore")
    pathids = setdiff(pathids, opts.ignore);
end

Rdata = strings;
Rdata(1, 1) = "setwd('" + replace(opts.outdir, '\', '/') + "')";
Rdata(2, 1) = "df<-read.table(file ='KEGGview_TMPLIST.txt', header = T, sep = '\t')";
Rdata(3, 1) = 'packages <- c("pathview", "clusterProfiler")';
Rdata(4, 1) = "if(length(setdiff(packages, rownames(installed.packages()))) != 0){BiocManager::install(setdiff(packages, rownames(installed.packages())))}";
Rdata(5, 1) = "require(pathview)";
Rdata(6, 1) = "require(clusterProfiler)";
Rdata(7, 1) = "library(""org.Hs.eg.db"", character.only = TRUE)";
% if fromType ~= "ENTREZID"
%     Rdata(8, 1) = "ids<-bitr(df$id, " + "fromType = '" + fromType + ...
%         "', toType = ""ENTREZID"", OrgDb=org.Hs.eg.db)";
%     Rdata(9, 1) = "dedup_ids = ids[!duplicated(ids$" + fromType + "),]";
%     Rdata(10, 1) = "df = df[df$id %in% dedup_ids$" + fromType + ",]";
%     Rdata(11, 1) = "df$ENTREZ = dedup_ids$ENTREZID";
% else
Rdata(8, 1) = 'colnames(df)[1] = "ENTREZID"';
% Rdata(9, 1) = "df = df[!duplicated(df$" + fromType + "),]";
Rdata(9, 1) = "df = df[!duplicated(df$ENTREZID),]";

% end

Rdata(12, 1) = "kegg_gene_list <- df$value";
Rdata(13, 1) = "names(kegg_gene_list) <- df$ENTREZ";
Rdata(14, 1) = "kegg_gene_list<-na.omit(kegg_gene_list)";

if isfield(opts, 'cpd_data')
    Rdata(15, 1) = "df2 <-read.table(file ='KEGGview_CPDLIST.txt', header = T, sep = '\t')";
    Rdata(16, 1) = "df2 = df2[!duplicated(df2$id),]";
    Rdata(17, 1) = "kegg_cpd_list <- df2$value";
    Rdata(18, 1) = "names(kegg_cpd_list) <- df2$id";
    Rdata(19, 1) = "kegg_cpd_list<-na.omit(kegg_cpd_list)";
else
    Rdata(15, 1) = "kegg_cpd_list = NULL";
end

pathids = string(pathids);
ct = numel(Rdata) + 1;
for i = 1:numel(pathids)
    Rdata(ct, 1) = "tryCatch({";
        
    if isfield(opts, 'gp')
        gtab = opts.gp(opts.gp.path == pathids(i), :);
        if isempty(gtab)
            Rdata(ct + 1, 1) = "tmpG = NULL"; 
        else
            idx = ismember("hsa:" + gene_list.id, gtab.gene{1});
            gtab = gene_list.id(idx);
            
            if isempty(gtab)
                Rdata(ct + 1, 1) = "tmpG = NULL"; 
            else
                genevalRange = [min(gene_list.value(idx)), max(gene_list.value(idx))];
                genevalRange = [floor(genevalRange(1)*10)/10, ceil(genevalRange(2)*10)/10];
                gtab = createRvec(gtab, 10, false);
                gtab(1) = "tmpG = " + gtab(1);
                for m = 1:numel(gtab)
                    Rdata(ct + 1, 1) = gtab(m);
                    ct = ct + 1;
                end
                Rdata(ct + 1, 1) = "tmpG = subset(kegg_gene_list, names(kegg_gene_list) %in% tmpG)"; 
                
            end
        end
        ct = ct + 1;

        if isfield(opts, 'gc')
            ctab = opts.gc(opts.gc.path == pathids(i), :);
            if isempty(ctab)
                Rdata(ct + 1, 1) = "tmpC = NULL"; 
            else
                idx = ismember(cpd_list.id, ctab.cpd{1});
                ctab = cpd_list.id(idx);
                
                if isempty(ctab)
                    Rdata(ct + 1, 1) = "tmpC = NULL"; 
                else
                    cpdValRange = [min(cpd_list.value(idx)), max(cpd_list.value(idx))];
                    cpdValRange = [floor(cpdValRange(1)*10)/10, ceil(cpdValRange(2)*10)/10];
                    ctab = createRvec(ctab, 10, true);
                    ctab(1) = "tmpC = " + ctab(1);
                    for m = 1:numel(ctab)
                        Rdata(ct + 1, 1) = ctab(m);
                        ct = ct + 1;
                    end
                    Rdata(ct + 1, 1) = "tmpC = subset(kegg_cpd_list, names(kegg_cpd_list) %in% tmpC)"; 
                end
            end
        else
            Rdata(ct + 1, 1) = "tmpC = NULL"; 
            
        end
        ct = ct + 1;

        Rdata(ct + 1, 1) = "pth <- pathview(gene.data=tmpG," + ...
            "cpd.data=tmpC, pathway.id='" + pathids(i) + "', species = ""hsa""" + ...
            ",limit = list(gene = " + createRvec(genevalRange) + ", cpd = " + createRvec(cpdValRange) + "), " + ...
            "low=list(gene = 'blue', cpd='yellow'),high=list(gene = 'red', cpd='green'))";
    else
        Rdata(ct + 1, 1) = "pth <- pathview(gene.data=kegg_gene_list," + ...
            "cpd.data=kegg_cpd_list, pathway.id='" + pathids(i) + "', species = ""hsa""" + ...
            ",limit = list(gene = " + createRvec(genevalRange) + ", cpd = " + createRvec(cpdValRange) + "), " + ...
            "low=list(gene = 'blue', cpd='yellow'),high=list(gene = 'red', cpd='green'))";
    end

    Rdata(ct + 2) = "}, error = function(e) {";
    Rdata(ct + 3) = 'cat("An error occurred:", conditionMessage(e), "\n")})';
    Rdata(ct + 4) = "";
    ct = numel(Rdata) + 1;
end

MATLAB2Rconnector('KEGGview.r', "code", Rdata);
delete(fullfile(opts.outdir, 'KEGGview_TMPLIST.txt'))

for i = 1:numel(pathids)
    delete(fullfile(opts.outdir, pathids(i) + ".png"))
    delete(fullfile(opts.outdir, pathids(i) + ".xml"))
end
end % END

%% subfunctions ===========================================================
function [pathids, opts] = findPathIds(id, tp, opts)

if strcmpi(tp, "ensembl")
    tid = EnsemblREST(id, "lookup", "refGenome", opts.genome);
    [~, idx] = ismember(id, tid.id);
    tid = tid.name(idx(idx > 0));
elseif strcmpi(tp, "symbol")
    tid = id;
end

% find KEGG gene ids
hsa = webread("https://rest.kegg.jp/list/hsa", weboptions("Timeout", 1e5));
hsa = splitlines(hsa);
hsa = string(hsa);
hsa(hsa == "") = [];
hsa = split(hsa, char(9));
hsa = table(hsa(:, 1), hsa(:, 2), hsa(:, 3), ...
    arrayfun(@(x)split(x,[",",";"]).strip, hsa(:, 4),"uni",false), ...
    'VariableNames', {'id', 'type', 'pos', 'name'});

if ~strcmpi(tp, "entrezid")
    idx = cellfun(@(x)ismember(tid, x), hsa.name, 'uni', false);
    empidx = ~cellfun(@any, idx);
    hsa(empidx, :) = []; kegg_gene = hsa.id;
else
    idx = ismember(hsa.id, "hsa:" + string(id));
    kegg_gene = hsa.id(idx);
end

% map genes to pathways
g2p = webread("https://rest.kegg.jp/link/pathway/hsa", weboptions("Timeout", 1e5));
g2p = string(splitlines(g2p)); g2p(g2p == "") = [];
g2p = split(g2p);

% minimum number of genes
gp = array2table(g2p, "VariableNames", ["gene", "path"]);
gp = groupsummary(gp, "path", @(x)join(x, ";"), "gene");
gp.gene = arrayfun(@(x)split(x, ";"), gp.(3), uni=false);
if isfield(opts, "ratio_gene")
    c = cellfun(@(x)sum(ismember(x, kegg_gene))/numel(x), gp.gene, uni=false);
    cutoff = opts.ratio_gene;
else
    c = cellfun(@(x)sum(ismember(x, kegg_gene)), gp.gene, uni=false);
    cutoff = opts.minN_gene;
end
c = vertcat(c{:});
idx = c > cutoff;
if isfield(opts, 'pathids')
    idx = idx | ismember(gp.path, "path:" + opts.pathids);
end
pathids = gp.path(idx);
opts.gp = gp;

% map compounds to pathways
if isfield(opts, 'cpd_data')
    c2p_1 = webread("https://rest.kegg.jp/link/cpd/path", weboptions("Timeout", 1e5));
    c2p_1 = string(splitlines(c2p_1)); c2p_1(c2p_1 == "") = [];
    c2p_1 = split(c2p_1);
    c2p_2 = webread("https://rest.kegg.jp/link/gl/path", weboptions("Timeout", 1e5));
    c2p_2 = string(splitlines(c2p_2)); c2p_2(c2p_2 == "") = [];
    c2p_2 = split(c2p_2);
    c2p_3 = webread("https://rest.kegg.jp/link/dr/path", weboptions("Timeout", 1e5));
    c2p_3 = string(splitlines(c2p_3)); c2p_3(c2p_3 == "") = [];
    c2p_3 = split(c2p_3);
    c2p = [c2p_1; c2p_2; c2p_3];
    gc = array2table(c2p, "VariableNames", ["path", "cpd"]);
    gc.cpd = extractAfter(gc.cpd, ":");
    gc = groupsummary(gc, "path", @(x)join(x, ";"), "cpd");
    gc.cpd = arrayfun(@(x)split(x, ";"), gc.(3), uni=false);
    kegg_cpd = opts.cpd_data.id;
    if isfield(opts, "ratio_cpd")
        c = cellfun(@(x)sum(ismember(x, kegg_cpd))/numel(x), gc.cpd, uni=false);
        cutoff = opts.ratio_cpd;
    else
        c = cellfun(@(x)sum(ismember(x, kegg_cpd)), gc.cpd, uni=false);
        cutoff = opts.minN_cpd;
    end
    c = vertcat(c{:});
    gc.path = replace(gc.path, ":map", ":hsa");
    pathids = union(pathids, gc.path(c > cutoff));
    opts.gc = gc;
    opts.gc(~ismember(opts.gc.path, pathids), :) = [];
    opts.gc.path = regexprep(opts.gc.path, "^path:", "");
end
opts.gp(~ismember(opts.gp.path, pathids), :) = [];
opts.gp.path = regexprep(opts.gp.path, "^path:", "");

% pathids = unique(g2p(ismember(g2p(:, 1), kegg_gene), 2));
pathids = regexprep(pathids, "^path:", "");

end