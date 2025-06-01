function olink_proteomics_parser(file, index)
% @22AUG2023: parses UKBB Olink proteomics data in long format.

arguments
    file {mustBeFile} % olink_data from Data Portal
    index {mustBeFile} % Data-Coding 143
end

if isempty(gcp("nocreate"))
    parpool("Processes", 40);
end

wd = fullfile(pwd, "olink_proteome");
if ~isfolder(wd), mkdir(wd); end

% read index: uniprot <-> protein names
index = bfilereader(index, "header", true);

df = tabularTextDatastore(file, "TextType", "string", ...
    "VariableNamingRule", "preserve");

% read protein_id
matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
pid = df;
df = tall(df);
pid.SelectedVariableNames = {'protein_id'};
pid = gather(tall(pid));
pid = pid.protein_id;
upid = unique(pid);

for k = 1:numel(upid)
    tt = tic;
    fprintf("reading protein NPX %d of %d ", k, numel(upid))
    idx = pid == upid(k);
    tab = df(idx, :);
    tab = gather(tab);

    % find protein name
    pname = index.meaning(index.coding == upid(k));
    pname_short = extractBefore(pname, ";");

    ins = unique(tab.ins_index);
    for j = 1:numel(ins)
        tmp = tab(tab.ins_index == ins(j), :);
        ds = struct;
        ds.eid = tmp.eid;
        ds.rawUKB = tmp.result;
        ds.numericFlag = true;
        ds.tag = pname;
        ds.termMeaning = '';
        ds.info.basket = "olink_data";
        ds.info.date = datetime("now");
        ds.info.df = "30900";
        ds.info.dfraw = "x30900_" + ins(j) + "_0";
        UKB_STRUCT_ALL = ds;
        save(fullfile(wd, pname_short + "_" + ins(j) + ".mat"), "UKB_STRUCT_ALL")

    end

    tt = toc(tt);
    fprintf( "(done in %.0f sec)\n", tt)
end
matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);

end % END