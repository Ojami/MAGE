function getOLINKbatch(opts)
%@29JULY2024: 
% This function generates OLINK batch variables for adjustments in pQTL
% analyses (see https://www.nature.com/articles/s41586-023-06592-6).

arguments
    opts.path {mustBeFolder} = fullfile(fileparts(which("getOLINKbatch.m")), "UKB_PHENO")
end

% read Olink plate IDs
tmp = phenoParser(query="", ...
    df="30901", ...
    instance=0, ...
    save=false, surv=false);
df = struct(eid=tmp.eid{1}, plate=tmp.rawUKB{1});
clear tmp
df = struct2table(df);

% plate ID to batch mapper
url = "https://biobank.ctsu.ox.ac.uk/crystal/ukb/auxdata/olink_batch_number.dat";
topts = detectImportOptions(url, ...
            FileType="text", ...
            VariableNamingRule="preserve");
topts = setvartype(topts, topts.VariableNames, "string");
pb = readtable(url, topts); 

if all(~pb.PlateID.startsWith("0"))
    pb.PlateID = "0" + pb.PlateID;
end

% sanity check
assert(all(ismember(df.plate, pb.PlateID)), "Some unkonw plates found in the data from UKB-RAP, check your data!")

[~, idx] = ismember(df.plate, pb.PlateID);
df.batch = pb.Batch(idx);

% write it to a struct pheno
df.batch = double(df.batch);
UKB_STRUCT_ALL = struct("eid", df.eid, "rawUKB", df.batch, ...
    "numericFlag", true, "termMeaning", '', "tag", "OLINK_BATCH");

save(fullfile(opts.path, "OLINK_BATCH.mat"), "UKB_STRUCT_ALL")
fprintf("wrote to %s\n", fullfile(opts.path, "OLINK_BATCH.mat"))

end % END