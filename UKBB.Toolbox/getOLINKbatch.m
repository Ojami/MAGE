function getOLINKbatch(plateFile, opts)
%@29JULY2024: use UKB-RAP to get plate IDs (df 30901) using table exporter.
% This function generates OLINK batch variables for adjustments in pQTL
% analyses (see https://www.nature.com/articles/s41586-023-06592-6).

arguments
    plateFile {mustBeFile} % from UKB-RAP
    opts.path {mustBeFolder} = fullfile(fileparts(which("getOLINKbatch.m")), "UKB_PHENO")
end

df = readtable(plateFile, FileType="text", TextType="string", ...
    VariableNamingRule="preserve");
df = df(:, ["eid", "30901-0.0"]);
df = renamevars(df, "30901-0.0", "plate");

% plate ID to batch mapper
pb = readtable("https://biobank.ctsu.ox.ac.uk/crystal/ukb/auxdata/olink_batch_number.dat", ...
    FileType="text", TextType="string", VariableNamingRule="preserve");

% sanity check
assert(all(ismember(df.plate, pb.PlateID)), "Some unkonw plates found in the data from UKB-RAP, check your data!")

[~, idx] = ismember(df.plate, pb.PlateID);
df.batch = pb.Batch(idx);

% write it to a struct pheno
UKB_STRUCT_ALL = struct("eid", df.eid, "rawUKB", df.batch, ...
    "numericFlag", true, "termMeaning", '', "tag", "OLINK_BATCH");

save(fullfile(opts.path, "OLINK_BATCH.mat"), "UKB_STRUCT_ALL")
fprintf("wrote to %s\n", fullfile(opts.path, "OLINK_BATCH.mat"))

end % END