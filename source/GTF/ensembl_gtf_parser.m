function ensembl_gtf_parser(file)
% Parses gene set file (GTF) from Ensembl database.
% Can be downloaded from http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
% for GRCh37: https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/

arguments
    file {mustBeFile}
end

% file = "Homo_sapiens.GRCh38.107.gtf.gz";
gunzip(file);
file = regexprep(file, ".gz$", "");

% header info for GTF format: http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/README
hdr = ["seqname", "source", "feature", "start", "stop", "score", "strand", ...
    "frame", "attribute"];
colIdx = [1, 3, 4, 5, 7, 9];

tab = tabularTextDatastore(file, "TextType" ,"string", "CommentStyle", "#", "FileExtensions", ".gtf");
tab.SelectedVariableNames = tab.SelectedVariableNames(colIdx);
tab.SelectedFormats{1} = '%q';
tab = tall(tab);
tab.Properties.VariableNames = hdr(colIdx);
tab.gene_id = extractBetween(tab.attribute, 'gene_id "', '"');
tab.gene_biotype = extractBetween(tab.attribute, 'gene_biotype "', '"');

tab.gene_name = tab.gene_biotype;
tab.gene_name(:) = missing;
idx = contains(tab.attribute, 'gene_name "');
tab.gene_name(idx) = extractBetween(tab.attribute(idx), 'gene_name "', '"');

tab.transcript_id = tab.gene_biotype;
tab.transcript_id(:) = missing;
idx = contains(tab.attribute, 'transcript_id "');
tab.transcript_id(idx) = extractBetween(tab.attribute(idx), 'transcript_id "', '"');

tab.exon_number = tab.gene_biotype;
tab.exon_number(:) = missing;
idx = contains(tab.attribute, 'exon_number "');
tab.exon_number(idx) = extractBetween(tab.attribute(idx), 'exon_number "', '"');
tab.exon_number = double(tab.exon_number);

tab.attribute = [];

features = gather(unique(tab.feature));
for i = 1:numel(features)
    fprintf("(%d of %d)-feature: %s\n", i, numel(features), features(i))
    res = gather(tab(tab.feature == features(i), :));
    res = unique(res, "rows", "stable");
    res.feature = [];
    res = table2struct(res, "ToScalar", true);
    res = rmmissing(res, 2, "MinNumMissing", height(res));
    save(file + "." + features(i) + ".mat", "-struct", "res")
end
delete(file)

end % END