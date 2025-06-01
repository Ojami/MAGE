function backID = biomaRt(geneList, flag)
% uses biomaRt R package to do ID conversion between Ensembl and Entrez

% INPUT:
%   ensemblList: a string array or cell string containing Ensembl gene IDs.
%   flag       : either "entrez" or "ensembl"

if ~any(strcmp(flag, {'ensembl', 'entrez'}))
    fprintf('ERROR: only ensembl or entrez is allowed!\n')
    return
end
geneList = cellstr(geneList);
writecell(geneList, 'MATLAB_biomaRt.txt')

R = ({});
currDir = pwd;
currDir = strrep(currDir, '\', '/');
R{1} = ['setwd("', currDir, '")'];
R{2} = 'library(biomaRt)';
R{3} = 'en_genes <- read.delim("MATLAB_biomaRt.txt", header = F, as.is = T)';
R{4} = 'en_genes <- en_genes[,1]';
R{5} = 'mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))';

if strcmp(flag, 'ensembl')
    R{6} = ['genes <- getBM(filters="ensembl_gene_id", ',...
        'attributes=c("ensembl_gene_id", "entrezgene_id"), ',...
        'values= en_genes, mart=mart)'];
elseif strcmp(flag, 'entrez')
    R{6} = ['genes <- getBM(filters="entrezgene_id", ',...
        'attributes=c("ensembl_gene_id", "entrezgene_id"), ',...
        'values= en_genes, mart=mart)'];
end
R{7} = 'write.table(genes, file = "GENE_IDs.txt", sep = "\t", quote = F,row.names = F)';

fid = fopen('MATLAB_biomaRt.r', 'w');
fprintf(fid,'%s\n',R{:});
fclose(fid);
MATLAB2Rconnector('MATLAB_biomaRt.r')
delete('MATLAB_biomaRt.txt')

backID = readtable("GENE_IDs.txt", 'Delimiter', '\t');
delete('GENE_IDs.txt')
end