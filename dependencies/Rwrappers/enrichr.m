function enrichr(genes, opts)
% a wrapper for enrichR R package: https://github.com/wjawaid/enrichR

arguments
    genes {mustBeText}
    opts.getdbs (1,1) logical = false % if true fetches list of available databases and modifies this function 'db' field
    opts.db {mustBeMember(opts.db, ["Genome_Browser_PWMs","TRANSFAC_and_JASPAR_PWMs","Transcription_Factor_PPIs","ChEA_2013","Drug_Perturbations_from_GEO_2014","ENCODE_TF_ChIP-seq_2014","BioCarta_2013","Reactome_2013","WikiPathways_2013","Disease_Signatures_from_GEO_up_2014","KEGG_2013","TF-LOF_Expression_from_GEO","TargetScan_microRNA","PPI_Hub_Proteins","GO_Molecular_Function_2015","GeneSigDB","Chromosome_Location","Human_Gene_Atlas","Mouse_Gene_Atlas","GO_Cellular_Component_2015","GO_Biological_Process_2015","Human_Phenotype_Ontology","Epigenomics_Roadmap_HM_ChIP-seq","KEA_2013","NURSA_Human_Endogenous_Complexome","CORUM","SILAC_Phosphoproteomics","MGI_Mammalian_Phenotype_Level_3","MGI_Mammalian_Phenotype_Level_4","Old_CMAP_up","Old_CMAP_down","OMIM_Disease","OMIM_Expanded","VirusMINT","MSigDB_Computational","MSigDB_Oncogenic_Signatures","Disease_Signatures_from_GEO_down_2014","Virus_Perturbations_from_GEO_up","Virus_Perturbations_from_GEO_down","Cancer_Cell_Line_Encyclopedia","NCI-60_Cancer_Cell_Lines","Tissue_Protein_Expression_from_ProteomicsDB","Tissue_Protein_Expression_from_Human_Proteome_Map","HMDB_Metabolites","Pfam_InterPro_Domains","GO_Biological_Process_2013","GO_Cellular_Component_2013","GO_Molecular_Function_2013","Allen_Brain_Atlas_up","ENCODE_TF_ChIP-seq_2015","ENCODE_Histone_Modifications_2015","Phosphatase_Substrates_from_DEPOD","Allen_Brain_Atlas_down","ENCODE_Histone_Modifications_2013","Achilles_fitness_increase","Achilles_fitness_decrease","MGI_Mammalian_Phenotype_2013","BioCarta_2015","HumanCyc_2015","KEGG_2015","NCI-Nature_2015","Panther_2015","WikiPathways_2015","Reactome_2015","ESCAPE","HomoloGene","Disease_Perturbations_from_GEO_down","Disease_Perturbations_from_GEO_up","Drug_Perturbations_from_GEO_down","Genes_Associated_with_NIH_Grants","Drug_Perturbations_from_GEO_up","KEA_2015","Gene_Perturbations_from_GEO_up","Gene_Perturbations_from_GEO_down","ChEA_2015","dbGaP","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Chem_Pert_down","GTEx_Tissue_Expression_Down","GTEx_Tissue_Expression_Up","Ligand_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_up","Ligand_Perturbations_from_GEO_up","MCF7_Perturbations_from_GEO_down","MCF7_Perturbations_from_GEO_up","Microbe_Perturbations_from_GEO_down","Microbe_Perturbations_from_GEO_up","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up","L1000_Kinase_and_GPCR_Perturbations_down","L1000_Kinase_and_GPCR_Perturbations_up","Reactome_2016","KEGG_2016","WikiPathways_2016","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","Kinase_Perturbations_from_GEO_down","Kinase_Perturbations_from_GEO_up","BioCarta_2016","HumanCyc_2016","NCI-Nature_2016","Panther_2016","DrugMatrix","ChEA_2016","huMAP","Jensen_TISSUES","RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO","MGI_Mammalian_Phenotype_2017","Jensen_COMPARTMENTS","Jensen_DISEASES","BioPlex_2017","GO_Cellular_Component_2017","GO_Molecular_Function_2017","GO_Biological_Process_2017","GO_Cellular_Component_2017b","GO_Molecular_Function_2017b","GO_Biological_Process_2017b","ARCHS4_Tissues","ARCHS4_Cell-lines","ARCHS4_IDG_Coexp","ARCHS4_Kinases_Coexp","ARCHS4_TFs_Coexp","SysMyo_Muscle_Gene_Sets","miRTarBase_2017","TargetScan_microRNA_2017","Enrichr_Libraries_Most_Popular_Genes","Enrichr_Submissions_TF-Gene_Coocurrence","Data_Acquisition_Method_Most_Popular_Genes","DSigDB","GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018","TF_Perturbations_Followed_by_Expression","Chromosome_Location_hg19","NIH_Funded_PIs_2017_Human_GeneRIF","NIH_Funded_PIs_2017_Human_AutoRIF","Rare_Diseases_AutoRIF_ARCHS4_Predictions","Rare_Diseases_GeneRIF_ARCHS4_Predictions","NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions","NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions","Rare_Diseases_GeneRIF_Gene_Lists","Rare_Diseases_AutoRIF_Gene_Lists","SubCell_BarCode","GWAS_Catalog_2019","WikiPathways_2019_Human","WikiPathways_2019_Mouse","TRRUST_Transcription_Factors_2019","KEGG_2019_Human","KEGG_2019_Mouse","InterPro_Domains_2019","Pfam_Domains_2019","DepMap_WG_CRISPR_Screens_Broad_CellLines_2019","DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019","MGI_Mammalian_Phenotype_Level_4_2019","UK_Biobank_GWAS_v1","BioPlanet_2019","ClinVar_2019","PheWeb_2019","DisGeNET","HMS_LINCS_KinomeScan","CCLE_Proteomics_2020","ProteomicsDB_2020","lncHUB_lncRNA_Co-Expression","Virus-Host_PPI_P-HIPSTer_2020","Elsevier_Pathway_Collection","Table_Mining_of_CRISPR_Studies","COVID-19_Related_Gene_Sets","MSigDB_Hallmark_2020","Enrichr_Users_Contributed_Lists_2020","TG_GATES_2020","Allen_Brain_Atlas_10x_scRNA_2021","Descartes_Cell_Types_and_Tissue_2021","KEGG_2021_Human","WikiPathway_2021_Human","HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression","GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021","MGI_Mammalian_Phenotype_Level_4_2021","CellMarker_Augmented_2021","Orphanet_Augmented_2021","COVID-19_Related_Gene_Sets_2021","PanglaoDB_Augmented_2021","Azimuth_Cell_Types_2021","PhenGenI_Association_2021","RNAseq_Automatic_GEO_Signatures_Human_Down","RNAseq_Automatic_GEO_Signatures_Human_Up","RNAseq_Automatic_GEO_Signatures_Mouse_Down","RNAseq_Automatic_GEO_Signatures_Mouse_Up","GTEx_Aging_Signatures_2021","HDSigDB_Human_2021","HDSigDB_Mouse_2021","HuBMAP_ASCTplusB_augmented_2022","FANTOM6_lncRNA_KD_DEGs","MAGMA_Drugs_and_Diseases","PFOCR_Pathways","Tabula_Sapiens","ChEA_2022","Diabetes_Perturbations_GEO_2022","LINCS_L1000_Chem_Pert_Consensus_Sigs","LINCS_L1000_CRISPR_KO_Consensus_Sigs","Tabula_Muris","Reactome_2022","SynGO_2022","GlyGen_Glycosylated_Proteins_2022","IDG_Drug_Targets_2022","KOMP2_Mouse_Phenotypes_2022","Metabolomics_Workbench_Metabolites_2022","Proteomics_Drug_Atlas_2023","The_Kinase_Library_2023","GTEx_Tissues_V8_2023","GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023","PFOCR_Pathways_2023","GWAS_Catalog_2023","GeDiPNet_2023","MAGNET_2023"])} = "Genome_Browser_PWMs"
    opts.output {mustBeTextScalar} = "enrichr_results"
    opts.maxtry = 20 % maximum number of times to try to reconnect if connection fails

    % for plot
    opts.plot = true
    opts.res = 300
    opts.width = 7
    opts.height = 8
    opts.y {mustBeMember(opts.y, ["Count", "Ratio"])} = "Ratio"
    opts.padj = 0.05 % filter rows based on this value
end

if opts.getdbs || ~isfield(opts, "db")
    getEnrichrDB;
    return
end

[hm, name, ext] = fileparts(opts.output);
if isempty(hm) || hm == "", hm = pwd; end
if ~isfolder(hm), mkdir(hm); end
opts.output = replace(fullfile(hm, name + ext), filesep, "/");
opts.output = regexprep(opts.output, [".xls$", ".xlsx$"], "");

fname = getRandomName("enrich", 4);
writematrix(string(genes), fname + ".txt", "QuoteStrings", "none")

rcode(1) = "httr::timeout(1e6)";
rcode(2) = "require('enrichR')";
rcode(3) = 'dbs <- c("' + join(opts.db, '","') + '")';
rcode(4) = "genes = unlist(strsplit(readLines('" + fname + ".txt'), ','))";
rcode(5) = "enum <- FALSE";
rcode(6) = "for (i in 1:" + opts.maxtry + "){";
rcode(7) = "tryCatch( { enriched <- enrichr(genes, dbs) }";
rcode(8) = "       , error = function(e) {enum <<- TRUE})";
rcode(9) = "if(!enum){break}}";
rcode(10) = "writexl::write_xlsx(enriched, '" + opts.output + ".xlsx')";

if opts.plot
    rcode(11) = "ns = names(enriched)";
    rcode(12) = "for (i in seq_along(ns)) {";
    rcode(13) = "pl = plotEnrich(enriched[[i]], showTerms = 20, " + ...
        "numChar = 40, y = '" + opts.y + "', orderBy = 'P.value', title = '')";
    rcode(14) = "png(paste0('" + opts.output + "', ns[i], '.png'), width = " + opts.width + ...
        ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
    rcode(15) = "print(pl)";
    rcode(16) = "dev.off()}";
end

MATLAB2Rconnector(fname + ".r", code=rcode', delr=true);
delete(fname + ".txt")

if opts.padj < 1 % filter output
    opts.output = opts.output.replace("/", filesep) + ".xlsx";
    sheets = sheetnames(opts.output);
    for k = 1:numel(sheets)
        tmp = readtable(opts.output, Sheet=sheets(k), TextType="string", ...
            VariableNamingRule="preserve");
        tmp(tmp.("Adjusted.P.value") > opts.padj, :) = [];
        writetable(tmp, opts.output, Sheet=sheets(k), WriteMode="overwritesheet")
    end
end

end % END

%% subfunctions ===========================================================
function getEnrichrDB

r(1, 1) = "httr::timeout(1e6)";
r(2, 1) = "require('enrichR')";
r(3, 1) = "enum <- FALSE";
r(4, 1) = "for (i in 1:10){";
r(5, 1) = "tryCatch( { dbs = listEnrichrDbs()$libraryName}";
r(6, 1) = "       , error = function(e) {enum <<- TRUE})";
r(7, 1) = "if(!enum){break}}";
r(8, 1) = "data.table::fwrite(as.list(dbs), 'enrichr_dbs.txt')";

MATLAB2Rconnector("enrichr_db.r", code=r, delr=true);
dbs = readmatrix("enrichr_dbs.txt", "OutputType", "string");
delete("enrichr_dbs.txt")

fn = readlines("enrichr.m");
idx = find(startsWith(fn, whitespacePattern + "opts.db"), 1);
hdr = extractBefore(fn(idx), "opts.db") + "opts.db";
hdr = hdr + ' {mustBeMember(opts.db, ["' + join(dbs, '","') + ...
    '"])} = "' + dbs(1) + '"';
fn(idx) = hdr;
writelines(fn, which("enrichr.m"))

end
