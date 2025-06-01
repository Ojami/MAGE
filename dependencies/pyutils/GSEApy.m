function [res, rank_res] = GSEApy(opts)
% a wrapper for GSEApy python package:
% https://gseapy.readthedocs.io/en/latest/gseapy_example.html
% Oveis Jamialahmadi, GU, July 2023.

arguments
    opts.getdbs (1,1) logical = false % if true fetches list of available databases and modifies this function 'lib' field

    opts.rankedGenes {mustBeA(opts.rankedGenes, "table")} % must be table with two columns (order doesn't matter), one is gene name (e.g. 'genes' and the other one is rank metric e.g. log2FC)
    opts.genes {mustBeText, mustBeVector} % list of input genes, e.g. ["PNPLA3", "MTTP", "APOB"]
    opts.background {mustBeText, mustBeVector} % list of background genes
    opts.lib {mustBeMember(opts.lib, ["ARCHS4_Cell-lines","ARCHS4_IDG_Coexp","ARCHS4_Kinases_Coexp","ARCHS4_TFs_Coexp","ARCHS4_Tissues","Achilles_fitness_decrease","Achilles_fitness_increase","Aging_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_up","Allen_Brain_Atlas_10x_scRNA_2021","Allen_Brain_Atlas_down","Allen_Brain_Atlas_up","Azimuth_2023","Azimuth_Cell_Types_2021","BioCarta_2013","BioCarta_2015","BioCarta_2016","BioPlanet_2019","BioPlex_2017","CCLE_Proteomics_2020","CORUM","COVID-19_Related_Gene_Sets","COVID-19_Related_Gene_Sets_2021","Cancer_Cell_Line_Encyclopedia","CellMarker_2024","CellMarker_Augmented_2021","ChEA_2013","ChEA_2015","ChEA_2016","ChEA_2022","Chromosome_Location","Chromosome_Location_hg19","ClinVar_2019","DGIdb_Drug_Targets_2024","DSigDB","Data_Acquisition_Method_Most_Popular_Genes","DepMap_CRISPR_GeneDependency_CellLines_2023","DepMap_WG_CRISPR_Screens_Broad_CellLines_2019","DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019","Descartes_Cell_Types_and_Tissue_2021","Diabetes_Perturbations_GEO_2022","DisGeNET","Disease_Perturbations_from_GEO_down","Disease_Perturbations_from_GEO_up","Disease_Signatures_from_GEO_down_2014","Disease_Signatures_from_GEO_up_2014","DrugMatrix","Drug_Perturbations_from_GEO_2014","Drug_Perturbations_from_GEO_down","Drug_Perturbations_from_GEO_up","ENCODE_Histone_Modifications_2013","ENCODE_Histone_Modifications_2015","ENCODE_TF_ChIP-seq_2014","ENCODE_TF_ChIP-seq_2015","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","ESCAPE","Elsevier_Pathway_Collection","Enrichr_Libraries_Most_Popular_Genes","Enrichr_Submissions_TF-Gene_Coocurrence","Enrichr_Users_Contributed_Lists_2020","Epigenomics_Roadmap_HM_ChIP-seq","FANTOM6_lncRNA_KD_DEGs","GO_Biological_Process_2013","GO_Biological_Process_2015","GO_Biological_Process_2017","GO_Biological_Process_2017b","GO_Biological_Process_2018","GO_Biological_Process_2021","GO_Biological_Process_2023","GO_Biological_Process_2025","GO_Cellular_Component_2013","GO_Cellular_Component_2015","GO_Cellular_Component_2017","GO_Cellular_Component_2017b","GO_Cellular_Component_2018","GO_Cellular_Component_2021","GO_Cellular_Component_2023","GO_Cellular_Component_2025","GO_Molecular_Function_2013","GO_Molecular_Function_2015","GO_Molecular_Function_2017","GO_Molecular_Function_2017b","GO_Molecular_Function_2018","GO_Molecular_Function_2021","GO_Molecular_Function_2023","GO_Molecular_Function_2025","GTEx_Aging_Signatures_2021","GTEx_Tissue_Expression_Down","GTEx_Tissue_Expression_Up","GTEx_Tissues_V8_2023","GWAS_Catalog_2019","GWAS_Catalog_2023","GeDiPNet_2023","GeneSigDB","Gene_Perturbations_from_GEO_down","Gene_Perturbations_from_GEO_up","Genes_Associated_with_NIH_Grants","Genome_Browser_PWMs","GlyGen_Glycosylated_Proteins_2022","HDSigDB_Human_2021","HDSigDB_Mouse_2021","HMDB_Metabolites","HMS_LINCS_KinomeScan","HomoloGene","HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression","HuBMAP_ASCTplusB_augmented_2022","HumanCyc_2015","HumanCyc_2016","Human_Gene_Atlas","Human_Phenotype_Ontology","IDG_Drug_Targets_2022","InterPro_Domains_2019","Jensen_COMPARTMENTS","Jensen_DISEASES","Jensen_DISEASES_Curated_2025","Jensen_DISEASES_Experimental_2025","Jensen_TISSUES","KEA_2013","KEA_2015","KEGG_2013","KEGG_2015","KEGG_2016","KEGG_2019_Human","KEGG_2019_Mouse","KEGG_2021_Human","KOMP2_Mouse_Phenotypes_2022","Kinase_Perturbations_from_GEO_down","Kinase_Perturbations_from_GEO_up","L1000_Kinase_and_GPCR_Perturbations_down","L1000_Kinase_and_GPCR_Perturbations_up","LINCS_L1000_CRISPR_KO_Consensus_Sigs","LINCS_L1000_Chem_Pert_Consensus_Sigs","LINCS_L1000_Chem_Pert_down","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up","Ligand_Perturbations_from_GEO_down","Ligand_Perturbations_from_GEO_up","MAGMA_Drugs_and_Diseases","MAGNET_2023","MCF7_Perturbations_from_GEO_down","MCF7_Perturbations_from_GEO_up","MGI_Mammalian_Phenotype_2013","MGI_Mammalian_Phenotype_2017","MGI_Mammalian_Phenotype_Level_3","MGI_Mammalian_Phenotype_Level_4","MGI_Mammalian_Phenotype_Level_4_2019","MGI_Mammalian_Phenotype_Level_4_2021","MGI_Mammalian_Phenotype_Level_4_2024","MSigDB_Computational","MSigDB_Hallmark_2020","MSigDB_Oncogenic_Signatures","Metabolomics_Workbench_Metabolites_2022","Microbe_Perturbations_from_GEO_down","Microbe_Perturbations_from_GEO_up","MoTrPAC_2023","Mouse_Gene_Atlas","NCI-60_Cancer_Cell_Lines","NCI-Nature_2016","NIBR_DRUGseq_2025_down","NIBR_DRUGseq_2025_up","NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions","NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions","NIH_Funded_PIs_2017_Human_AutoRIF","NIH_Funded_PIs_2017_Human_GeneRIF","NURSA_Human_Endogenous_Complexome","OMIM_Disease","OMIM_Expanded","Old_CMAP_down","Old_CMAP_up","Orphanet_Augmented_2021","PFOCR_Pathways","PFOCR_Pathways_2023","PPI_Hub_Proteins","PanglaoDB_Augmented_2021","Panther_2015","Panther_2016","PerturbAtlas","Pfam_Domains_2019","Pfam_InterPro_Domains","PheWeb_2019","PhenGenI_Association_2021","Phosphatase_Substrates_from_DEPOD","ProteomicsDB_2020","Proteomics_Drug_Atlas_2023","RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO","RNAseq_Automatic_GEO_Signatures_Human_Down","RNAseq_Automatic_GEO_Signatures_Human_Up","RNAseq_Automatic_GEO_Signatures_Mouse_Down","RNAseq_Automatic_GEO_Signatures_Mouse_Up","Rare_Diseases_AutoRIF_ARCHS4_Predictions","Rare_Diseases_AutoRIF_Gene_Lists","Rare_Diseases_GeneRIF_ARCHS4_Predictions","Rare_Diseases_GeneRIF_Gene_Lists","Reactome_2013","Reactome_2015","Reactome_2016","Reactome_2022","Reactome_Pathways_2024","Rummagene_kinases","Rummagene_signatures","Rummagene_transcription_factors","SILAC_Phosphoproteomics","SubCell_BarCode","SynGO_2022","SynGO_2024","SysMyo_Muscle_Gene_Sets","TF-LOF_Expression_from_GEO","TF_Perturbations_Followed_by_Expression","TG_GATES_2020","TRANSFAC_and_JASPAR_PWMs","TRRUST_Transcription_Factors_2019","Table_Mining_of_CRISPR_Studies","Tabula_Muris","Tabula_Sapiens","TargetScan_microRNA","TargetScan_microRNA_2017","The_Kinase_Library_2023","The_Kinase_Library_2024","Tissue_Protein_Expression_from_Human_Proteome_Map","Tissue_Protein_Expression_from_ProteomicsDB","Transcription_Factor_PPIs","UK_Biobank_GWAS_v1","Virus-Host_PPI_P-HIPSTer_2020","VirusMINT","Virus_Perturbations_from_GEO_down","Virus_Perturbations_from_GEO_up","WikiPathway_2021_Human","WikiPathway_2023_Human","WikiPathways_2013","WikiPathways_2015","WikiPathways_2016","WikiPathways_2019_Human","WikiPathways_2019_Mouse","WikiPathways_2024_Human","WikiPathways_2024_Mouse","dbGaP","huMAP","lncHUB_lncRNA_Co-Expression","miRTarBase_2017"])} = "ARCHS4_Cell-lines"
    opts.organism {mustBeTextScalar} = "human"
    opts.verbose (1,1) logical = true

    % GSEA params
    opts.threads (1,1) double = 40
    opts.min_size (1,1) double = 10
    opts.max_size (1,1) double = 500
    opts.permutation (1,1) double = 5000
    opts.seed (1,1) double = 123

    % output
    opts.plot (1,1) logical = false % only for GSEA
    opts.write (1,1) logical = true % to write output tables to xlsx files
    opts.output {mustBeTextScalar} = fullfile(pwd, "GSEApy_results", "gsea") % path with file name to output where files should be saved
    opts.qvalue (1,1) double = 0.05 % q-value cutoff for output
end

% terminate(pyenv)
pd = py.importlib.import_module("pandas");
gp = py.importlib.import_module("gseapy");

if opts.getdbs || ~isfield(opts, "lib")
    getGSEApyDB(gp);
    return
end

% check necessary argument: genes or rankedGenes
if ~any(isfield(opts, ["genes", "rankedGenes"]))
    error("eith 'genes' or 'rankedGenes' should be set!")
elseif isfield(opts, "genes")
    idx = ismissing(opts.genes) | opts.genes == "";
    opts.genes(idx) = [];
    if isempty(opts.genes)
        fprintf("no gene remaiend after removing missing names!\n")
    elseif iscolumn(opts.genes)
        opts.genes = opts.genes';
    end
end

opts.output = string(opts.output);
opts.output = regexprep(opts.output, [".xlsx$", ".xls$", ".txt$"], "");
pth = fileparts(opts.output);
if pth ~= "" && ~isempty(pth)
    if ~isfolder(pth), mkdir(pth); end
else
    pth = pwd;
end

% GSEA --------------------------------------------------------------------
if isfield(opts, "rankedGenes")
    if width(opts.rankedGenes) ~= 2
        error("rankedGenes table must have only two columns!")
    end
    
    % find gene/rank cols
    dt = varfun(@class, opts.rankedGenes, "OutputFormat", "cell");
    dt = string(dt);
    dt(dt == "cell") = "string";
    [~, idx] = ismember(["string", "double"], dt);
    if any(idx < 1)
        error("columns should be string/cell and double!")
    end
    opts.rankedGenes = opts.rankedGenes(:, idx);
    opts.rankedGenes = rmmissing(opts.rankedGenes);
    opts.rankedGenes(string(opts.rankedGenes.(1)) == "", :) = [];
    if ~isempty(opts.rankedGenes)

        % convert to py.dic()
        rankedstr = struct;
        for k = 1:width(opts.rankedGenes)
            rankedstr.("x" + k) = opts.rankedGenes.(k)';
            if isstring(rankedstr.("x" + k))
                rankedstr.("x" + k) = cellstr(rankedstr.("x" + k));
            end
        end
        opts.rankedGenes = py.dict(rankedstr);
        opts.rankedGenes = pd.DataFrame(opts.rankedGenes);
        opts.rankedGenes = opts.rankedGenes.set_index("x1");

        pre_res = gp.prerank(rnk=opts.rankedGenes, ...
                         gene_sets=py.list(cellstr(opts.lib)),...
                         threads=py.int(opts.threads),...
                         min_size=py.int(opts.min_size),...
                         max_size=py.int(opts.max_size),...
                         permutation_num=py.int(opts.permutation), ...
                         outdir=py.None, ...
                         seed=py.int(opts.seed), ...
                         verbose=opts.verbose);
        rank_res = df2tab(pre_res.res2d);
        rank_res.Name = [];
        gsets = string(pre_res.gene_sets);
        terms = rank_res.Term; % for plot
        rank_res.Gene_set(:) = "";
        for k = 1:numel(gsets)
            idx = startsWith(rank_res.Term, gsets(k) + "__");
            if any(idx)
                rank_res.Gene_set(idx) = gsets(k);
                rank_res.Term(idx) = erase(rank_res.Term(idx), textBoundary("start") + gsets(k) + "__");
            end
        end
        rank_res.Gene_set(rank_res.Gene_set == "") = "prerank";
        rank_res = movevars(rank_res, "Gene_set", Before=1);
        if ~any(rank_res.("FDR q-val") <= opts.qvalue) && opts.verbose
            fprintf("no significant term was found in GSEA!\n")
        end
        
        % only plot significant GSEA results
        if opts.plot
            tmp = rank_res;
            tmp(tmp.("FDR q-val") > opts.qvalue, :) = [];
            if ~isempty(tmp)
                % if numel(terms) < 5
                %     ax = pre_res.plot(terms=py.list(cellstr(terms)')); % , figsize=py.tuple([8,9]));
                %     fig = ax.get_figure;
                %     fig.savefig(opts.output + ".svg", dpi=400)
                % else
                    for k = 1:numel(terms)
                        pyrun("gseaplot(rank_metric=pre_res.ranking, " + ...
                            " term=terms[0], ofname='" + ...
                            replace(fullfile(pth, terms(k)), filesep, "/") + ".svg'," + ...
                            " **pre_res.results[terms[0]])", ...
                            pre_res=pre_res, ...
                            terms=cellstr(terms(k)), ...
                            gseaplot=gp.gseaplot)
                    end
                % end
            end
        end
    end
else
    rank_res = table;
end


% ORA (enrichr) -----------------------------------------------------------
if isfield(opts, "genes") && ~isempty(opts.genes)
    if ~isfield(opts, "background")
        opts.background = py.None;
    else
        idx = ismissing(opts.background) | opts.background == "";
        opts.background(idx) = [];
        if iscolumn(opts.background)
            opts.background = opts.background';
        end
        opts.background = py.list(cellstr(opts.background));
    end
    
    enr = gp.enrichr(gene_list=py.list(cellstr(opts.genes)), ...
                     gene_sets=py.list(cellstr(opts.lib)),...
                     organism=string(opts.organism), ...
                     background=opts.background, ...
                     outdir=py.None);
    enr = enr.results;
    res = df2tab(enr);
else
    res = table;
end

% write results to output
if opts.write
    if ~isempty(rank_res)
        writeGSEAtable(rank_res, opts.output + ".GSEA.xlsx", opts.qvalue)
    end

    if ~isempty(res)
        writeGSEAtable(res, opts.output + ".enrichr.xlsx", opts.qvalue)
    end
end

end % END

%% subfcuntions ===========================================================
function res = df2tab(enr)

cols = string(enr.columns.to_numpy.tolist);
vals = cell(enr.values.tolist);
if isempty(vals)
    res = table;
    return;
end
vals = cellfun(@string, vals, uni=false);
vals = vertcat(vals{:});
res = array2table(vals, VariableNames=cols);
dcols = cols(cols.lower.contains(["p-val", "odds", "score", "q-val"]) ...
    | ismember(cols, ["ES", "NES"]));
res = convertvars(res, dcols, @double);
res(:, cols(cols.lower.startsWith("old"))) = [];

end

%% ------------------------------------------------------------------------
function getGSEApyDB(gp)

lib_names = string(gp.get_library_name());
fn = readlines("GSEApy.m");
idx = find(startsWith(fn, whitespacePattern + "opts.lib"), 1);
hdr = extractBefore(fn(idx), "opts.lib") + "opts.lib";
hdr = hdr + ' {mustBeMember(opts.lib, ["' + join(lib_names, '","') + ...
    '"])} = "' + lib_names(1) + '"';
fn(idx) = hdr;
writelines(fn, which("GSEApy.m"))

end

%% ------------------------------------------------------------------------
function writeGSEAtable(res, output, qvalue)

if any(colnames(res) == "FDR q-val")
    col = "FDR q-val";
else
    col = "Adjusted P-value";
end

res(res.(col) > qvalue, :) = [];
if isempty(res), return; end

gset = unique(res.Gene_set);
for k = 1:numel(gset)
    tmp = res(res.Gene_set == gset(k), :);
    writetable(tmp, output, Sheet=matlab.lang.makeValidName(truncateStr(gset(k))));
end

end
