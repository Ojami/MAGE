function geneWrapper(qcpath, wesann, opts)
% geneWrapper creates necessary files to feed SAIGE/REGENIE gene-based
% analysis. To ease the process, the function first uses a set of
% pregenerated QC files ('vmiss', 'smiss', and 'hardy' generated from
% PLINK2 --missing --hardy midp), and exludes all monomorphic variants, and
% those with variant/sample missingness above a cut-off or HWE P below hwp
% cutoff. It's possible to exlude variants ('exlude') based on QC.
% Moreover, variants with a missingness above 'missingness' cutoff are
% excluded. Next, geneWrapper fetches the subset of these variants that are
% present in annotated ('wesann') file. This file is the output of variant
% annotation tools (see runVEP function). There are 4 different ways in
% which the annotated variants can be fetched from 'wesann' file ('model'
% option):
%   "lof": only loss-of-function variants. The current definition includes
%          only splice donors/acceptors, frameshifts and stop gained
%          variants.
%   "lof.missesne": in addition to LoF, missense variants are also fetched.
%   "vep.general": in addition to annotated LoF variants, fetches all
%                  variants present in dbNSFP (pre-processed files located
%                  in 'dbnsfpHome' directory). These variants can be
%                  filtered using 'categories' option which can be used to
%                  define multiple categories (e.g. variants with CADD > 20
%                  and REVEL > 0.75, SIFT and ClinVar pathogenic variants,
%                  ...). Alternatively (not recommended for REGENIE), user
%                  can directly use filtering arguments ('cadd', 'revel',
%                  'sift', 'polyphen', 'clinvar', 'rule') which generates
%                  only 1 single category.
%   "vep.missense": similar to "vep.general" but only filters variants that
%                   are present in annotated file. This filtering option
%                   generates a smaller set compared to "vep.general" as
%                   there exist some variants being predicated as
%                   deleterious (e.g. with SIFT or Polyphen), but are not
%                   missense (e.g. located on UTR regions).
% Following filtering variants, the function only keeps protein coding
% genes using Ensembl REST API (biotype field). Alternatively, 'geneset'
% argument can be used to explicitly determine which genes should be kept.
% Afterwards, variants with more than 1 gene will be assigned to both genes
% (e.g. is snp1 gene is geneA;geneB in the annotation file, snp1 will be
% tested twice, once for geneA and once for geneB). In the last step, all
% genes with 1 variant are removed (cannot be tested in SKAT family tests).
% However, this step should be optional in future adjustments to the
% function. 
% 
% @12JULY2024: 'categories' now must be a table with variables:
%               - name: category name, e.g. M1
%               - rule: which rule to be applied: it can be a string
%               ("or"/"and") or a numeric string corresponding to the
%               number consensus methods (e.g. "2" means at least 2 methods
%               should predict a deleterious mutation). 

%               - method: algorithms(and their filter) to be used,
%               supporting algorithms are (case insensitive) sift, revel,
%               clinvar, polyphen, cadd, PrimateAI, LRT, MutationTaster,
%               M-CAP, and AlphaMissense. These have been derived from
%               dbNSFP (see dbNSFPparser) database, and only CADD and REVEL
%               needs a numeric value cutoff. filtering criteria to be used
%               as a string scalar in paranthesis (for CADD and REVEL it
%               should be a numeric as a string scalar, e.g. "30"). The
%               rest except polyphen are ignored, as dbNSFP prediction
%               category will be used, for instance, for MutationTaster,
%               "A" and "D" are considered as deleterious mutations. For
%               polyphen, it can be either "r,relaxed" (D:probably damaging
%               or P:possibly damaging) or "s,strict" (D:probably
%               damaging). If left empty, "strict" is assumed by default.
%               This list should be comma separated.
%               EXAMPLE: categories =
%                   name      rule            method            
%               -----------------------------------------------------
%                    M1       "or"        "sift,cadd(30)"          
%                    M2       "2"   "polyphen(r),revel(0.5),m-cap"
% 
% Example:
%   geneWrapper(freq, wes, geneset=genes, exclude=qc, categories=categories)
% Oveis Jamialahmadi, University of Gothenburg, March 2022. 
% 
% @06NOV2024: start_lost and stop_lost are now considered in LoF: https://www.nature.com/articles/s41586-021-04103-z

arguments
    qcpath {mustBeFolder} % path to QC files from PLINK 
    wesann {mustBeFile, mustBeTextScalar}
    % opts.N (1,1) double % sample size used for files 'qcpath', sample size is used for missingness (if empty, uses info in acount files)
    opts.missingness (1,1) double {mustBeInRange(opts.missingness, 0, 100)} = 15 % variants/sample with missigness rate > cutoff (%) will be removed
    opts.hwp (1,1) double = 1e-15 % variants with an HWE P below this cutoff will be removed
    opts.geneset {mustBeText, mustBeVector} = "" % list of genes to keep (default: all protein coding genes). Must be in format of symbol(gene id) with no space, e.g. A1BG(ENSG00000121410)
    opts.exclude {mustBeVector, mustBeFile} % variants to exclude (e.g. those markers failing some QC filters).
    opts.model {mustBeTextScalar, mustBeMember(opts.model, ["lof", "lof.missense", "vep.general", "vep.missense"])} = "vep.general";
    opts.chr {mustBeText, mustBeVector} = "" % default: extract it from qcpath (same size of qcpath)
    opts.merge (1,1) logical = false % merge groupFiles for all chromosomes into a single file
    
    % 'vep' specific options

    % annotation categories for REGENIE gene-based tests, each element
    % should be filtering rule based on cadd/revel/polyphen/sift/clinvar
    % e.g.
    % ["(s1,or)cadd_30,revel_0.5,sift_r"
    %  "(s2,and)cadd_25&revel_0.5,sift_s,polyphen_s,clinvar"],
    % for which the first element translets into a category named as s1
    % with variants having either cadd >= 30, revel >= 0.5 or relaxed sift.
    % Second element is translated into category s2 with cadd >= 25,
    % revel >= 0.5, strict sift and polyphen, and clinvar pathogenic. Non
    % provided tools for each category are treated as non effective. 
    % This option overrides 'cadd', 'revel',
    % 'sift', 'polyphen' and 'clinvar' and 'rule' arguments.
    % 
    % @12JULY2024: 'categories' must be a table now. see the note above.
    % opts.categories {mustBeText, mustBeVector} 
    opts.categories {mustBeA(opts.categories, 'table')} 

    % mask definition for REGENIE gene-based tests. Names should be
    % according to category names in 'categories' except for "LoF" and
    % "missense", which are included in annotation file by default.
    % Example: mask = ["LoF", "LoF,missnese", "LoF,s1", "s1,s2"]. If not
    % provided, automatically will be generated based on 'categories'; e.g.
    % for category names "s1" and "s2" in 'categories', mask would be
    % ["LoF", "LoF,missense", "LoF,s1", "LoF,s2"]. both "LoF" and
    % "LoF,missense" are always included in this case. Note that, REGENIE
    % v3.0 doesn't permit duplicates in annotation files (i.e. no
    % overlapping categories). Thus, it's important that mask file is
    % defined as needed. For instance a mask defined as ["LoF",
    % "LoF,missnese", "LoF,CADD25"] is not specific for "missense" since
    % there is an overlap beteween "CADD25" and "missense". The proper
    % inclusive mask in this case would be ["LoF", "LoF,missnese,CADD25",
    % "LoF,CADD25"], so that "LoF,missense,CADD25" covers all missense
    % variants. This issue doesn't exist when 'mask' option is left empty,
    % and geneWrapper automatically defines 'mask'.
    opts.mask {mustBeText, mustBeVector} 
    opts.cadd (1,1) double = 0; % filter variants with CADD-phred above this
    opts.revel (1,1) double = 0; % filter variants with REVEL above this
    opts.sift {mustBeMember(opts.sift, ["strict", "relaxed", "none"])} = "none";
    opts.polyphen {mustBeMember(opts.polyphen, ["strict", "relaxed", "none"])} = "none"; % probably_damaging is "D" and possibly_damaging "P" in dbNSFP 
    opts.clinvar (1,1) logical = false; % if true, matches ClinVar flags of "Likely_pathogenic" and "pathogenic"
    opts.rule {mustBeMember(opts.rule, ["and", "or"])} = "and"; % how to merge filters
    opts.dbnsfpHome {mustBeTextScalar} = fullfile(fileparts(which("geneWrapper.m")), "step2.GENE", "annotation", "dbNSFP") 

    opts.tool {mustBeTextScalar, mustBeMember(opts.tool, ["saige", "regenie"])} = "regenie" % output will be in required format of either SAIGE-GENE+ or REGENIE v3+
    opts.debug (1,1) logical = false % show full details
end

% parse QC files 
filter_ids = parseWESqc(qcpath, opts);

if opts.tool == "saige"
    annfile = @(x) sprintf("geneWrapper.groupFile.%s.txt", x);
elseif opts.tool == "regenie"
    annfile = @(x) sprintf("geneWrapper.annotation.%s.txt", x);
    setfile = @(x) sprintf("geneWrapper.set.%s.txt", x);
    maskfile = @(x) sprintf("geneWrapper.mask.%s.txt", x);
    if opts.merge % not supported for regenie due to the possibility of having non-overlapping sets
        opts.merge = false;
        fprintf('warning:with tool as regenie, merge flag is ignored!\n')
    end
end

% loop over qc files
if all(opts.chr == "")
    opts.chr = filter_ids.chr;
end
tic
[groupfile, setlist, mask] = deal(cell(numel(filter_ids.chr), 1));
for i = 1:numel(filter_ids.chr)
    fprintf('processing file %d of %d\n', i, numel(filter_ids.chr))
    optsin = opts;
    optsin.chr = opts.chr(i);

    % if numel exclude == 1, 1 file exsits for all variants (not
    % recommended here, this can instead be used directly in REGENIE)
    if isfield(opts, 'exclude') && numel(opts.exclude) > 1
        optsin.exclude = opts.exclude(i);
    end

    if opts.merge
        [groupfile{i}, setlist{i}, mask{i}] = innerWrapperSingle(filter_ids.filter{i}, wesann, optsin);
    else
        [groupfile{1}, setlist{1}, mask{1}] = innerWrapperSingle(filter_ids.filter{i}, wesann, optsin);
        if opts.tool == "saige"
            writecell(groupfile{1}, annfile("chr" + i), 'FileType', 'text', 'QuoteStrings', false);
        elseif opts.tool == "regenie"
            tag = "chr" + i;
            writetable(groupfile{1}, annfile(tag), 'Delimiter', '\t', 'WriteVariableNames', false)
            writetable(setlist{1}, setfile(tag), 'Delimiter', '\t', 'WriteVariableNames', false)
            writematrix(mask{1}, maskfile(tag), 'FileType', 'text', 'QuoteStrings', false, 'Delimiter', '\t')
        end
    end
    fprintf('================================================\n\n')
end

if opts.merge
    groupfile(cellfun(@isempty, groupfile)) = [];
    setlist(cellfun(@isempty, setlist)) = [];
    groupfile = vertcat(groupfile{:});
    setlist = vertcat(setlist{:});

    if opts.tool == "saige"
        fprintf('Writing %d genes to %s\n', numel(groupfile), annfile)
        writecell(groupfile, annfile("merged"), 'FileType', 'text', 'QuoteStrings', false);
        fprintf('Done!\n')
    elseif opts.tool == "regenie"
        writetable(groupfile, annfile("merged"), 'Delimiter', '\t', 'WriteVariableNames', false)
        writetable(setlist, setfile("merged"), 'Delimiter', '\t', 'WriteVariableNames', false)
    end
end
toc

end % END

%% subfunctions ===========================================================
function [wesann, setlist, mask] = innerWrapperSingle(filter_id, wesannFile, opts)

% DEPRECATED
% % read freq (.count) file generated with plink --freq count 
% [~, ~, fileExt] = fileparts(countfile);
% fprintf('reading and processing count file %s...', countfile)
% countfile = datastore(countfile,'VariableNamingRule', 'preserve',...
%     'Type', 'tabulartext','FileExtensions', fileExt, 'TextType', 'string');
% countfile.SelectedVariableNames = {'ID', 'REF', 'ALT', 'ALT_CTS', 'OBS_CT'};
% countfile.SelectedFormats = {'%s', '%s', '%s', '%f', '%f'};
% countfile = readall(countfile);
% fprintf('\b\b Done (%d variants).\n', height(countfile))

if isfield(opts, 'exclude') % exclude these varaints 
    exvariants = readmatrix(opts.exclude, 'FileType', 'text', ...
        'NumHeaderLines', 0, 'ExpectedNumVariables', 1, ...
        'OutputType', 'string', 'Delimiter', '');
    rmidx = ismember(filter_id, exvariants);
    if any(rmidx)
        % fprintf('%d variants were excluded due to exclusion criteria.\n',
        % sum(rmidx)) %@17OCT2024
        filter_id = union(exvariants, filter_id);
    end
end

% % Remove variants with AAC = 0
% aac_zero = filter_id.ALT_CTS < 1;
% if any(aac_zero)
%     fprintf('%d monomorphic variants were removed (!AAC).\n', sum(aac_zero))
%     filter_id(aac_zero, :) = [];
% end
% 
% % remove variants with high missingness
% if isfield(opts, 'N')
%     n = opts.N;
%     if (opts.N > max(filter_id.OBS_CT./2))
%         warning('off', 'backtrace')
%         warning('%s does not correspond to count files!', opts.N)
%         warning('on', 'backtrace')
% 
%         n = max(filter_id.OBS_CT./2); % approximate sample size
%     end
% else
%     n = max(filter_id.OBS_CT./2); % approximate sample size
% end
% missingness = (n - filter_id.OBS_CT./2)./n;
% remidx = missingness >= (opts.missingness./100);
% if any(remidx)
%     fprintf('%d variants were removed due to missingness > %.0f%%.\n', sum(remidx), opts.missingness)
%     filter_id(remidx, :) = [];
% end

% Check variants with MAF > 0.5. For these variants ref and alt alleles
% should be swapped. 
% ref2Alt_ratio = countfile.ALT_CTS./countfile.OBS_CT;
% filter_id.ALT_CTS = []; filter_id.OBS_CT = [];
% ref2Alt_swap = ref2Alt_ratio > 0.5;

% SAIGE does flipping internally 
% if sum(ref2Alt_swap)
%     fprintf('alternate (minor) allele for %d variants will be flipped...\n', sum(ref2Alt_swap))
%     swap_ref_allele = ref_allele(ref2Alt_swap);
%     swap_alt_allele = alt_allele(ref2Alt_swap);
%     ref_allele(ref2Alt_swap) = swap_alt_allele;
%     alt_allele(ref2Alt_swap) = swap_ref_allele;
% end

% if opts.chr == ""
%     opts.chr = string(extractBefore(filter_id.ID(1), ':'));
% end

% reading annotations -----------------------------------------------------

% TODO: can be set as an optional argument 
% LOF filters: https://www.nature.com/arles/s41586-020-2853-0
% A more stringent definition does not contain start_/stop_lost as
% the one used in LOFTEE
% @06NOV2024: start_lost and stop_lost are now considered in LoF: https://www.nature.com/articles/s41586-021-04103-z

disp("@06NOV2024: start_lost and stop_lost are now considered in LoF: https://www.nature.com/articles/s41586-021-04103-z")
lof_filters = ["stop_gained"; "splice_donor_variant"; ...
        "frameshift_variant"; "splice_acceptor_variant"; ...
        "start_lost"; "stop_lost"];

fprintf('scheme: %s\nfiltering using dbNSFP annotations...\n', opts.model)
bim.chr = opts.chr;
if ismember(opts.model, ["vep.general", "lof"])
    all_filters = lof_filters;
elseif opts.model == "lof.missense"
    % only LoF + missense (no inframe indels) --> shoud be deprecated since
    % vep.general/vep.missense with 'categorical' option creates flexible
    % gene sets
    all_filters = union(lof_filters, "missense_variant");
else % "vep.missense": all variants with HIGH/MODERATE impact flag in Ensembl VEP
    all_filters = ["stop_gained"; "start_lost"; "splice_donor_variant"; ...
        "frameshift_variant"; "splice_acceptor_variant"; "stop_lost";...
        "transcript_ablation"; "missense_variant"; ...
        "protein_altering_variant"; "inframe_insertion"; "inframe_deletion"];
end

opts.chr = "^" + opts.chr + ":";
all_filters = join(all_filters, '|');
extCol = ["Uploaded_variation", "SYMBOL", "Gene", "Consequence"];
% wesann = bfilereader(wesannFile, 'pattern', [opts.chr, all_filters],...
%     'patternCol', ["Uploaded_variation", "Consequence"],...
%     'multiCol', true, 'extractCol', extCol,...
%     'parallel', true, 'verbose', 'off', 'header', true);

%@07APR2025: a bug in latest bfilereader with parallel flag makes it
%slower.
wesann = bfilereader(wesannFile, 'pattern', [opts.chr, all_filters],...
    'patternCol', ["Uploaded_variation", "Consequence"],...
    'multiCol', true, 'extractCol', extCol,...
    'parallel', false, 'verbose', 'off', 'header', true);

colnames = wesann.Properties.VariableNames;
[~, colnames_idx] = ismember(extCol, colnames);
wesann = wesann(:, colnames_idx);
wesann.Properties.VariableNames = {'id', 'gene', 'geneid', 'consequence'};
wesann.gene = strrep(wesann.gene, '.;', ''); % a bug in variant annotator
wesann.geneid = strrep(wesann.geneid, '.;', '');

if opts.debug
    cons = unique(wesann.consequence).join(newline);
    fprintf('consequence summary in annotated file:\n%s\n', cons)
    fprintf('\n')
end

if contains(opts.model, "vep.") % filter rest of variants using functional prediction tools via dbNSFP dataset
    if opts.model == "vep.general"
        bim.snp = []; % filter all variants in dbNSFP dataset
    else
        lofidx = contains(wesann.consequence, lof_filters);
        bim.snp = wesann.id(~lofidx); % filter non-LoF, i.e. missense variants
        cols = ["consequence", "gene", "geneid"];
        for k = 1:numel(cols)
            bim.(cols(k)) = wesann.(cols(k))(~lofidx);
        end
        wesann(~lofidx, :) = []; % keep LoF variants
    end
    
    % check if categories should be created for different filtering schemea
    if isfield(opts, 'categories')
        [vep, opts] = parseAnnotatedCategories(bim, opts);
    
    else % two categories: missense and LoF
        vep = callVEP(bim, 'method', 'dbnsfp', 'readall', false, ...
            'filter', true, 'cadd', opts.cadd, 'revel', opts.revel,...
            'clinvar', opts.clinvar, 'sift', opts.sift, 'polyphen', ...
            opts.polyphen, 'rule', opts.rule,...
            'dbnsfpHome', opts.dbnsfpHome, 'verbose', false);
        vep = vep(:, [1, 3, 4]);
        vep.Properties.VariableNames = {'id', 'gene', 'geneid'};
    end
    vep = unique(vep, "rows", "stable");
    
    % add consequence column to vep table 
    if isempty(bim.snp) % fetch annotations for these ids
        get_filters = ["stop_gained"; "start_lost"; "splice_donor_variant"; ...
            "frameshift_variant"; "splice_acceptor_variant"; "stop_lost";...
            "transcript_ablation"; "missense_variant"; ...
            "protein_altering_variant"; "inframe_insertion"; "inframe_deletion"];
        get_filters = setdiff(get_filters, lof_filters);
        get_filters = join(get_filters, '|');
        extCol = ["Uploaded_variation", "SYMBOL", "Gene", "Consequence"];
        bim = bfilereader(wesannFile, 'pattern', [opts.chr, get_filters],...
            'patternCol', ["Uploaded_variation", "Consequence"],...
            'multiCol', true, 'extractCol', extCol,...
            'parallel', false, 'verbose', 'off', 'header', true);
        colnames = bim.Properties.VariableNames;
        [~, colnames_idx] = ismember(extCol, colnames);
        bim = bim(:, colnames_idx);
        bim.Properties.VariableNames = {'snp', 'gene', 'geneid', 'consequence'};
        bim.gene = strrep(bim.gene, '.;', ''); % a bug in variant annotator
        bim.geneid = strrep(bim.geneid, '.;', '');
    end
    
    [idx1, idx2] = ismember(vep.id, bim.snp); idx2(idx2 < 1) = [];
    vep.consequence = strings(height(vep), 1);
    vep.consequence(idx1) = bim.consequence(idx2);
    vep.gene(idx1) = bim.gene(idx2);
    vep.geneid(idx1) = bim.geneid(idx2);

    if opts.debug
        cons = unique(vep.consequence).join(newline);
        fprintf('consequence summary in dbNSFP:%s\n', cons)
        fprintf('\n')
    end

    % append those missense variants in bim (i.e. in annotated dataset)
    % that are absent from vep (i.e. dbNSFP database) to vep table. This
    % can happen only with "vep.general" or "vep.missense"
    if isstruct(bim)
        bim = rmfield(bim, "chr");
        bim = struct2table(bim);
    end
    bim(ismember(bim.snp, vep.id), :) = [];
    bim(~contains(bim.consequence, "missense_variant"), :) = [];
    if ~isempty(bim.snp)
        if opts.debug
            fprintf('%d missense variants absent from dbNSFP but present in annotated file were added\n', numel(bim.snp))
        end
        bim.categories = strings(height(bim), 1);
        snpCol = ismember(bim.Properties.VariableNames, "snp");
        bim.Properties.VariableNames(snpCol) = "id";
        bim = bim(:, vep.Properties.VariableNames);
        vep = [vep; bim];
    end
    
    % add categories (missense nad LoF are added automatically)
    lof_idx = contains(wesann.consequence, lof_filters);
    wesann.consequence(lof_idx) = "LoF";
    wesann.consequence(~lof_idx) = ""; % to be excluded
    
    miss_idx = contains(vep.consequence, "missense_variant");
    vep.consequence(miss_idx) = "missense";
    vep.consequence(~miss_idx) = ""; % to be excluded (e.g. inframe indels)
    
    % add other categories is present
    if isfield(opts, 'categories')
        vep.consequence = vep.consequence + "," + vep.categories;
        vep.consequence = regexprep(vep.consequence, '(^,$|,$|^,)', '');
        vep.categories = [];
    end

    vep = unique(vep, "rows", "stable");
    
    %@22JUNE2023. we no longer merge overlapping variants from VEP (see
    % runVEP) since this may cause wrong (although only for few variants)
    % functional (e.g. LoF) assignments at later steps. As example
    % "1:151067336:T:C" variant overlaps CDC42SE1 and MLLT11, one is
    % splice_acceptor_variant and the other is missense_variant. So, it
    % creates a dilemma for later tools to dissect which is which (are both
    % LoF if the tool decides to keep the most severe one? no!). Hence, we
    % leave them as they are (i.e. those variants are repeated in the final
    % file, once for each gene). However, we merge multiple consequences
    % (due to VEP annotations from different trnascripts) for each variant.
    % 
    % So, to avoid deleting these overlapping variants (e.g. they are both
    % in vep and wesann tables, one is LoF the other one is missense, but
    % they lie on different genes, and therefore should be kept). Note, we
    % should reconcile (and expand if necessary) variants in 'vep' that
    % overlap multiple genes. As noted above, with the new update of runVEP
    % function, they are already expanded, but variants in 'vep' table come
    % from dbNSFP, so there may exist such variants for with multiple
    % gene/geneid (separated with ; or ,). However, such variants, if
    % overlapping with 'wesann' variants, should be exluded (i.e. if not in
    % wessannFile missense category, then we assume they're LoF because
    % they appear in 'wesann', hence one appearance is enough).
    vep.t1 = vep.id + ":" + vep.gene;
    vep.t2 = vep.id + ":" + vep.geneid;
    wesann.t1 = wesann.id + ":" + wesann.gene;
    wesann.t2 = wesann.id + ":" + wesann.geneid;
    multi_idx = vep.gene.contains([";", ","]) | vep.geneid.contains([";", ","]);
    rmidx = ismember(vep.t1, wesann.t1) | ismember(vep.t2, wesann.t2) ...
        | (multi_idx & ismember(vep.id, wesann.id));
    vep(rmidx, :) = [];
    vep(:, ["t1", "t2"]) = [];
    wesann(:, ["t1", "t2"]) = [];

    % [f1, f2] = ismember(wesann.id, vep.id); f2(f2<1) = [];
    % vtab = vep(f2, :);
    % wtab = wesann(f1, :);
    % vtab.t1 = vtab.id + ":" + vtab.gene;
    % vtab.t2 = vtab.id + ":" + vtab.geneid;
    % wtab.t1 = wtab.id + ":" + wtab.gene;
    % wtab.t2 = wtab.id + ":" + wtab.geneid;
    % multi_idx = vtab.gene.contains([";", ","]) | vtab.geneid.contains([";", ","]);
    % rmidx = ismember(vtab.t1, wtab.t1) | ismember(vtab.t2, wtab.t2) | multi_idx;

    % vep(ismember(vep.id, wesann.id), :) = [];
    wesann = [wesann; vep];
    clear vep
else % "lof" and "lof.missense"
    lof_idx = contains(wesann.consequence, lof_filters);
    miss_idx = contains(wesann.consequence, "missense_variant");
    wesann.consequence(miss_idx) = "missense";
    wesann.consequence(lof_idx) = "LoF"; % overwrite missense+lof
    wesann.consequence(~(miss_idx | lof_idx)) = "";
end

if opts.debug
    fprintf('%d variants after merging annotated and filtering\n', height(wesann))
    fprintf('%d variants are not LoF, missense or filtered\n', sum(wesann.consequence == ""))
end
wesann(wesann.consequence == "", :) = []; % exclude other variants (e.g. inframe indels)

% remove variants failed QC (external, missingness and zero AAC) ----------
fprintf('removing non QCed variants...')
% [id_idx1, id_idx2] = ismember(wesann.id, filter_id);
% wesann = wesann(id_idx1, :);
% filter_id = filter_id(id_idx2(id_idx1));

%@17OCT2024: due to some changes, ref/alt are extracted from wesann ids
%directly
id_idx1 = ismember(wesann.id, filter_id);
wesann(id_idx1, :) = [];
tmp_id = split(wesann.id, ":");
wesann.ref = tmp_id(:, 3);
wesann.alt = tmp_id(:, 4);
clear tmp_id filter_id

% wesann.ref = filter_id.REF;
% wesann.alt = filter_id.ALT;
% if ~all(wesann.id == filter_id)
%     error('Cannot match annotated and allele count files!')
% else
%     clear filter_id
% end
fprintf('\b\b Done (%d variants failing QC).\n', sum(id_idx1))

% checking biotype of genes: only keep protein_coding ---------------------
if isempty(opts.geneset) || any(opts.geneset == ""); opts.geneset = []; end
geneid = unique(wesann.geneid);
geneid = arrayfun(@(x)split(x, ';'), geneid, 'uni', false);
geneid = unique(vertcat(geneid{:}));

if isempty(opts.geneset) % otherwise keep genes in geneset irrespective of biotypess
    fprintf('fetching gene biotypes from Ensembl REST...\n')
    steps = unique([1:900:numel(geneid), numel(geneid) + 1]);
    if steps(end)-steps(end-1) < 2 % uses GET, so output is struct
        % shift steps
        steps(end-1) = steps(end-1) - 1;
    end
    raw = cell(numel(steps)-1, 1);
    for i = 1:numel(steps)-1
        rawtmp = EnsemblREST(geneid(steps(i):steps(i+1)-1), 'lookup', 'refGenome', '38', 'parseFlag', false);
        rawtmp = struct2cell(rawtmp);
        rawtmp(cellfun(@isempty, rawtmp)) = [];
        raw{i} = rawtmp;
        clear rawtmp
    end
    raw = vertcat(raw{:});
    
    biotype = strings(numel(raw), 3);
    for i = 1:numel(raw)
        if isfield(raw{i}, 'display_name')
            biotype(i, :) = string({raw{i}.id, raw{i}.biotype, raw{i}.display_name});
        else
            biotype(i, :) = string({raw{i}.id, raw{i}.biotype, missing});
        end
    end
    clear raw
    
    % check remaining gene ids with genome version 37
    geneid = setdiff(geneid, biotype(:, 1));
    steps = unique([1:900:numel(geneid), numel(geneid) + 1]);
    if ~isempty(geneid)
        biotyper = cell(numel(steps) - 1, 1);
        for i = 1:numel(steps)-1
            geneRange = geneid(steps(i):steps(i+1)-1);
            if numel(geneRange) < 2
                geneRange(2) = geneRange;
            end
            biotyper{i} = EnsemblREST(geneRange, 'lookup', 'refGenome', '37');
        end
        biotyper(cellfun(@isempty, biotyper)) = [];
        biotyper = vertcat(biotyper{:});
        if ~isempty(biotyper) && ~istable(biotyper)
            biotyper = cellfun(@(x){x.id, x.biotype}, biotyper, 'uni', false);
            biotyper = string(vertcat(biotyper{:}));
            biotype = [biotype; [biotyper, strings(size(biotyper, 1), 1)]]; % don't get gene names from '37'
        elseif istable(biotyper)
            biotype = [biotype; [biotyper.id, biotyper.biotype, biotyper.name]];
        end
    end
    removeGenes = setdiff(geneid, biotype(:, 1)); % remaining genes (not available in Ensembl 37 or 38) so should be removed
    nonCodingIdx = biotype(:, 2) ~= "protein_coding";
    removeGenes = union(removeGenes, biotype(nonCodingIdx, 1));
    biotype(nonCodingIdx, :) = [];
    fprintf('\b\b Done.\n')

else

    opts.geneset = regexp(opts.geneset, '(.*)[(](.*)[)]', 'tokens');
    rmidx = cellfun(@isempty, opts.geneset); % those with only symbol or gene id
    opts.geneset(rmidx, :) = [];
    opts.geneset = vertcat(opts.geneset{:});
    opts.geneset = vertcat(opts.geneset{:});
    idx = ismember(opts.geneset(:, 2), geneid);
    opts.geneset(~idx, :) = []; % only keep genes present in WES dataset
    removeGenes = setdiff(geneid, opts.geneset(:, 2)); % exclude genes not present in geneset
    biotype = [opts.geneset(:, 2), ...
        repmat("protein_coding", size(opts.geneset, 1), 1), ...
        opts.geneset(:, 1)]; % gene id | biotype | gene symbol
end

if ~isempty(removeGenes)
    fprintf('%d deprecated or non protein coding genes will be removed.\n', numel(removeGenes))
    % one issue here is overlapping variants for which joined symbols and
    % ids do not correspond, i.e. symbol for variant X maybe A1;B1, but id
    % maybe Bid;Aid.
    % First try to identify unique (non-overlapping) variants to be removed
    rem_idx = ismember(wesann.geneid, removeGenes);
    wesann(rem_idx, :) = []; % remove these unique variants 
    
    % Next, check what has been left: overlapping variants 
    % these genes and their variants will be removed after reconciling
    % overlapping variants in next step.
end

wesann.gene = []; % remove gene names. Working with ENSG ids is less ambiguous.
% Find overlapping genes
f_overlap = contains(wesann.geneid, ';');
fprintf('%d overlapping vairants found.\n', sum(f_overlap))
fprintf('expanding overlapped variants...')

overlapped_genes = wesann(f_overlap, :);
wesann(f_overlap, :) = [];

overlapped_genes_unique = unique(overlapped_genes.geneid);
[fixed_genes, fixed_id, fixed_cons, fixed_ALT, fixed_REF] = ...
    deal(cell(numel(overlapped_genes_unique), 1));

for ii = 1:numel(overlapped_genes_unique)
    broken_genes = split(overlapped_genes_unique(ii), ';');
    overlap_idx = overlapped_genes.geneid == overlapped_genes_unique(ii);
    overlapping_id = overlapped_genes.id(overlap_idx);
    overlapping_category = overlapped_genes.consequence(overlap_idx);
    overlapping_alt = overlapped_genes.alt(overlap_idx);
    overlapping_ref = overlapped_genes.ref(overlap_idx);
    fixed_genes{ii} = sort(repmat(broken_genes, numel(overlapping_id), 1));
    fixed_id{ii} = repmat(overlapping_id, numel(broken_genes), 1);
    fixed_cons{ii} = repmat(overlapping_category, numel(broken_genes), 1);
    fixed_ALT{ii} = repmat(overlapping_alt, numel(broken_genes), 1);
    fixed_REF{ii} = repmat(overlapping_ref, numel(broken_genes), 1);
end

fixed_tab = table(fixed_id, fixed_genes, fixed_cons, fixed_REF, fixed_ALT);
fixed_tab = varfun(@(x)vertcat(x{:}), fixed_tab);
fixed_tab.Properties.VariableNames = {'id', 'geneid', 'consequence', 'ref', 'alt'};

wesann = [wesann; fixed_tab];
clear fixed_tab fixed_id fixed_genes fixed_cons fixed_REF fixed_ALT
fprintf('\b\b Done.\n')

% Remove empty gene ids
f_dot = ismember(wesann.geneid, [".", " "]);
wesann(f_dot, :) = [];

% check if any genes left to be removed
if ~isempty(removeGenes)
    wesann(ismember(wesann.geneid, removeGenes), :) = [];
end

% get final gene symbols --------------------------------------------------
fprintf('updating gene names...')
biotype(:, 2) = []; % don't need biotype anymore
missingGenes = biotype(ismissing(biotype(:, 2)), 1);
if ~isempty(missingGenes)
    % Ensembl REST doesn't contain gene names for these ids, use biotools.fr API instead
    missingGenes = EnsemblREST(missingGenes, 'biotoolsfr'); 
    if ~isempty(missingGenes)
        [~, idx] = ismember(missingGenes.id, biotype(:, 1)); 
        biotype(idx, 2) = missingGenes.name;
    end
end
% replace missing gene names
missingidx = ismissing(biotype(:, 2));
biotype(missingidx, 2) = "";
[idx1, idx2] = ismember(wesann.geneid, biotype(:, 1));
wesann.gene = strings(size(wesann, 1), 1);
wesann.gene(idx1) = biotype(idx2, 2);
if opts.tool == "regenie"
    wesann.gene_list = wesann.gene + "(" + wesann.geneid + ")";
else
    wesann.gene_list = wesann.gene + "_" + wesann.geneid;
    wesann.gene_list = regexprep(wesann.gene_list, '^_', '');
end
wesann.gene = []; wesann.geneid = [];
fprintf('\b\b Done.\n')

% Remove genes with 1 variant (2nd round after removing uncalled variants)
[Geneunique, ~, gene_idx] = unique(wesann.gene_list, 'stable');
count = histcounts(gene_idx,[1:numel(Geneunique), inf]);
remove_singles = ismember(wesann.gene_list, Geneunique(count == 1));
if any(remove_singles)
    fprintf('%d genes with 1 variant were removed \n', sum(remove_singles))
    wesann(remove_singles, :) = [];
end

% sort variants 
wesann.pos = regexp(wesann.id, ".+?[:](.+?)[:]", 'tokens');
wesann.pos = double(string(vertcat(wesann.pos{:})));
[~, idx] = sort(wesann.pos);
wesann = wesann(idx, :);

[mask, setlist] = deal([]); % only for REGENIE
%--------------------------------------------------------------------------
if opts.tool == "saige"
    fprintf('updating variant ids for SAIGE-GENE GroupFile input...')

    % Remove REF:ALT from variant IDs
    wesann.id = regexprep(wesann.id,'(\w+|\d+)[:](\d+)[:](\w+)[:](\w+)', '$1:$2');
    wesann.id = regexprep(wesann.id, '[.].*', '');
    
    if any(contains(wesann.id, '.'))
        fprintf('ERROR: fractional part error!\n')
        return
    end
    
    % Create variant IDs for group file
    wesann.id = strcat(wesann.id, '_', wesann.ref, '/', wesann.alt);
    fprintf('\b\b Done.\n')
    
    %----------------------------------------------------------------------
    fprintf('grouping variants for each gene...\n')
    wesann(:, ["alt", "ref", "consequence", "pos"]) = [];
    wesann = groupsummary(wesann, 'gene_list', @(x)strjoin(x, '\t'), 'id');
    wesann = compose('%s\t%s', wesann.gene_list, wesann.(3));

else % REGENIE
    % REGENIE doesn't allow duplicates in annotation file, i.e. only one
    % category can be assigned to each variant (note that each variant can
    % be assigned to > 1 gene). Therefore, in creating mask, all
    % overlapping categories can be included to have a full annotation
    % category. Example: to have LoF,missense mask when having another
    % category of CADD25 (i.e. a subcategory of missense), one should
    % include LoF,missense,CADD25 to cover all missense variants. Note that
    % this only applies to either mutually exclusive categories or when a
    % category is subset of another category (e.g. missense variants with
    % CADD>25 are a subset of missense). Otherwise, for each pair of
    % overlapping categories, only the bigger (with more varinats) category
    % is kept and the other one is removed.

    wesann(:, ["alt", "ref"]) = [];

    % find gene pos: mean position of first/last variants 
    gpos = groupsummary(wesann, 'gene_list', @(x)round((min(x)+max(x))/2), 'pos');
    [idx1, idx2] = ismember(wesann.gene_list, gpos.(1));
    wesann.pos(idx1) = gpos.(3)(idx2(idx1));


    % create set list file
    fprintf('creating set list...')
    setlist = groupsummary(wesann, 'gene_list', @(x)strjoin(x, ','), 'id');
    setlist.GroupCount = [];
    setlist.Properties.VariableNames{2} = 'id'; 
    setlist.chr = repmat(erase(opts.chr, ["^", ":"]), height(setlist), 1);
    [~, idx2] = ismember(setlist.gene_list, wesann.gene_list);
    setlist.pos = wesann.pos(idx2);
    setlist = setlist(:, ["gene_list", "chr", "pos", "id"]);
    setlist = sortrows(setlist, 'pos');
    fprintf('\b\b Done.\n')

    % prepare annotation file
    fprintf('creating annoation file...')    
    wesann.pos = [];
%     wesann.consequence = replace(wesann.consequence, ["&", "|"], ["AND", "OR"]);
    idx = wesann.consequence.contains(","); 
    cats = arrayfun(@(x) split(x, ','), unique(wesann.consequence(idx)), 'uni', false); % count only overlapping categories (i.e. those with ',')
    cats = unique(vertcat(cats{:}));
    
    for i = 1:numel(cats)
        wesann.(cats(i)) = contains(wesann.consequence, ...
            regexpPattern("(,|^)" + cats(i) + "(,|$)"));
    end
    
    wesann.consequence = replace(wesann.consequence, ",", "_");

    % grps = arrayfun(@(x)split(x, ","), wesann.consequence, 'uni', false);
    % count = wesann.consequence.count(",") + 1;
    % wesann = varfun(@(m)repelem(m, count), wesann);
    % wesann.Properties.VariableNames = regexprep(wesann.Properties.VariableNames, '^Fun_', '');
    % wesann.consequence = vertcat(grps{:});
    % wesann = wesann(:, ["id", "gene_list", "consequence"]);

    % define masks 
    if isfield(opts, 'mask') && numel(wesann) < 2 % if there are different variant sets, override 'mask' option
        mask = opts.mask;
        mask_names = "mask" + (1:numel(mask))';
    else
        % see desc above under annotation file
        mask = cats;
        mask_names = ["mask_lof"; "mask_lof_" + mask];
        for j = 1:numel(mask) % create inclusive masks: all overlapping masks
            mask(j) = join(unique(wesann.consequence(wesann.(mask(j)))), ",");
        end
        mask = ["LoF"; "LoF," + mask];
    end

    if isrow(mask); mask = mask'; end
    mask = [mask_names, mask];

    % rearrange annotation table
    wesann = wesann(:, ["id", "gene_list", "consequence"]);
    fprintf('\b\b Done.\n\n')

end

end % END

%% Subfunctions
function [vep, opts] = parseAnnotatedCategories(bim, opts)
%@12JULY2024: 'categories' now must be a table (see notes above on how to
%prepare this table, or see call_geneWrapper to find examples)

% opts.categories = string(opts.categories);
% opts.categories = replace(opts.categories, whitespacePattern, '');
% opts.categories = regexp(opts.categories, '[(](.+?)[,](.+?)[)](.+)', 'tokens');
% opts.categories = vertcat(opts.categories{:});
% opts.categories = vertcat(opts.categories{:});
% opts.categories = lower(opts.categories);

cols = colnames(opts.categories);
assert(all(ismember( ["name", "rule", "method"], cols)), "geneWrapper:'name', 'rule' and 'method' columns are needed!")

% case insensitive for rule and method
opts.categories.rule = string(opts.categories.rule);
opts.categories = convertvars(opts.categories, ["method", "rule"], @lower);

%@12JULY2024: check rules
idx = opts.categories.rule == "and";
if any(idx)
    opts.categories.rule(idx) = string(opts.categories.method(idx).count(",") + 1); % all methods 
end
idx = opts.categories.rule == "or";
if any(idx)
    opts.categories.rule(idx) = "1"; % all methods 
end

opts.categories.rule = double(opts.categories.rule);
assert(all(~isnan(opts.categories.rule)), "geneWrapper:only 'and'/'or' or numeric values are allowed for filtering rules!")

% check for method/rule duplicates and throw a warning 
tmpid = opts.categories.method + ":" + opts.categories.rule;
tmpid = duplicates(tmpid);
if ~isempty(tmpid)
    fprintf("\tWARNING:there are duplicates in method/rule definitions, keeping only 1!\n")
    [~, idx] = unique(opts.categories(:, ["method", "rule"]), "rows", "stable");
    opts.categories = opts.categories(idx, :);
end

% if ~all(ismember(opts.categories(:, 2), ["and", "or"]))
%     error('geneWrapper:only "and"/"or" is allowed for filtering rules!')
% end

% first get full variants in dbNSFP dataset
vep = callVEP(bim, 'method', 'dbnsfp', 'readall', false, ...
        'filter', false, 'dbnsfpHome', opts.dbnsfpHome, 'verbose', false);
vep = vep(:, [1, 3, 4]);
vep.Properties.VariableNames = {'id', 'gene', 'geneid'};

% sz = size(opts.categories, 1);
sz = height(opts.categories);
for k = 1:sz
    % filtering = opts.categories(k, :);
    filfor.name = split(opts.categories.method(k), ",");
    
    %@12JULY2024
    idx = contains(filfor.name, ["(", ")"]);
    filter_vals = extractBetween(filfor.name(idx), "(", ")");
    filfor.name(idx) = erase(filfor.name(idx), "(" + filter_vals + ")");
    filfor.val = strings(numel(filfor.name), 1);
    filfor.val(idx) = filter_vals;
    % for sift and polyphen only
    spCheck = ["sift", "polyphen"];
    for j = 1:numel(spCheck)
        idx = filfor.name == spCheck(j);
        if any(idx)
            if filfor.val(idx).startsWith("r")
                filfor.val(idx) = "relaxed";
            else
                filfor.val(idx) = "strict";
            end
        end
    end

    filfor.name = replace(filfor.name, "-", "_"); % for m-cap
    
    % find double cutoffs
    filfor.filter = double(~ismissing(double(filfor.val)));

    % find named arguments 
    idx = ~filfor.filter & filfor.val ~= "" & ~ismissing(filfor.val);
    filfor.filter(idx) = 2;


    % filtools.name = ["cadd", "revel", "sift", "polyphen", "clinvar"];
    % filtools.val = ["double", "double", "strict,relaxed,none", ...
    %     "strict,relaxed,none", "logical"];
    % for l = 1:numel(filtools.name)
    %     ispresent = startsWith(filfor.name, filtools.name(l));
    % 
    %     switch filtools.val(l)
    %         case 'double'
    %             filfor.(filtools.name(l)) = 0;
    %             if any(ispresent)
    %                 filfor.(filtools.name(l)) = filfor.name(ispresent).extractAfter("_").double;
    %             end
    %         case "logical"
    %             filfor.(filtools.name(l)) = false;
    %             if any(ispresent); filfor.(filtools.name(l)) = true; end
    %         otherwise % categorical
    %             filfor.(filtools.name(l)) = "none";
    %             if any(ispresent)
    %                 filfor.(filtools.name(l)) = filfor.name(ispresent).extractAfter("_");
    %                 filtools_val = filtools.val(l).split(",");
    %                 ispresent = startsWith(filtools_val, filfor.(filtools.name(l)));
    %                 filfor.(filtools.name(l)) = string(filtools_val(ispresent));
    %             end
    %     end
    % 
    %     opts.filterset.(filtools.name(l))(k) = filfor.(filtools.name(l));
    % end
    % 
    opts.filterset.tag(k) = opts.categories.name(k);
    opts.filterset.rule(k) = opts.categories.rule(k);

    if ismember(opts.filterset.tag(k), ["lof", "missense"])
        fprintf('ignoring filter "%s", since it''s automatically added by geneWrapper!\n', opts.filterset.tag(k))
        opts.filterset.tag(k) = "";
        continue
    end

    % for callVEP
    vopts = struct;
    vopts.method = "dbnsfp";
    vopts.readall = false;
    vopts.filter = true;
    vopts.dbnsfpHome = opts.dbnsfpHome;
    vopts.verbose = false;
    vopts.ruleN = opts.categories.rule(k);

    for j = 1:numel(filfor.name)
        v = filfor.val(j);
        if filfor.filter(j) == 2
            vopts.(filfor.name(j)) = v;
        elseif filfor.filter(j) == 1
            vopts.(filfor.name(j)) = double(v);
        else
            vopts.(filfor.name(j)) = true;
        end

    end
    vopts = namedargs2cell(vopts);
    vep_filtered = callVEP(bim, vopts{:});

    % vep_filtered = callVEP(bim, 'method', 'dbnsfp', 'readall', false, ...
    %     'filter', true, 'cadd', opts.filterset.cadd(k), ...
    %     'revel', opts.filterset.revel(k), ...
    %     'clinvar', opts.filterset.clinvar(k), ...
    %     'sift', opts.filterset.sift(k), ...
    %     'polyphen', opts.filterset.polyphen(k), ...
    %     'rule', opts.filterset.rule(k),...
    %     'dbnsfpHome', opts.dbnsfpHome, 'verbose', false);
    filter_idx = ismember(vep.id, vep_filtered.id);
    thiscategory = strings(height(vep), 1);
    thiscategory(filter_idx) = opts.filterset.tag(k);
    if k > 1
        vep.categories = join([vep.categories, thiscategory], ",");
    else
        vep.categories = thiscategory;
    end        
end

% %@26OCT2024: test for including autoencoder detected variants
% disp("this is a test for autoencoder variants: autoencoder_full_knn_ids.mat")
% ac_ids = load("D:\MATLAB.dependencies\MAGE\REGENIE\step2.GENE\annotation\dbNSFP\autoencoder_full_knn_ids.mat").autoencoder_full_knn_ids;
% idx = ismember(vep.id, ac_ids);
% thiscategory = strings(height(vep), 1);
% thiscategory(idx) = "Autoencoder";
% vep.categories = join([vep.categories, thiscategory], ",");

vep.categories = regexprep(vep.categories, '(\,+)', ',');
vep.categories = regexprep(vep.categories, '(^,$|,$|^,)', '');

end % END

%% ------------------------------------------------------------------------
function filter_id = parseWESqc(qcpath, opts)

% filter_id = load("filter_id.mat").filter_id;
% 
% chr = natsort(union(filter_id.vmiss_chr, filter_id.hardy_chr));
% for k = 1:numel(chr)
%     idx1 = filter_id.vmiss_chr == chr(k);
%     idx2 = filter_id.hardy_chr == chr(k);
% 
%     filter_id.filter{k, 1} = union(filter_id.vmiss{idx1}, filter_id.hardy{idx2});
% end
% filter_id.chr = chr;
% filter_id = rmfield(filter_id, ["vmiss", "vmiss_chr", "hardy", "hardy_chr"]);
% 
% return

files = getfilenames(qcpath, ["smiss", "vmiss", "hardy"], fullpath=true);
sfiles = files.smiss;
vfiles = files.vmiss;
hfiles = files.hardy;
filter_id = struct;

% filter samples 
if ~isempty(sfiles)
    if isfile("exclude_wes_eid.txt")
        disp("exclude_wes_eid.txt already exists, skipping sample QC step")
    
    else

        tt = tic;
        fprintf("reading %d smiss files for sample qc...", numel(sfiles))
        filter_eid = ({});
        for k = 1:numel(sfiles)
            df = bfilereader(sfiles(k), filterCol="F_MISS", ...
                filter=opts.missingness/100, ...
                operator=">=", extractCol="IID");

            % df = readtable(sfiles(k), "TextType", "string", "FileType", "text", "VariableNamingRule", "preserve");
            % df = df(df.F_MISS >= opts.missingness/100, "IID");

            if ~isempty(df)
                filter_eid{k} = df.IID;
                clear df
            end
        end
    
        tt = toc(tt);
        tt = duration(0, 0, tt);
        fprintf('\b\b done in %s (hr:min:sec)\n', tt)
    
        filter_eid(cellfun(@isempty, filter_eid)) = [];
        filter_eid = unique(vertcat(filter_eid{:}));
        if ~isempty(filter_eid)
            fprintf("found %d samples with a missingness above %.2f %%\n", ...
                numel(filter_eid), opts.missingness/100)
            disp("samples to be excluded are written to exclude_wes_eid.txt file")
            writematrix(filter_eid, "exclude_wes_eid.txt")
        end
    end

else
    disp("no smiss file for sample QC was found!")

end

% filter variant by HWE midp
if isempty(hfiles)
    error("no hardy files (HWE P) for variant qc was found!")
else
    tt = tic;
    fprintf("reading %d hardy files for variant qc\n", numel(hfiles))
    for k = 1:numel(hfiles)
        df = bfilereader(hfiles(k), filterCol="MIDP", ...
            filter=opts.hwp, ...
            operator="<", extractCol="ID");
        if ~isempty(df)
            ids = df.ID;

            % check for monomorphic variants
            df = bfilereader(hfiles(k), ...
                header=true, extractCol=["ID", "O(HET_A1)", "TWO_AX_CT"], ...
                parallel=false);
            df.("O(HET_A1)") = double(df.("O(HET_A1)"));
            idx1 = ismissing(df.("O(HET_A1)")) | df.("O(HET_A1)") == 0;
            idx2 = ismissing(df.TWO_AX_CT) | df.TWO_AX_CT == 0;
            df = df(idx1 & idx2, :);

            if ~isempty(df), ids = union(ids, df.ID); end

            [~, ff, ext] = fileparts(hfiles(k));
            fprintf("\t%s: %d IDs to be excluded.\n", ff+ext, numel(ids))
            filter_id.hardy{k, 1} = ids;

            % add chrom
            df = bfilereader(hfiles(k), ...
                "extractCol", "#CHROM", "extractRow", 1, ...
                "header", true);

            filter_id.hardy_chr(k, 1) = string(df.("#CHROM"));
            clear df
        end
    end
    
    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('done in %s (hr:min:sec)\n', tt)
end

% filter variants
if isempty(vfiles)
    error("no vmiss file for variant QC was found!")
else

    tt = tic;
    fprintf("reading %d vmiss files for variant qc\n", numel(vfiles))
    for k = 1:numel(vfiles)
        df = bfilereader(vfiles(k), filterCol="F_MISS", ...
            filter=opts.missingness/100, ...
            operator=">=", extractCol="ID");

        if ~isempty(df)
            [~, ff, ext] = fileparts(vfiles(k));
            fprintf("\t%s: %d IDs to be excluded.\n", ff+ext, height(df))
            filter_id.vmiss{k, 1} = df.ID;

            % add chrom
            df = bfilereader(vfiles(k), ...
                "extractCol", "#CHROM", "extractRow", 1, ...
                "header", true);

            filter_id.vmiss_chr(k, 1) = string(df.("#CHROM"));
            clear df
        end
    end

    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('done in %s (hr:min:sec)\n', tt)
end


% merge hardy and vmiss IDs
chr = natsort(union(filter_id.vmiss_chr, filter_id.hardy_chr));
for k = 1:numel(chr)
    idx1 = filter_id.vmiss_chr == chr(k);
    idx2 = filter_id.hardy_chr == chr(k);

    filter_id.filter{k, 1} = union(filter_id.vmiss{idx1}, filter_id.hardy{idx2});
end
filter_id.chr = chr;
filter_id = rmfield(filter_id, ["vmiss", "vmiss_chr", "hardy", "hardy_chr"]);

end

