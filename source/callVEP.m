function [res, r] = callVEP(bim, opts)
% callVEP calls EnsemblREST to annotate input variants using online VEP.
% inputs must be a struct with fields: snp, chr, pos, a1 and a2, which a2
% is minor allele. returns a table "res" containing all information for 
% each input variant.
% 
% Oveis Jamialahmadi, Sahlgrenska Akademy, May 2021.
% 
% @08/17/2021 'clinvar' flag was added.
% @08/17/2021 'method' option was added. Two options are
%             available: "vep", which annotates using REST Ensembl API,
%             while "dbnsfp" uses local files parsed with dbNSFPparser
%             function. Note that with "dbnsfp" option, some filtering
%             methods may not be available, and depends on available
%             predictions available within the parsed files. Current
%             available columns for dbNSFP4.2a are rs_dbSNP, SIFT4G_pred,
%             Polyphen2_HVAR_pred, REVEL_score, CADD_phred and
%             clinvar_clnsig.
% @25/08/2021  bug fixed in parsing the 'spliceai' annotation.
% @18/10/2021  LOFTEE column was added to the output 'res' table.
% @30/10/2021  'loftee' flag was added (default: false).
% @23/12/2021  a bug was fixed for single query.
% @19/04/2022  maximum POST input size now has been changed to 200
%              (before it was 100).
% @01/11/2022 now returns domains and domains_db in the parsed table.
%             However, filtering yet should be implemented.
% @02/12/2022 some bugs with writevep were fixed.
% 
% @12JULY2024: now supports further methods (with dbNSFP only): lrt,
% mutationtaster, m_cap and alphamissense
% 
% @12JULY2024: 'ruleN' option was added to define more complex merging
% rules among filters. This flag overrides 'rule' and for now only works
% with dbNSFP filtering. 'ruleN' is a numeric value corresponding to the
% number of consensus among queried methods. e.g. if 'ruleN' is 1, it's
% equivalent to 'rule' = "or". 
% 
% TODO: constructing variant ID from 'bim' requires the native (unadjusted)
% a1 and a2 order (as in original BIG/BIM). This currently is not checked,
% especially that getbulkgeno function reformats 'bim' so that a2 is always
% minor. 03/10/2022 

arguments
    bim {mustBeA(bim, 'struct')}
    opts.method {mustBeMember(opts.method, ["vep", "dbnsfp"])} = "dbnsfp" % "vep" annotates using REST Ensembl API, while "dbnsfp" uses local files parsed with dbNSFPparser function.
    opts.genomeRef {mustBeMember(opts.genomeRef, ["37", "38"])} = "38" % genome reference
    opts.verbose (1,1) logical = true
    opts.filter (1,1) logical = false % to fiter non-LoF variants for additional predictors
    opts.loftee (1,1) logical = false % keeps only HC variants.
    opts.cadd (1,1) double = 0; % filter variants with CADD-phred above this
    opts.revel (1,1) double = 0; % filter variants with REVEL above this
    opts.sift {mustBeMember(opts.sift, ["strict", "relaxed", "none"])} = "none";
    opts.polyphen {mustBeMember(opts.polyphen, ["strict", "relaxed", "none"])} = "none"; % probably_damaging is "D" and possibly_damaging "P" in dbNSFP 
    opts.primateai (1,1) logical = false; % if primateAI == "D"
    opts.spliceai (1,1) double = 0; % filter variants with spliceAI value above this. Hint: 0.2 (high recall) | 0.5 (recommended) | 0.8 (high precision)
    opts.clinvar (1,1) logical = false; % if true, matches ClinVar flags of "Likely_pathogenic" and "pathogenic"
    opts.rule {mustBeMember(opts.rule, ["and", "or"])} = "and"; % how to merge filters
    opts.backup (1,1) logical = true; % save raw annotations, so it can be used for future use (e.g. with different filters). In this case, raw annotations will be saved to a mat file within callVEP directory. 
    
    % dbNSFP options
    opts.dbnsfpHome {mustBeTextScalar} % dbNSFP parsed annotations directory. Files should be available for each chromosome separately. only used when method is set to "dbnsfp"
    opts.dbnsfpFile {mustBeTextScalar} = "" % dbNSFP file pattern, this is helpful if there are non dbNSFP files in 'dbnsfpHome' dir. This can be part of the file name for all chromosomes, e.g. "dbNSFP4.2a_variant.chr" matches dbNSFP4.2a_variant.chr1.txt to dbNSFP4.2a_variant.chr22.txt. If left empty (default), only dbNSFP parsed files (txt) are present in the dir.
    opts.parallel (1,1) logical = false % Only used with "dbnsfp". If true uses bfilereader in parallel, otherwise (default) uses readtable (use parallel only if parsed files are big enough).
    opts.readall (1,1) logical = true % by default (true) fetches all columns from the dbNSFP parsed file (useful in gwasrunner). If false, then only filtering columns will be read (useful in saigeWrapper).

    % new optsions @12JULY2024
    opts.lrt (1,1) logical = false
    opts.mutationtaster (1,1) logical = false
    opts.m_cap (1,1) logical = false
    opts.alphamissense (1,1) logical = false;
    opts.ruleN (1,1) double
end

checkFields = fieldnames(bim);

switch opts.method
    case "vep"
        if any(~ismember({'snp', 'chr', 'pos', 'a1', 'a2'}, checkFields))
            error('callVEP:inputs', 'input struct must contain snp, chr, pos, a1 and a2 fields!')
        end
        [res, r] = getVEP(bim, bim.chr + " " + bim.pos + " . " + bim.a2 + " " + bim.a1 + " . . .", opts);
    case "dbnsfp"
        % if bim.snp is empty, then all variants will be filtered.
        if any(~ismember({'snp', 'chr'}, checkFields))
            error('callVEP:inputs', 'input struct must contain snp and chr fields!')
        end
        r = [];
        res = getdbNSFP(unique(bim.chr), bim.snp, opts);
    otherwise
        error('callVEP:unknownmethod', 'only vep or dbnsfp is allowed')
end

if opts.verbose
    fprintf('\b\b Done.\n')
    toc
end

end % END

%% subfunctions ===========================================================
function r = getrawVEP(vcf, opts)
chunksize = 200;
ct = [1 chunksize];
flag = true;
r = [];
while flag
    if ct(2) > numel(vcf)
        ct(2) = numel(vcf);
        flag = false;
    end
    if opts.verbose
        fprintf('reading ids in range: [%d - %d]...', ct(1), ct(2))
    end
    r = [r; EnsemblREST(vcf(ct(1):ct(2)), 'vep', 'refGenome', opts.genomeRef)];
    ct = ct + chunksize;
    if opts.verbose
        fprintf('\b\b Done.\n')
    end
    if ct(1) >= numel(vcf)
        flag = false;
    end
end
end

%% ------------------------------------------------------------------------
function idx = createEmptyVec(n, rule)
if rule == "or"
    idx = false(n, 1);
elseif rule == "and"
    idx = true(n, 1);
end
end

%% ------------------------------------------------------------------------
function writeBIM(bim, vcf, remvcf, wd)
if ~isempty(remvcf) % only write vcf/vcfid for new unannotated variants
    rmidx = ~ismember(vcf, remvcf);
    flds = fieldnames(bim);
    for i = 1:numel(flds)
        bim.(flds{i})(rmidx) = [];
    end
    vcf(rmidx) = [];
end

uchr = unique(bim.chr);
chr = bim.chr;
flds = fieldnames(bim);
bim.vcf = vcf; clear vcf
bim.vcfid = "chr" + bim.chr + "_" + bim.pos + "_" + bim.a2 + "_" + bim.a1;
bim = rmfield(bim, flds);
flds = fieldnames(bim);

for i = 1:numel(uchr)
    idx = chr == uchr(i);
    bimtmp = struct;
    for j = 1:numel(flds)
        bimtmp.(flds{j}) = bim.(flds{j})(idx);
    end
    fileName = fullfile(wd, "callVEP.chr" + uchr(i) + ".mat");
    
    if exist(fileName, 'file')
        bimfile = load(fileName, 'vcf', 'vcfid');
        rmidx = ismember(bimtmp.vcf, bimfile.vcf);
        for j = 1:numel(flds) % remove overlapped vcf ids
            bimtmp.(flds{j})(rmidx) = [];
            wbim = bimtmp.(flds{j});
            if isempty(wbim); wbim = []; end
            bimtmp.(flds{j}) = [wbim; bimfile.(flds{j})];
        end
        save(fileName, '-struct', 'bimtmp', '-append')
    else
        save(fileName, '-struct', 'bimtmp')
    end
    
end
end

%% ------------------------------------------------------------------------
function writeann(rALL, bim, vcf, wd)
bim.vcf = vcf; clear vcf
bim.vcfid = "chr" + bim.chr + "_" + bim.pos + "_" + bim.a2 + "_" + bim.a1;
uchr = unique(bim.chr);

if isstruct(rALL), rALL = {rALL}; end % only 1 varaint 
rvcf = string(cellfun(@(x)x.input, rALL, 'uni', false));
rchr = string(cellfun(@(x)x.seq_region_name, rALL, 'uni', false));

for i = 1:numel(uchr)
    
    % load 'vcf' and 'vcfid' 
    rname = fullfile(wd, "callVEP.chr"+uchr(i)+".mat");
    bimfile = load(rname, 'vcf', 'vcfid');
    
    % get only uchr(i) subset of r 
    ridx = rchr == string(uchr(i));
    rvcfchr = rvcf(ridx);
    rALLchr = rALL(ridx);
    rvcf(ridx) = [];
    rchr(ridx) = [];
    rALL(ridx) = []; 
    
    [ridx, ridx2] = ismember(rvcfchr, bimfile.vcf);
    rALLchr(~ridx) = []; rALLchrID = bimfile.vcfid(ridx2);
    if ~isempty(rALLchrID)
        rALLchr = cell2struct(rALLchr, rALLchrID);
        save(rname, '-struct', 'rALLchr', '-append')
    end
    
    if isempty(rALL) % nothing left to process
        return
    end
end
end

%% ------------------------------------------------------------------------
function [r, remvcf] = readBackup(bim, vcf, annfiles)
% read annotations structs from backup files (matfiles) and further checks
% if there are unannotated vcf in input query, and if so, returns them in
% remvcf string.

wd = fileparts(annfiles{1}); % target directory: callVEP.backup
bim.vcf = vcf; clear vcf
bim.vcfid = "chr" + bim.chr + "_" + bim.pos + "_" + bim.a2 + "_" + bim.a1;
uchr = unique(bim.chr);
remvcf = []; r = [];
for i = 1:numel(uchr)
    idx = bim.chr == uchr(i);
    rname = fullfile(wd, "callVEP.chr" + uchr(i) + ".mat");
    
    fileIdx = ismember(annfiles, rname);
    if any(fileIdx)
        bimfile = load(annfiles{fileIdx}, 'vcf', 'vcfid');

        sdiff = setdiff(bim.vcf(idx), bimfile.vcf);
        if isempty(sdiff)
            sdiff = [];
            query_vcfids = [];
        else
            % for missing raw data but existing vcfid/vcf in cache
            query_vcfids = bim.vcfid(ismember(bim.vcf, sdiff));
        end
        remvcf = [remvcf; sdiff]; % unannotated variants (should be fetched)
        
        fetchidx = bimfile.vcfid(ismember(bimfile.vcf, bim.vcf(idx)));
        if ~isempty(fetchidx) % check if has vcfid/vcf but missing raw data
            warning('off')
            rtmp = load(rname, fetchidx{:});
            check_vcfids = bim.vcfid(idx);

            % exclude already missing/unannotated query variants for this chr
            if ~isempty(query_vcfids), check_vcfids = setdiff(check_vcfids, query_vcfids); end
            sdiff = setdiff(check_vcfids, fieldnames(rtmp));
            if ~isempty(sdiff)
                sdiff = bim.vcf(ismember(bim.vcfid, sdiff));
                remvcf = [remvcf; sdiff];

                % remove vcfid/vcf from cache mat file (for which there is
                % no raw vep data)
                idx = ismember(bimfile.vcf, sdiff);
                bimfile.vcf(idx) = []; bimfile.vcfid(idx) = [];
                save(rname, '-struct', 'bimfile', '-append')
            end
            rtmp = struct2cell(rtmp);
            warning('on')
            if isempty(rtmp) % Rest API broke: so bim exists without any annotation
                remvcf = [remvcf; bim.vcf(idx)];
            else
                r = [r; rtmp];
            end
        end
    else % unannotated
        remvcf = unique([remvcf; bim.vcf(idx)]);
    end
end

remvcf = unique(remvcf, "stable");
end

%% ------------------------------------------------------------------------
function res = filterPred(res, opts)
    cols = string(res.Properties.VariableNames);

    %@12JULY2024: for dbNSFP only
    if opts.method == "dbnsfp"
        cols = intersect(cols, ["CADD-phred", "REVEL", "SIFT", "Polyphen", ...
            "ClinVar", "LRT", "PrimateAI", "MutationTaster", "M_CAP", ...
            "AlphaMissense"]);
        idx = false(height(res), numel(cols));
        for k = 1:numel(cols)
            if (cols(k) == "CADD-phred") && opts.cadd > 0
                tmpidx = res.("CADD-phred") >= opts.cadd;
            elseif (cols(k) == "REVEL") && opts.revel > 0
                tmpidx = res.("REVEL") >= opts.revel;
            elseif cols(k) == "SIFT"
                tmpidx = res.SIFT.contains("D");
            elseif cols(k) == "Polyphen"
                if opts.polyphen == "relaxed" % either D or P
                    tmpidx = res.Polyphen.contains(["D", "P"]);
                else
                    tmpidx = res.Polyphen.contains("D");
                end
            elseif cols(k) == "ClinVar"
                tmpidx = res.ClinVar.lower.contains("pathogenic") & ~res.ClinVar.lower.contains("conflicting");
            elseif any(cols(k) == ["PrimateAI", "LRT", "M_CAP"])
                tmpidx = res.(cols(k)).contains("D");
            elseif cols(k) == "MutationTaster"
                tmpidx = res.MutationTaster.contains(["A", "D"]);
            elseif cols(k) == "AlphaMissense"
                tmpidx = res.AlphaMissense.contains("P");
            else
                error("callVEP:unrecognized filtering algorithm!")
            end
            idx(:, k) = tmpidx;

        end

        idx = sum(idx, 2);
        if isfield(opts, "ruleN")
            idx = idx >= opts.ruleN;
        elseif opts.rule == "or"
            idx = idx >= 1;
        else % and
            idx = idx == numel(cols);
        end

        res(~idx, :) = [];

        return
    end

    if opts.cadd > 0 && any(cols == "CADD-phred")
        idx_cadd = res.("CADD-phred") >= opts.cadd;
    else
        idx_cadd = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.revel > 0 && any(cols == "REVEL")
        idx_revel = res.REVEL >= opts.revel;
    else
        idx_revel = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.spliceai > 0 && any(cols == "spliceai")
        idx_spliceai = res.spliceai >= opts.spliceai;
    else
        idx_spliceai = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.sift ~= "none" && any(cols == "SIFT")
        if strcmp(opts.method, 'dbnsfp') % for dbNSFP, D is damaging, T is tolerated
            idx_sift = contains(res.SIFT, "D"); % relaxed and strict are treated the same in this case
        else
            if opts.sift == "relaxed" % also low confidence damaging
                idx_sift = contains(res.SIFT, "deleterious");
            else % damaging but not low_confidence (though low_confidence can be present for other transcripts)
                idx_sift = ~cellfun(@isempty,regexp(res.SIFT, 'deleterious(?!_)'));
            end
        end
    else % skip SIFT
        idx_sift = createEmptyVec(size(res, 1), opts.rule);
    end

    if opts.loftee && any(cols == "LOFTEE")
        idx_loftee = contains(res.LOFTEE, "HC");
    else % skip SIFT
        idx_loftee = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.polyphen ~= "none" && any(cols == "Polyphen")
        if strcmp(opts.method, 'dbnsfp') % for dbNSFP, D is probably damaging, P is possibly damaging, B is benign
            if opts.polyphen == "relaxed" % either D or P
                idx_polyphen = contains(res.Polyphen, ["D", "P"]);
            else % only D (not anymore, see below)
                % @12JULY2024: changed to presence of "D" at least once
                % among trnscripts (see UKBB WES 500 K paper). This is more
                % relaxed than picking only "D"
                idx_polyphen = contains(res.Polyphen, "D"); % & ~contains(res.Polyphen, ["B", "P"]);
            end
        else
            if opts.polyphen == "relaxed" % anything with damaging
                idx_polyphen = contains(res.Polyphen, "damaging");
            else % only probably_damaging
                idx_polyphen = contains(res.Polyphen, "probably_damaging");
            end
        end
    else % skip PolyPhen
        idx_polyphen = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.primateai && any(cols == "primateai_pred")
        idx_primateai = res.primateai_pred == "D";
    else
        idx_primateai = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.clinvar && any(cols == "ClinVar")
        idx_clinvar = contains(lower(res.ClinVar), "pathogenic") & ~contains(lower(res.ClinVar), "conflicting");
    else
        idx_clinvar = createEmptyVec(size(res, 1), opts.rule);
    end
    
    if opts.rule == "or"  
        idx =  idx_cadd | idx_revel | idx_sift | idx_polyphen | idx_primateai | idx_spliceai | idx_clinvar | idx_loftee;
    elseif opts.rule == "and"
        idx =  idx_cadd & idx_revel & idx_sift & idx_polyphen & idx_primateai & idx_spliceai & idx_clinvar & idx_loftee;
    end
    
    res(~idx, :) = [];
end

%% ------------------------------------------------------------------------
function [res, r] = getVEP(bim, vcf, opts)

if opts.backup % save/use/load backup files for variants in this vcf
    wd = fullfile(fileparts(which('callVEP.m')), 'callVEP.backup', "GRCh" + opts.genomeRef);
    if ~exist(wd, 'dir')
        mkdir(wd)
    end
    matfiles = {dir(fullfile(wd, '*.mat')).name}.';
    matfiles(~startsWith(matfiles, 'callVEP')) = [];
    if ~isempty(matfiles)
        matfiles = fullfile(wd, matfiles);
        opts.backupann = false;
        [r, remvcf] = readBackup(bim, vcf, matfiles);
        if ~isempty(remvcf)
            opts.backupbim = true;
        else
            opts.backupbim = false;
        end
    else % else: no backup file exists in callVEP.backup directory (first time use)
        opts.backupbim = true;
        remvcf = []; opts.backupann = true; % no backedup annotation files exists
    end
    
    if opts.backupbim % write vcf/bim to a mat file
        writeBIM(bim, vcf, remvcf, wd)
    end
end

if opts.verbose
    fprintf('getting VEP raw for %d ids...\n', numel(vcf))
    tic
end

if opts.backup
    if ~isempty(remvcf) % there are unannotated vcf that cannot be found in backup files
        r_rem = getrawVEP(remvcf, opts);
        if isstruct(r_rem), r_rem = {r_rem}; end
        r = [r; r_rem];
        writeann(r_rem, bim, vcf, wd)
    elseif opts.backupann % annotation files don't exist (first time use)
        r = getrawVEP(vcf, opts);
        writeann(r, bim, vcf, wd)
    end
else
    r = getrawVEP(vcf, opts);
end

if ~iscell(r); r = {r}; end % single query

if opts.verbose 
    fprintf('processing VEP raw info...')
end

res = table(strings(numel(r), 1), strings(numel(r), 1), ...
    strings(numel(r), 1), strings(numel(r), 1), nan(numel(r), 1),...
    strings(numel(r), 1), strings(numel(r), 1), nan(numel(r), 1),...
    strings(numel(r), 1), strings(numel(r), 1),strings(numel(r), 1),...
    strings(numel(r), 1), strings(numel(r), 1), nan(numel(r), 1),...
    strings(numel(r), 1), strings(numel(r), 1), strings(numel(r), 1),...
    strings(numel(r), 1), strings(numel(r), 1), strings(numel(r), 1), ...
    'VariableNames', {'id', 'snp', 'consequence', 'gene',...
    'CADD-phred', 'SIFT', 'Polyphen', 'REVEL', 'ClinVar', 'clin_sig',...
    'hgvsc', 'hgvsp', 'primateai_pred', 'spliceai', 'spliceai_pred',...
    'LOFTEE', 'domains', 'domains_db', 'interpro_domain', 'AlphaMissense'});

for i = 1:numel(r)
    res.id(i) = string(r{i}.input);
    clinsig = [];
    if isfield(r{i}, 'colocated_variants')
        if numel(r{i}.colocated_variants) > 1
            try
                res.snp(i) = join(string(cellfun(@(x)x.id, r{i}.colocated_variants, 'uni', false)), ';');
            catch % struct
                res.snp(i) = join(string({r{i}.colocated_variants.id}), ';');
            end
        else
            res.snp(i) = string(r{i}.colocated_variants.id);
        end
        
        tmp_coloc = r{i}.colocated_variants;
        for j = 1:numel(r{i}.colocated_variants)
            
            if iscell(tmp_coloc)
                if isfield(r{i}.colocated_variants{j}, 'clin_sig_allele')
                    clinsig = [clinsig; string(r{i}.colocated_variants{j}.clin_sig_allele)];
                end
            elseif isfield(r{i}.colocated_variants, 'clin_sig_allele')
                clinsig = [clinsig; string(r{i}.colocated_variants.clin_sig_allele)];
            end
        end
    end
    
    res.consequence(i) = string(r{i}.most_severe_consequence);
    
    [cadd, sift, polyphen, revel, clinvar, hgvsc, hgvsp, primateai,...
        spliceai, spliceai_pred, gene, domains, domains_db, interpro_domain,...
        loftee, alphamissense] = deal([]);
    
    try
        if isstruct(r{i}.transcript_consequences) % only 1 consequence
            tmp = r{i}.transcript_consequences;
            r{i} = rmfield(r{i}, 'transcript_consequences');
            r{i}.transcript_consequences{1} = tmp;
        end
    catch % no transcript_consequences exists
        r{i}.transcript_consequences{1} = '';
    end
    
    for j = 1:numel(r{i}.transcript_consequences)
        if isfield(r{i}.transcript_consequences{j}, 'consequence_terms') &&...
                isfield(r{i}.transcript_consequences{j}, 'gene_symbol')
            try
                if any(string(r{i}.transcript_consequences{j}.consequence_terms) == res.consequence(i))
                    gene = [gene; string(r{i}.transcript_consequences{j}.gene_symbol)];
                end
            catch

                tmpcons = {r{i}.transcript_consequences{j}(1).consequence_terms};
                if any(string(vertcat(tmpcons{:})) == res.consequence(i))
                    tmpGene = {r{i}.transcript_consequences{j}(1).gene_symbol};
                    tmpGene = string(vertcat(tmpGene{:}));
                    gene = [gene; join(unique(tmpGene), ';')];
                end

            end
        end
        
        if isfield(r{i}.transcript_consequences{j}, 'lof')
            loftee = [loftee; string(r{i}.transcript_consequences{j}.lof)];
        end
        if isfield(r{i}.transcript_consequences{j}, 'cadd_phred')
            cadd = [cadd; r{i}.transcript_consequences{j}.cadd_phred];
        end
        if isfield(r{i}.transcript_consequences{j}, 'sift_prediction')
            sift = [sift; string(r{i}.transcript_consequences{j}.sift_prediction)];
        end
        if isfield(r{i}.transcript_consequences{j}, 'polyphen_prediction')
            polyphen = [polyphen; string(r{i}.transcript_consequences{j}.polyphen_prediction)];
        end
        if isfield(r{i}.transcript_consequences{j}, 'revel_score')
            revel = [revel; string(r{i}.transcript_consequences{j}.revel_score)];
        end
        if isfield(r{i}.transcript_consequences{j}, 'clinvar_clnsig')
            clinvar = [string(clinvar); r{i}.transcript_consequences{j}.clinvar_clnsig];
        end

        if isfield(r{i}.transcript_consequences{j}, 'alphamissense')
            alphamissense = [alphamissense; string(r{i}.transcript_consequences{j}.alphamissense.am_class)];
        end

        if isfield(r{i}.transcript_consequences{j}, 'domains')
            try
                tmpDomains = cellfun(@(x)x.name, r{i}.transcript_consequences{j}.domains, 'uni', false);
                tmpDomains_db = cellfun(@(x)x.db, r{i}.transcript_consequences{j}.domains, 'uni', false);
            catch
                try
                    tmpDomains = string({r{i}.transcript_consequences{j}.domains.name}.');
                    tmpDomains_db = string({r{i}.transcript_consequences{j}.domains.db}.');
                catch
                    tmpDomains = ""; tmpDomains_db = "";
                end
            end
            domains = [string(domains); string(tmpDomains)];
            domains_db = [string(domains_db); string(tmpDomains_db)];
        end

        if isfield(r{i}.transcript_consequences{j}, 'interpro_domain')
            interpro_domain = [string(interpro_domain); string(unique({r{i}.transcript_consequences{j}.interpro_domain}))];
        end
        if isfield(r{i}.transcript_consequences{j}, 'hgvsc')
            hgvsc = [string(hgvsc); r{i}.transcript_consequences{j}.hgvsc];
        end
        if isfield(r{i}.transcript_consequences{j}, 'hgvsp')
            hgvsp = [string(hgvsp); r{i}.transcript_consequences{j}.hgvsp];
        end
        if isfield(r{i}.transcript_consequences{j}, 'primateai_pred')
            primateai = [string(primateai); r{i}.transcript_consequences{j}.primateai_pred];
        end
        if isfield(r{i}.transcript_consequences{j}, 'spliceai')
            if numel(r{i}.transcript_consequences{j}) > 1
                r{i}.transcript_consequences{j} = r{i}.transcript_consequences{j}(1);
            end
            if strcmp(r{i}.transcript_consequences{j}.consequence_terms, res.consequence(i))
                spliceaitmp = r{i}.transcript_consequences{j}.spliceai;
                % info: https://raw.githubusercontent.com/ensembl-variation/VEP_plugins/master/SpliceAI.pm
                spliceaitmp = [spliceaitmp.DS_AG;spliceaitmp.DS_AL;spliceaitmp.DS_DG;spliceaitmp.DS_DL];
                spliceaitmp_letter = ["AG", "AL", "DG", "DL"];
                [spliceaitmp, spliceaitmp_idx] = max(spliceaitmp);
                spliceai = [spliceai; spliceaitmp];
                spliceai_pred = [string(spliceai_pred); "(" + spliceaitmp_letter(spliceaitmp_idx) + ")"];
            end
        end
        
    end
    
    if ~isempty(loftee)
        res.("LOFTEE")(i) = join(unique(loftee, 'stable'), ';');
    end
    if ~isempty(cadd)
        res.("CADD-phred")(i) = unique(cadd);
    end
    if ~isempty(polyphen)
        res.Polyphen(i) = join(unique(polyphen, 'stable'), ';');
    end
    if ~isempty(sift)
        res.SIFT(i) = join(unique(sift, 'stable'), ';');
    end
    if ~isempty(revel)
        try
            revel = extract(string(revel),  regexpPattern('\d+\.?\d*') | digitsPattern);
            res.REVEL(i) = max(double(unique(revel(:))));
        catch
            res.REVEL(i) = nan;
        end
    end
    if ~isempty(clinvar)
        res.ClinVar(i) = unique(clinvar);
    end
    if ~isempty(clinsig)
        res.clin_sig(i) = unique(clinsig);
    end

    if ~isempty(alphamissense)
        res.AlphaMissense(i) = join(unique(alphamissense, 'stable'), ';');
    end
    
    if ~isempty(hgvsc)
        res.hgvsc(i) = join(unique(hgvsc, 'stable'), ';');
    end
    if ~isempty(hgvsp)
        res.hgvsp(i) = join(unique(hgvsp, 'stable'), ';');
    end
    if ~isempty(primateai)
        res.primateai_pred(i) = join(unique(primateai, 'stable'), ';');
    end
    if ~isempty(spliceai)
        try
            spliceai = unique(spliceai, 'stable');
            [~, spliceaiIdx] = max(spliceai);
            spliceai = spliceai(spliceaiIdx);
            res.spliceai(i) = spliceai;
            res.spliceai_pred(i) = join(unique(spliceai_pred, 'stable'), ';');
        catch
            % leave it empty
        end
    end
    if ~isempty(domains)
        res.domains(i) = join(unique(domains, "stable"), ';');
        res.domains_db(i) = join(unique(domains_db, "stable"), ';');
    end
    if ~isempty(interpro_domain)
        res.interpro_domain(i) = join(unique(interpro_domain), ';');
    end
    if ~isempty(gene)
        res.gene(i) = join(unique(gene, 'stable'), ';');
    end
end

if all(res.domains == "")
    res.domains = []; res.domains_db = [];
end
if all(res.interpro_domain == "")
    res.interpro_domain = [];
end

[~, f] = ismember(vcf, res.id);
res = res(f, :);
res.id = bim.snp;

lof_filters = ["stop_gained"; "splice_donor_variant"; ...
            "frameshift_variant"; "splice_acceptor_variant"];
lof_missense_filters = ["stop_gained"; "start_lost"; "splice_donor_variant"; ...
    "frameshift_variant"; "splice_acceptor_variant"; "stop_lost";...
    "transcript_ablation"; "missense_variant"; ...
    "protein_altering_variant"; "inframe_insertion"; "inframe_deletion"];

[res.LoF, res.LoFMissense] = deal(false(height(res), 1));
res.LoF(contains(res.consequence, lof_filters)) = true;
res.LoFMissense(contains(res.consequence, lof_missense_filters)) = true;

if opts.filter % additional filters for non-LoF variants
    resKeep = res(~res.LoF, :); % non-LoF variants
    res(~res.LoF, :) = [];
    if opts.loftee && any(ismember(res.Properties.VariableNames, 'LOFTEE'))
        res(~contains(res.LOFTEE, "HC"), :) = [];
    end
    
    res = [res; filterPred(resKeep, opts)];
end
end

%% ------------------------------------------------------------------------
function res = getdbNSFP(chr, ids, opts)
if ~exist(opts.dbnsfpHome, 'dir')
    error('callVEP:getdbNSFP', '%s is not a valid directory!', opts.dbnsfpHome)
end

dbfiles = getfilenames(opts.dbnsfpHome, "txt").txt;
if opts.dbnsfpFile ~= ""
    dbfiles(~contains(dbfiles, opts.dbnsfpFile)) = [];
end

if isempty(dbfiles)
    error('callVEP:getdbNSFP', 'no text file was found in %s!', opts.dbnsfpHome)
end

if numel(dbfiles) > 1 % multiple files are available
    dbfiles = natsort(dbfiles);
    try
        % chrPos = checkfile{1} ~= checkfile{2};
        % chrPos = @(x) sprintf([checkfile{1}(1:find(chrPos)-1),'%s', ...
        %     checkfile{1}(find(chrPos)+1:numel(chrPos))], x);
        chrPos = strshare(dbfiles, 'pattern', true);
        chrPos = @(x) sprintf(strrep(chrPos, '%d', '%s'), x);
    catch
        error('callVEP:getdbNSFP', 'cannot find a file pattern, try setting dbnsfpFile option')
    end
    
    chr = string(chr);
    tmpdb = string;
    for i = 1:numel(chr)
        tmpdb(i) = chrPos(chr(i));
        if ~any(dbfiles == tmpdb(i)) && chr(i) == "23"
            tmpdb(i) = chrPos("X");
        elseif ~any(dbfiles == tmpdb(i)) && chr(i) == "24"
            tmpdb(i) = chrPos("Y");
        end
    end
    dbfiles = intersect(tmpdb, dbfiles);
    if isempty(dbfiles)
        error('callVEP:getdbNSFP', 'cannot find dbNSFP file(s) for query chr(s) within %s', opts.dbnsfpHome)
    end
end
dbfiles = fullfile(opts.dbnsfpHome, dbfiles);

% check columns to be fetched
extCols = ["id", "rs_dbSNP", "genename", "Ensembl_geneid"];
if opts.readall
    extCols = [extCols, "CADD_phred", "REVEL_score", "SIFT4G_pred", ...
        "Polyphen2_HVAR_pred", "clinvar_clnsig", "LRT_pred", ...
        "MutationTaster_pred", "M-CAP_pred", "PrimateAI_pred", ...
        "AlphaMissense_pred"];
else
    if opts.cadd > 0
        extCols = [extCols, "CADD_phred"];
    end
    if opts.revel > 0
        extCols = [extCols, "REVEL_score"];
    end
    if opts.sift ~= "none"
        extCols = [extCols, "SIFT4G_pred"];
    end
    if opts.polyphen ~= "none"
        extCols = [extCols, "Polyphen2_HVAR_pred"];
    end
    if opts.clinvar
        extCols = [extCols, "clinvar_clnsig"];
    end
    if opts.lrt 
        extCols = [extCols, "LRT_pred"];
    end
    if opts.primateai
        extCols = [extCols, "PrimateAI_pred"];
    end
    if opts.mutationtaster
        extCols = [extCols, "MutationTaster_pred"];
    end
    if opts.m_cap
        extCols = [extCols, "M-CAP_pred"];
    end
    if opts.alphamissense
        extCols = [extCols, "AlphaMissense_pred"];
    end
end

res = cell(numel(dbfiles), 1);
for i = 1:numel(dbfiles) % loop over dbfiles and read annotations
    if opts.parallel
        if isempty(ids) % fetch all variants
            res{i} = bfilereader(dbfiles(i), 'extractCol', extCols, ...
                'header', true, 'return', 'rawTable', 'verbose', 'off',...
                'parallel', false);
        else
            res{i} = bfilereader(dbfiles(i), 'extractCol', extCols, ...
                'pattern', "\b"+ids+"\b", 'patternCol', 'id', 'header', true, 'return', ...
                'rawTable', 'verbose', 'off', 'parallel', false);
        end
    else
        fileOpts = detectImportOptions(dbfiles(i), ...
            'PreserveVariableNames', true, 'TextType', 'string', ...
            'NumHeaderLines', 0);
        fileOpts.SelectedVariableNames = extCols;
        restmp = readtable(dbfiles(i), fileOpts);
        if ~isempty(ids)
            restmp(~ismember(restmp.id, ids), :) = [];
        end
        res{i} = restmp;
        clear restmp
    end
end

res = vertcat(res{:});

if opts.cadd > 0 || opts.readall
    res.CADD_phred = double(res.CADD_phred);
end
if opts.revel > 0 || opts.readall
    res.REVEL_score = double(res.REVEL_score);
end

% make column names consistent with vep module
cols = res.Properties.VariableNames;
cols(strcmp(cols, 'CADD_phred')) = {'CADD-phred'};
cols(strcmp(cols, 'REVEL_score')) = {'REVEL'};
cols(strcmp(cols, 'SIFT4G_pred')) = {'SIFT'};
cols(strcmp(cols, 'Polyphen2_HVAR_pred')) = {'Polyphen'};
cols(strcmp(cols, 'clinvar_clnsig')) = {'ClinVar'};
cols(strcmp(cols, 'LRT_pred')) = {'LRT'};
cols(strcmp(cols, 'PrimateAI_pred')) = {'PrimateAI'};
cols(strcmp(cols, 'MutationTaster_pred')) = {'MutationTaster'};
cols(strcmp(cols, 'M-CAP_pred')) = {'M_CAP'};
cols(strcmp(cols, 'AlphaMissense_pred')) = {'AlphaMissense'};
res.Properties.VariableNames = cols;

if opts.filter
    res = filterPred(res, opts);
end
end