function res = MRlap(efile, ofile, opts)
% a wrapper for MRlap R package: https://github.com/n-mounier/MRlap
% Oveis Jamialahmadi, University of Gothenburg, 08 NOV 2024.

arguments
    efile {mustBeFile} % exposure summary stat file (can be gzipped)
    ofile {mustBeFile} % outcome summary stat file (can be gzipped)

    % MRlap arguments
    % The path to the folder in which the LD scores used in the analysis
    % are located. Expects LD scores formated as required by the original
    % LD score regression software.
    opts.ld {mustBeFolder}

    % The path to a file of SNPs with alt, ref alleles and rsid used to
    % allign alleles across traits
	opts.hm3 {mustBeFile}

    opts.exposure_name {mustBeTextScalar} = "exposure"
    opts.outcome_name {mustBeTextScalar} = "outcome"
	opts.do_pruning (1,1) logical = true
	opts.user_SNPsToKeep {mustBeText, mustBeVector}
	opts.MR_threshold (1,1) double = 5e-08
	opts.MR_pruning_dist (1,1) double = 500
	opts.MR_pruning_LD (1,1) double = 0
	opts.MR_reverse (1,1) double = 0.001
	opts.MR_plink {mustBeFile}
	opts.MR_bfile {mustBeFile}
	opts.save_logfiles (1,1) logical = true

    % developer options (non-native)
    opts.remove_loci {mustBeA(opts.remove_loci, "table")} % loci to be removed from efile or ofile, must be a table with snp/pos/chr columns (case insensitive)
    opts.remove_from {mustBeTextScalar, mustBeMember(opts.remove_from, ["exposure", "outcome"])} = "outcome" % choose the smallest file in size for faster computation
    opts.r2 (1,1) double = 0.01 % R2 for removing in LD pairs.
    opts.bgenhome {mustBeFolder}
end

tmp_name = getRandomName("", 5);

files = [efile; ofile];
for k = 1:numel(files)
    [~, ~, ext] = fileparts(files(k));
    
    if k == 1 % exposure
        name = "exposure";
    else
        name = "outcome";
    end

    name = name + tmp_name; 

    if ~isfile(name + ext)
        copyfile(files(k), name + ext)
    end
    files(k) = name + ext;
end

efile = files(1); ofile = files(2);

% change colnames
ecol_old = bfilereader(efile, summary="firstline");
ocol_old = bfilereader(ofile, summary="firstline");

ecol_new = createGWASheader(ecol_old);
ocol_new = createGWASheader(ocol_old);
ecol_new(ecol_new == "BP") = "pos";
ocol_new(ocol_new == "BP") = "pos";

% check if modifying column names is necessary
cc = ["snp", "a1", "a2", "n", "beta", "se", "chr", "pos"];
if any(~ismember(cc, ecol_old.lower))
    changeColnames(efile, ecol_old, ecol_new, infile=true)
end

if any(~ismember(cc, ocol_old.lower))
    changeColnames(ofile, ocol_old, ocol_new, infile=true)
end

% check if some loci should be excluded
if isfield(opts, "remove_loci")
    
    tt = tic;
    cols = colnames(opts.remove_loci);
    new_cols = cols.lower;
    assert(all(ismember(["snp", "chr", "pos"], new_cols)), "snp/chr/pos must be present in 'remove_loci' table!")
    idx = ismember(new_cols, ["snp", "chr", "pos"]);
    cols(~idx) = []; new_cols(~idx) = [];
    opts.remove_loci = opts.remove_loci(:, cols);
    opts.remove_loci = renamevars(opts.remove_loci, cols, new_cols);


    % remove from user_SNPsToKeep
    if isfield(opts, "user_SNPsToKeep")
        opts.user_SNPsToKeep = setdiff(opts.user_SNPsToKeep, opts.remove_loci.snp);

        if isempty(opts.user_SNPsToKeep)
            fprintf("\tWARNING: no variant remained in user_SNPsToKeep after excluding remove_loci SNPs!\n")
            opts = rmfield(opts, "user_SNPsToKeep");
        end
    end

    % find variants in LD in European subset of UKBB
    rl = opts.remove_loci;
    
    % fetch variants within distance of 'MR_pruning_dist'
    rl.chr = string(rl.chr);
    rl.range = max(1, (rl.pos - opts.MR_pruning_dist*1e3)) + "-" + (rl.pos + opts.MR_pruning_dist*1e3);

    rng(123);
    ldsample = getQCEID(3, false); 
    ldsample = randsample(ldsample, 5e4); % 50,000 white Brits
    rl.remove = cell(height(rl), 1);
    for k = 1:height(rl)

        if opts.r2 == 0 % remove all variants
            s = struct;
            s.snp = rl.range(k);
            s.chr = string(rl.chr(k));
            check_ld = bgireader(s, table=true, verbose=false, bgenhome=opts.bgenhome);
            check_ld = vertcat(check_ld{:});
            rl.remove{k} = check_ld.snp;

        else
            check_ld = getInsampleLD(rl.range(k), rl.chr(k), ldsample=ldsample, ...
                parallel=true, chunk=70, range=true, bgenhome=opts.bgenhome);
            check_ld.ld = check_ld.ld.^2;
            idx = check_ld.snp == rl.snp(k);
            if ~any(idx)
                fprintf("\t\tfailed to find variant %s in geno files!\n", rl.snp(k))
            else
                check_ld.ld = check_ld.ld(:, idx);
                in_ld_idx = check_ld.ld >= opts.r2; % stringent criteria
                rl.remove{k} = check_ld.snp(in_ld_idx);
            end
        end

        fprintf("\t%d variants (R2 >= %.3f) to be removed for locus %s (%d of %d)\n",...
            numel(rl.remove{k}), opts.r2, rl.chr(k)+":"+rl.snp(k), k, height(rl))
        clear check_ld
    end
    
    opts.remove_loci = rl.remove;
    opts.remove_loci(cellfun(@isempty, opts.remove_loci)) = [];
    opts.remove_loci = unique(vertcat(opts.remove_loci{:}));
    clear rl

   if opts.remove_from == "outcome"
       rmfile = ofile;
   else
       rmfile = efile;
   end
    
   filterFileRows(rmfile, opts.remove_loci, "SNP", infile=true, keep=false)
   tt = toc(tt);
   tt = duration(0, 0, tt);
   fprintf('\tremoved %d SNPs in %s (hr:min:sec)\n', numel(opts.remove_loci), tt)

end


% write to an R script
rmfis = ["remove_from", "remove_loci", "r2", "bgenhome"];
for k = 1:numel(rmfis)
    if isfield(opts, rmfis(k))
        opts = rmfield(opts, rmfis(k));
    end
end


opts.ld = replace(opts.ld, filesep, "/");
opts.hm3 = replace(opts.hm3, filesep, "/");
R = string;
R(1) = "library(MRlap)";
R(numel(R) + 1) = "res = MRlap(exposure='" + efile + "', " + ...
    "outcome='" + ofile + "', " + struct2rfun(opts) + ")";

output = fullfile(pwd, "res" + tmp_name + ".txt");
R(numel(R) + 1) = "write.table(t(unlist(res)), '" + ...
    replace(output, filesep, "/") + "', sep='\t', quote=F, " + ...
    "row.names=F, col.names=T)";
MATLAB2Rconnector("mrlap" + tmp_name, code=R);

% read ldsc log files
lfiles = [opts.exposure_name + ".sumstats.gz_" + ...
    opts.outcome_name + ".sumstats.gz_ldsc.log", ...
    opts.exposure_name + "_munge.log", opts.outcome_name + "_munge.log"];

log_res = struct;
for k = 1:numel(lfiles)
    if ~isfile(lfiles(k)), continue; 
    end
    
    if k == 1
        fi = "lds";
    elseif k == 2
        fi = "exposure";
    else
        fi = "outcome";
    end

    log_res.(fi) = readlines(lfiles(k));
    delete(lfiles(k));

end

% read results
df = readtable(output, TextType="string", VariableNamingRule="preserve");
fi = unique(extractBefore(colnames(df), "."));
res = table;
for k = 1:numel(fi)
    idx = startsWith(colnames(df), fi(k) + ".");
    tmp = df(:, idx);
    cols = colnames(tmp);
    cols = cols.erase(textBoundary("start") + fi(k) + ".");
    tmp = renamevars(tmp, colnames(tmp), cols);

    % merge IVs
    idx = contains(colnames(tmp), ...
        textBoundary("start") + "IVs" + digitsPattern + textBoundary("end"));
    if any(idx)
        merged_ivs = join(tmp{:, idx}, "; ");
        tmp(:, idx) = [];
        tmp.IVs = merged_ivs;
    end

    if fi(k).lower == "mrcorrection"
        res = tmp;
    else % add to userdata
        res.Properties.UserData.(fi(k)) = tmp;
    end

end

% add logs
if ~isempty(fieldnames(log_res))
    res.Properties.UserData.log = log_res;
end

arrayfun(@delete, [efile, ofile, output])

%% UNDER DEVELOPMENT
% % read and keep only matching snps ----------------------------------------
% % under the assumption that all have the same coordinates (GRCh37 or 38)
% 
% % first for cirr1 (exposure 1) and HCC (outcome)
% matchGWAS(gfiles.cirr1, gfiles.hcc)

end % END

% %% subfunctions ===========================================================
% function matchGWAS(file1, file2)
% % keeps matching snps/chr:pos:a1:a2
% 
% [~, name1, ext1] = fileparts(file1);
% df1 = tabularTextDatastore(file1, FileExtensions=ext1, TextType="string",...
%     VariableNamingRule="preserve", TreatAsMissing=["NA", "TEST_FAIL"]);
% 
% [~, name2, ext2] = fileparts(file2);
% df2 = tabularTextDatastore(file2, FileExtensions=ext2, TextType="string",...
%     VariableNamingRule="preserve", TreatAsMissing=["NA", "TEST_FAIL"]);
% 
% df1.VariableNames = createGWASheader(df1.VariableNames);
% df2.VariableNames = createGWASheader(df2.VariableNames);
% 
% % MRlap names
% df1.VariableNames(ismember(df1.VariableNames, "BP")) = {'pos'};
% df2.VariableNames(ismember(df2.VariableNames, "BP")) = {'pos'};
% 
% df1 = tall(df1); df2 = tall(df2);
% 
% cols = intersect(colnames(df1), colnames(df2));
% df1 = df1(:, cols); df2 = df2(:, cols);
% 
% name = getRandomName("", 5);
% fastWriteTable(df1, output="exposure" + name + ".txt", delimiter="\t");
% fastWriteTable(df2, output="outcome" + name + ".txt", delimiter="\t");
% gzip("exposure" + name + ".txt"); delete("exposure" + name + ".txt")
% gzip("outcome" + name + ".txt"); delete("outcome" + name + ".txt")
% 
% 
% % df1.id = df1.CHR + ":" + df1.BP;
% % df2.id = df2.CHR + ":" + df2.BP;
% % 
% % % gather IDs (chr:pos:a1:a2)
% % map1 = gather(df1(:, ["id", "A1", "A2"]));
% % map2 = gather(df2(:, ["id", "A1", "A2"]));
% % 
% % [fidx1, fidx2] = ismember(map1.id, map2.id); fidx2(fidx2 < 1) = [];
% % map1 = map1(fidx1, :); map2 = map2(fidx2, :);
% % 
% % % match alleles
% % idx = map1.A1 ~= map2.A1;
% % if any(idx)
% %     idx = find(idx);
% %     idx2 = map1.A1(idx) == map2.A2(idx); % to be flipped
% %     idx_flip = idx(idx2);
% %     idx(idx2) = []; 
% % end
% % map1(idx, :) = []; map2(idx, :); % remove non-matching alleles
% 
% end % END
