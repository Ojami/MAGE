function res = MRrunner(opts)
% wrapper around TwoSampleMR and MendelianRandomization R pacakges.
% Oveis Jamialahmadi, University of Gothenburg, May 2023.

% Docs: https://mrcieu.github.io/TwoSampleMR
% MR-PRESSO: https://www.nature.com/articles/s42003-023-05408-7
% TODO:
% Isq <- function(y,s)
% {
% 	k <- length(y)
% 	w <- 1/s^2
% 	sum.w <- sum(w)
% 	mu.hat <- sum(y*w)/sum.w
% 	Q <- sum(w*(y-mu.hat)^2)
% 	Isq <- (Q - (k-1))/Q
% 	Isq <- max(0,Isq)
% 	return(Isq)
% }

% @27MARCH: 'harmonizedOnly' option was added to only hamornize variants.
% Doesn't do MR analysis and returns a table of harmonized variants for exp
% and outcome.
% 
% @05SEPT2024: https://github.com/MRCIEU/TwoSampleMR/issues/125: mr_ivw is
% indeed multiplicative random effects model and not fixed-effect
% (mr_ivw_fe)

% Notes: 
%   1- https://github.com/qingyuanzhao/mr.raps
%   
%   2-One important limitation of these methods is the assumption that all
%   valid IVs estimate the same causal effect. Particularly for complex
%   exposures such as BMI, it is possible that different genetic variants
%   have different ratio estimates not because they are invalid IVs, but
%   because there are different ways of intervening on BMI that lead to
%   different effects on the outcome. This can be remedied somewhat in
%   methods based on the IVW method by using a random-effects model (Bowden
%   et al., 2017), or in the contamination mixture method, where causal
%   effects evidenced by different sets of variants will lead to a
%   multimodal likelihood function, and potentially a confidence interval
%   that consists of more than one region.
%   https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22295
% 
%   3-To interpret distortion test of MR-PRESSO: https://www.nature.com/articles/s41598-020-76361-2
%   non-significant result means unbiased estimates after outlier
%   correction.
% 
%   4-To use MR-PRESSO with IVW in presence of heterogeneity see the
%   methods of this paper: https://www.nature.com/articles/s41586-021-03767-x
% 
%   5-Other papers on how to handle presence of heterogeneity:
%   - https://www.nature.com/articles/s41467-023-38234-w
%   - https://www.nature.com/articles/s41588-020-00757-z
%   - https://www.nature.com/articles/s41467-020-14389-8
%   - https://www.nature.com/articles/s41467-019-12536-4
%   - https://www.nature.com/articles/s41588-023-01596-4: MR–Egger
%     regression intercept and MR heterogeneity tests were conducted as
%     additional sensitivity analyses. In case of significant heterogeneity,
%     the MR–pleiotropy residual sum and outlier global test was used to
%     remove genetic variants based on their contribution to heterogeneity.
% 
% TO DO:
%   @07NOV2024: MRlap to correct for sample overlap needs full GWAS summary
%   stat and doesn't work with lead variants only.
        % df = dat[dat$mr_keep, ]
        % cols_e <- grep("\\.exposure$", names(df), value = TRUE) 
        % cols_o <- grep("\\.outcome$", names(df), value = TRUE) 
        % dfe <- df[, c("SNP", "exposure", cols_e)]
        % dfo <- df[, c("SNP", "outcome", cols_o)]
        % names(dfe) <- sub("\\.exposure$", "", names(dfe))
        % names(dfo) <- sub("\\.outcome$", "", names(dfo))
        % 
        % # rename columns
        % new_names <- c("SNP"="snp", "effect_allele"="a1", "other_allele"="a2", "samplesize"="N")
        % dfe = dfe %>% rename_with(~ new_names[.x], .cols = intersect(names(dfe), names(new_names)))
        % dfo = dfo %>% rename_with(~ new_names[.x], .cols = intersect(names(dfo), names(new_names)))
% 
arguments
    opts.package {mustBeMember(opts.package, ["TwoSampleMR", "MendelianRandomization", "MVMR"])} = "TwoSampleMR"
    opts.exposure {mustBeA(opts.exposure, "table")}% a table of exposure summary stats
    opts.outcome {mustBeA(opts.outcome, "table")}% a table of outcome summary stats
    opts.output {mustBeTextScalar} = "MRresults" % path for where to save the output

    % TwoSampleMR arguments
    opts.exposure_name {mustBeTextScalar} % exposure name
    opts.outcome_name {mustBeTextScalar} % outcome name
    opts.bt (1,1) logical = true % is outcome_name binary?
    opts.action (1,1) {mustBeMember(opts.action, [1, 2, 3])} = 2 % choose 1 if you want to keep all SNPs (e.g. when you are confident and already harmonized).

    opts.harmonizeOnly (1,1) logical = false

    % MVMR package: https://wspiller.github.io/MVMR/articles/MVMR.html
    opts.phenocov {mustBeA(opts.phenocov, "table")} % table of correlation coefficients between exposures (if all come from the same sample). Column names should match 'exposure' Pheno column
    
    % MR-PRESSO
    opts.NbDistribution (1,1) double = 1000 % try larger values with larger # IVs
end

% check input tables
ti = ["exposure", "outcome"];
for k = 1:numel(ti)
    if ~(isfield(opts, ti(k)) && istable(opts.(ti(k))))
        continue
    end

    cols = colnames(opts.(ti(k)));

    if ~isfield(opts, ti(k) + "_name")
        idx = ismember(lower(cols), ["pheno", "trait"]); % gwasrunner output
        if any(idx)
            opts.(ti(k) + "_name") = opts.(ti(k)).(find(idx,1))(1);
        else
            opts.(ti(k) + "_name") = ti(k);
        end
    end
    
    % check if the tables are from gwasrunner output
    if any(cols == "P-adj")
        if any(cols == "P")
            opts.(ti(k)).P = [];
        end
    end

    % check if ORs should be converted to beta (only for gwasrunner output)
    cols = colnames(opts.(ti(k)));
    if any(cols == "N.case") && ~any(cols.lower == "beta")
        idx = ~isnan(opts.(ti(k)).("N.case"));
        bidx = find(ismember(cols, ["β|OR", "OR", "β", "OR|β"]), 1);
        opts.(ti(k)).beta = opts.(ti(k)).(bidx);
        opts.(ti(k)).beta(idx) = log(opts.(ti(k)).beta(idx));

        if ti(k) == "outcome"
            if all(idx)
                opts.bt = true;
            elseif all(~idx)
                opts.bt = false;
            end
        end
    end
end

% initialize inputs
R = string;
if opts.bt
    opts.bt = "T";
else
    opts.bt = "F";
end

[wd, name] = fileparts(opts.output);
if wd == "", wd = pwd; end
if ~isfolder(wd), mkdir(wd); end
opts.output = fullfile(wd, name);
file = fullfile(wd, getRandomName("MRrunner", 5)); % random tmp file name

% check if there is only 1 IV (Wald test only)
if height(opts.exposure) == 1
    opts.wonly = true;
else
    opts.wonly = false;
end

%% ========================================================================
if any(ismember(opts.package, ["TwoSampleMR", "MVMR"]))

    R(1) = "library(TwoSampleMR)";
    if opts.package == "MVMR"
        R(2) = "library(MVMR)";
    end
    
    % exposure
    if isfield(opts, 'exposure')
        opts.exposure.Properties.VariableNames = createGWASheader(colnames(opts.exposure), method="TwoSampleMR");
    else
        error("fix this using https://mrcieu.github.io/TwoSampleMR/articles/exposure.html")
    end
    
    if opts.package == "MVMR"
        opts.exposure = renamevars(opts.exposure, "Pheno", "Phenotype");
    else
        opts.exposure.Phenotype(:) = opts.exposure_name;
    end
    writetable(opts.exposure, file + ".exp.txt", Delimiter="\t")
    R(3) = "exp_data = as.data.frame(data.table::fread('" + file.replace(filesep, "/") + ".exp.txt', sep = '\t'))";
    R(4) = "exp_data = format_data(exp_data, type = 'exposure')";
    
    % outcome
    if isfield(opts, 'outcome')
        opts.outcome.Properties.VariableNames = createGWASheader(colnames(opts.outcome), method="TwoSampleMR");
    else
        % ao <- available_outcomes()
        % head(subset(ao, select = c(trait, id)))
        % bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
        % chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')
        error("fix this using https://mrcieu.github.io/TwoSampleMR/articles/outcome.html")
    end

    opts.outcome.Phenotype(:) = opts.outcome_name;
    writetable(opts.outcome, file + ".outcome.txt", Delimiter="\t")
    R(5) = "out_data = as.data.frame(data.table::fread('" + file.replace(filesep, "/") + ".outcome.txt', sep = '\t'))";
    R(6) = "out_data = format_data(out_data, type = 'outcome')";

    % Harmonise data
    R(7) = "#harmonisation";
    if opts.package == "MVMR" && ~opts.harmonizeOnly
        R(8) = "exp_data$id.exposure = exp_data$exposure";

        % does not work for passing to MVMR package, ergo all exposure and
        % outcome should be harmonized beforehand
        R(9) = "dat = mv_harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, harmonise_strictness = " + opts.action + ")";
    else
        R(9) = "dat = harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=" + opts.action + ")";
    end
    R(10) = "dat = add_rsq(dat)";

    if opts.harmonizeOnly
        R(11) = "write.table(dat, '" + file.replace(filesep, "/") + ".harmonized.txt', quote = F, sep = '||', row.names = F, col.names = T)";
    else

        if opts.package == "MVMR" 
            R(13) = "df <- format_mvmr(BXGs = dat$exposure_beta, BYG = dat$outcome_beta, seBXGs = dat$exposure_se, seBYG = dat$outcome_se, RSID = rownames(dat$exposure_beta))";

            R(15) = "#Test for weak instruments";
            R(16) = "sres <- strength_mvmr(r_input = df, gencov = 0)"; % genecov must be read from input in case they are from same samples

            R(20) = "#Test for horizontal pleiotropy using conventional Q-statistic estimation";
            R(21) = "pres <- pleiotropy_mvmr(r_input = df, gencov = 0)";

            R(23) = "# Step 5: Estimate causal effects";
            R(24) = "res <- ivw_mvmr(r_input = df)";

            % R(26) = "# Step 6: Robust causal effect estimation";
            % R(27) = "res1 <- qhet_mvmr(df, mvmrcovmatrix, CI = T, iterations = 100)";

            R(29) = "# write to output";
            R(30) = "out = as.data.frame(cbind(res, t(sres)))"; % , res1)";
            R(31) = "out$Qstat = pres$Qstat";
            R(32) = "out$Qpval = pres$Qpval";
            R(33) = "out = cbind(Exposure = dat$expname$exposure, Outcome = dat$outname$outcome, out)";
            R(34) = "write.table(out, '" + file.replace(filesep, "/") + ".res_dat.txt', quote = F, sep = '||', row.names = F, col.names = T)";

        else

            % perform the MR
            R(11) = "mr_methods = mr_method_list()";
    
            %@30MAY2023: mr_raps throws some error, maybe a bug
            if opts.wonly
                R(12) = "mr_methods = 'mr_wald_ratio'";
            else
                R(12) = "mr_methods = union(mr_methods$obj[mr_methods$use_by_default], c('mr_raps', 'mr_penalised_weighted_median', 'mr_ivw_fe'))"; % can be changed later
                % R(12) = "mr_methods = mr_methods$obj[mr_methods$use_by_default]";
            end
            R(13) = "mr_out = mr(dat, method_list=mr_methods)";

            if ~opts.wonly
    
                R(15) = "#Heterogenetiy analysis";
                R(16) = "mr_het = mr_heterogeneity(dat)";
        
                R(18) = "#Horizontal pleiotropy";
                R(19) = "mr_ple = mr_pleiotropy_test(dat)";
                
                R(21) = "#Single SNP analysis";
                R(22) = "mr_single = mr_singlesnp(dat)"; % default method: wald test
        
                R(24) = "#Leave-one-out analysis";
                R(25) = "mr_loo = mr_leaveoneout(dat)";
            

                % save figures
                R(27) = "#Figures";
                R(28) = "p1 = mr_scatter_plot(mr_out, dat)";
                R(29) = "ggplot2::ggsave(p1[[1]], file = '" + opts.output.replace(filesep, "/") + "_scatter.png', width = 7, height = 7)";
        
                R(30) = "p2 = mr_forest_plot(mr_single)";
                R(31) = "ggplot2::ggsave(p2[[1]], file = '" + opts.output.replace(filesep, "/") + "_forest.png', width = 6, height = 9)";
        
                R(32) = "p3 = mr_leaveoneout_plot(mr_loo)";
                R(33) = "ggplot2::ggsave(p3[[1]], file = '" + opts.output.replace(filesep, "/") + "_loo.png', width = 6, height = 9)";
        
                R(34) = "p4 = mr_funnel_plot(mr_single)";
                R(35) = "ggplot2::ggsave(p4[[1]], file = '" + opts.output.replace(filesep, "/") + "_funnel.png', width = 7, height = 7)";
            end
    
            % directionality test
            % needs more evaluation: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html
            % R(37) = "#directionality";
            % R(38) = "mr_dir = directionality_test(dat)";
    
            % combine all results
            if opts.wonly
                R(41) = "write.table(mr_out, '" + file.replace(filesep, "/") + ".res.txt', quote = F, sep = '||', row.names = F, col.names = T)";
                R(42) = "write.table(dat, '" + file.replace(filesep, "/") + ".res_dat.txt', quote = F, sep = '||', row.names = F, col.names = T)";
            else
                R(40) = "res = combine_all_mrresults(mr_out, mr_het, mr_ple, " + ...
                    "mr_single, ao_slc=F, Exp=" + opts.bt + ")";
                R(41) = "write.table(res, '" + file.replace(filesep, "/") + ".res.txt', quote = F, sep = '||', row.names = F, col.names = T)";
        
                % MR-PRESSO
                R(43) = "#MR-PRESSO";
                R(44) = "mp = tryCatch({";
                R(45) = "MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome', " + ...
                    "BetaExposure = 'beta.exposure', SdOutcome = 'se.outcome'," + ...
                    " SdExposure = 'se.exposure', OUTLIERtest = TRUE, " + ...
                    "DISTORTIONtest = TRUE, data = dat, NbDistribution = " + opts.NbDistribution + "," + ...
                    "SignifThreshold = 0.05)";
                R(46) = "}, error = function(e) {NULL})";
                
                R(48) = "if (!is.null(mp)){";
                R(49) = "   if ('Distortion Test' %in% names(mp$`MR-PRESSO results`)){";
                R(50) = "     outSNPs = dat$SNP[mp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`]";
                R(51) = "     outSNPs = paste(outSNPs, collapse = ',')";
                R(52) = "     mp = cbind(mp$`Main MR results`, Global.P = mp$`MR-PRESSO results`$`Global Test`$Pvalue, Distortion.P = mp$`MR-PRESSO results`$`Distortion Test`$Pvalue, Outliers=outSNPs)";
                R(53) = "   }else{";
                R(54) = "     mp = cbind(mp$`Main MR results`, Global.P = mp$`MR-PRESSO results`$`Global Test`$Pvalue)}";
                R(55) = "   write.table(mp, '" + file.replace(filesep, "/") + ".res_mp.txt', quote = F, sep = '||', row.names = F, col.names = T)";
                R(56) = "}";
                R(57) = "write.table(dat, '" + file.replace(filesep, "/") + ".res_dat.txt', quote = F, sep = '||', row.names = F, col.names = T)";

            end
        end
    end

    MATLAB2Rconnector(file + ".r", code=R, delr=true, log=false);

    if opts.harmonizeOnly
        res = readtable(file + ".harmonized.txt", TextType="string", ...
            VariableNamingRule="preserve", NumHeaderLines=0, Delimiter="||", ...
            TreatAsMissing="NA");

        % delete residual files
        delete(file + ".harmonized.txt")
    else
        
        if opts.package == "TwoSampleMR"
            res = readtable(file + ".res.txt", TextType="string", ...
                VariableNamingRule="preserve", NumHeaderLines=0, Delimiter="||", ...
                TreatAsMissing="NA");
            res(:, ["id.exposure", "id.outcome"]) = [];

            %@05SEPT2024: see comment above
            res.Method(res.Method == "Inverse variance weighted") = "Inverse variance weighted (multiplicative random effects)";
            
            if isfile(file + ".res_mp.txt")
                mp = readtable(file + ".res_mp.txt", TextType="string", ...
                    VariableNamingRule="preserve", NumHeaderLines=0, Delimiter="||", ...
                    TreatAsMissing="NA");
                mp.Exposure(:) = res.exposure(1);
            end
    
            dat = readtable(file + ".res_dat.txt", TextType="string", ...
                VariableNamingRule="preserve", NumHeaderLines=0, Delimiter="||", ...
                TreatAsMissing="NA");
            dat(:, ["id.exposure", "id.outcome"]) = [];

            % equivalent to (N - 2).*(r2./(1-r2))
            % https://www.mendelianrandomization.com/images/MRbook_ch7_sample.pdf
            dat.F = (dat.("beta.exposure")./dat.("se.exposure")).^2;
    
            writetable(res, opts.output + ".xlsx", Sheet="MR", WriteMode="overwritesheet")
            if isfile(file + ".res_mp.txt")
                writetable(mp, opts.output + ".xlsx", Sheet="MR-PRESSO", WriteMode="overwritesheet")
                delete(file + ".res_mp.txt")
            end
            writetable(dat, opts.output + ".xlsx", Sheet="metadata", WriteMode="overwritesheet")
            delete(file + ".res.txt")
            delete(file + ".res_dat.txt")
            
        elseif opts.package == "MVMR"
            res = readtable(file + ".res_dat.txt", TextType="string", ...
                VariableNamingRule="preserve", NumHeaderLines=0, Delimiter="||", ...
                TreatAsMissing="NA");
            writetable(res, opts.output + ".xlsx", Sheet="MVMR", WriteMode="overwritesheet")
            delete(file + ".res_dat.txt")
        end

    end

    delete(file + ".exp.txt")
    delete(file + ".outcome.txt")

%% ====================================================================
elseif opts.package == "MendelianRandomization"
end


end % END