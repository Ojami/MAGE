function res = mediate(opts)
% Written by R2MATLABWrapper on 10-May-2023 performs mediation analysis
% using R packages:mediation and MeMoBootR (for serial mediation).
% 
% Oveis Jamialahmadi, University of Gothenburg, May 2023.
% 
% "mediation" package terminology:
% ADE: direct effect of X on Y after taking into account a mediation
%      (indirect) effect.
% ACME: the mediation effect (ACME) is the total effect minus the direct
%       effect.
% 
% @08MARCH2024: TO DO: add "robmed" and "manymome" R packages
% @08MARCH2024: When 'treat' is non-binary, control.value (default 0) and
%               treat.value (default 1) should be set:
%               cran.r-project.org/web/packages/mediation/mediation.pdf
% @08MARCH2024: 'interaction' flage was added (default: false) to add the
%               interaction of treat*mediator to the model. Also calculates
%               treatment-mediator interaction: AMCE(treat_value) -
%               ACME(control_value)
% @04APR2024: 'sens' flag was added (default: true). Runs sensitivity
%              analysis using medsens() function and returns rho at which
%              ACME for Control/case group is equal to 0. This tests
%              sequential ignorability assumption (possible existence of
%              unobserved pre-treatment covariates). If there exist
%              unobserved pre-treatment confounders which affect both the
%              mediator and the outcome, we expect that the sequential
%              ignorability assumption is violated and ρ is no longer zero.
%              
%              If the ρ value is small in magnitude, it means that the ACME
%              estimates obtained in the previous step would be reversed if
%              the errors for the mediation and outcome models are just
%              weakly correlated. In other words, a small value indicates
%              the ACME estimates may not be robust to unobserved
%              confounders. The usefulness of this measure is that it
%              quantifies the sensitivity to such confounding, the presence
%              of which, as mentioned earlier, is untestable with observed
%              data.
%              https://www.frontiersin.org/articles/10.3389/feduc.2022.886722/full
%              https://www.cambridge.org/core/journals/political-analysis/article/identification-and-sensitivity-analysis-for-multiple-causal-mechanisms-revisiting-evidence-from-framing-experiments/E8EDFE67C8776530FDB445BCA3EC5968


arguments
    % mediation::mediate
    opts.treat {mustBeTextScalar} % should be a valid name in 'tab'
	opts.mediator {mustBeTextScalar} % should be a valid name in 'tab'
    opts.outcome {mustBeTextScalar} % should be a valid name in 'tab'. If 'tab' is empty, this is equivalent to 'trait' argument of getPhenoCov func.
	opts.sims (1,1) double = 1000
	opts.boot (1,1) logical = true % use non-parametric bootstrapping
	opts.boot_ci_type {mustBeMember(opts.boot_ci_type, ["perc", "bca"])} = "perc" 
	opts.covariates
	opts.control
	opts.conf_level (1,1) double = 0.95
	opts.control_value (1,1) double = 0
	opts.treat_value (1,1) double = 1
	opts.long (1,1) logical = true
	opts.dropobs (1,1) logical = false
	opts.robustSE (1,1) logical = true 
	opts.cluster
	opts.group_out
	opts.use_speed (1,1) logical = false
    opts.sens (1,1) logical = true % to run medsens 
    
    % MeMoBootR::mediation2 (PROCESS model 6 for serial mediation)
    opts.mediator2 {mustBeTextScalar} % should be a valid name in 'tab'

	% non-native arguments for output
    opts.tab {mustBeA(opts.tab, "table")} % table of data required for mediation analysis. If left emtpy, uses 'getPhenoCov' to fetch the data tabe
    opts.method {mustBeMember(opts.method, ["mediation", "MeMoBootR"])} = "mediation" % package name to be used
    opts.interaction (1,1) logical = false % interaction between treat and mediator

    % getPhenoCov arguments (used when 'tab' is empty)
    opts.qc (1,1) double {mustBeMember(opts.qc, 0:6)} = 1 % see getQCEID function for more details
    opts.qcinc (:, 1) double {mustBeNumeric} = [] % custom eids to be included; in this case intersect of qc and qcinc will be used for analyses.
    opts.phenopath {mustBeFolder} 
    opts.transform {mustBeMember(opts.transform, ["none", "log", "log10", "irnt"])} = "irnt" % transformation function for 'outcome'
    opts.covar {mustBeText} = ["Age", "Sex", "BMI"]; % list of covariates, it automatically adds PC1-10 and array batch.
    opts.wes (1, 1) logical = false % for WES analysis (array batch is removed from covariates)
    opts.verbose (1,1) logical = true 
    opts.covarTab {mustBeA(opts.covarTab, "table")} % additional covariates to be added (a table with covariates + 'eid' column).
    opts.defaultcovar (1,1) logical = true % to automatically add PC1-10 and array batch to covariates.
    opts.snp
    opts.chr
    opts.model {mustBeMember(opts.model, ["add", "dom", "rec"])} = "add"
    opts.genohome {mustBeFolder} % must only contains unique BGEN files for each chr
    opts.dosage (1,1) logical = true
end

% check necessary inputs 'treat', 'mediator', 'mediator2', and 'outcome'
checkfis = ["treat", "mediator", "mediator2", "outcome"];
reqfis = strings(numel(checkfis), 1);
for k = 1:numel(checkfis)
    if checkfis(k) == "mediator2" && opts.method == "mediation"
        continue
    end

    if ~isfield(opts, checkfis(k))
        if checkfis(k) == "mediator2" % mediator2 is optional
            continue
        else
            error("%s argument is missing!", checkfis(k))
        end
    end
    reqfis(k) = checkfis(k);
end
reqfis(reqfis == "") = []; % final required inputs

% check input table
if isfield(opts, 'tab')
    % fix this
    error("must be implemented!")
else
    for k = 1:numel(reqfis)
        if reqfis(k) == "outcome", continue; end
        opts.covar = union(opts.covar, opts.(reqfis(k)));
    end
    gopts = opts;
    rfis = ["qc", "qcinc", "phenopath", "snp", "chr", "genohome", "model", ...
        "transform", "covar", "wes", "verbose", "defaultcovar", "covarTab",...
        "dosage"];
    fis = setdiff(fieldnames(gopts), rfis);
    gopts = rmfield(gopts, fis);
    
    if isfield(gopts, "snp")
        gopts.covar = setdiff(gopts.covar, gopts.snp);
    end
    
    if isfield(gopts, "covarTab")
        gopts.covar = setdiff(gopts.covar, colnames(gopts.covarTab));
    end
    
    gopts.legacy = false; % use 'trait' file names for the column header
    gopts.trait = opts.outcome;
    gopts = namedargs2cell(gopts);
    [xt, xc] = getPhenoCov(gopts{:});
    opts.tab = join(xt, xc, Keys="eid");
    opts.tab.eid = [];
    opts.covar = colnames(opts.tab);

    % remove missing values from snp
    if isfield(opts, "snp")
        opts.tab(opts.tab.(opts.snp) < 0, :) = [];
    end

    % apply the same transformation function of outcome to mediator
    opts.tab.(opts.mediator) = feval(opts.transform, opts.tab.(opts.mediator));
    clear xt xc
end

% initialize R specific arguments
ropts = rmfield(opts, [rfis, "method", "tab", "interaction", "sens"]);

% check if required inputs are present in the table
for k = 1:numel(reqfis)
     assert(any(colnames(opts.tab) == opts.(reqfis(k))))
end

% remove required fields from the covariate list 
for k = 1:numel(reqfis)
    opts.covar = setdiff(opts.covar, opts.(reqfis(k)));
end

% clean the table
 opts.tab = rmmissing(opts.tab);
 if numel(unique(opts.tab.(opts.outcome))) <= 2 
     opts.binary = true;
 else
     opts.binary = false;
 end

 % initialize R inputs
file = getRandomName("mediate_R", 5); % random tmp file name
cols = matlab.lang.makeValidName(cellstr(colnames(opts.tab)));
tab = table2array(opts.tab);

for k = 1:numel(reqfis)
    opts.(reqfis(k)) = matlab.lang.makeValidName(opts.(reqfis(k)));
    ropts.(reqfis(k)) = matlab.lang.makeValidName(ropts.(reqfis(k)));
end

save(file + ".mat", "cols", "tab")
clear tab

% write the R script
r(1) = "require(mediation)";
r(2) = "data = R.matlab::readMat('" + file + ".mat" + "')";
r(3) = "hdr = unlist(data$cols, use.names = F)";
r(5) = "data = as.data.frame(data$tab)";
r(6) = "colnames(data) = hdr";

if opts.binary
    r(7) = "data$" + opts.outcome + " = as.factor(data$" + opts.outcome + ")";
end

opts.int = "+";
if opts.interaction, opts.int = "*"; end

% set up regression models
if opts.method == "mediation"
    r(9) = "#set up regression models";
    r(10) = "med.fit = lm(" + opts.mediator + "~" + opts.treat + "+" + ...
        join(opts.covar, "+") + ", data=data)";

    if opts.binary
        r(11) = "out.fit = glm(" + opts.outcome + "~" + opts.treat + opts.int ...
            + opts.mediator + "+" + join(opts.covar, "+") + ...
            ", data=data, family = binomial())";
    else
        r(11) = "out.fit = lm(" + opts.outcome + "~" + opts.treat + opts.int ...
            + opts.mediator + "+" + join(opts.covar, "+") + ...
            ", data=data)";
    end

    underFis = ["model_m,model.m", "model_y,model.y",...
        "boot_ci_type,boot.ci.type", "conf_level,conf.level", ...
        "control_value,control.value", "treat_value,treat.value",...
        "group_out,group.out"];
    
    if isfield(ropts, "mediator2")
        ropts = rmfield(ropts, "mediator2"); % "MeMoBootR" arguments
    end

    r(12) = "out = mediate(med.fit, out.fit, " + struct2rfun(ropts, replace=underFis) + ")";
    r(13) = "x = summary(out)";
    
    % print.summary.mediate
    r(14) = "clp <- 100 * x$conf.level";
    if opts.binary || opts.interaction
        r(15) = "smat <- c(x$d0, x$d0.ci, x$d0.p)";
        r(16) = "smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))";
        r(17) = "smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))";
        r(18) = "smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))";
        r(19) = "smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))";
        r(20) = "smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))";
        r(21) = "smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))";
        r(22) = "smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))";
        r(23) = "smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))";
        r(24) = "smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))";
        r(25) = "smat <- cbind(c('ACME (control)', 'ACME (treated)'," + ...
            "'ADE (control)', 'ADE (treated)', 'Total Effect', " + ...
            "'Prop. Mediated (control)', 'Prop. Mediated (treated)', " + ...
            "'ACME (average)','ADE (average)','Prop. Mediated (average)'), smat)";
    else
        r(15) = "smat <- c(x$d1, x$d1.ci, x$d1.p)";
        r(16) = "smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))";
        r(17) = "smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))";
        r(18) = "smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))";
        r(19) = 'smat <- cbind(c("ACME", "ADE", "Total Effect", "Prop. Mediated"), smat)';
    end

    r(26) = 'colnames(smat) <- c("Term", "Estimate", paste(clp, "% CI Lower", sep=""), paste(clp, "% CI Upper", sep=""), "P")';
    r(27) = "write.table(smat, '" + file + ".txt', quote = F, row.names = F, sep = '||')";

    if opts.interaction % test of significance between ACME(0) and 1: nuell is there is no difference
        r(29) = "bint = test.TMint(out,conf.level=.95)";
        r(30) = "bint = c(bint$statistic, bint$conf.int, bint$p.value)";
        r(31) = 'names(bint) = c("Estimate", "95% CI Lower", "95% CI Upper", "P")';
        r(32) = "write.table(rbind(bint), file='" + file + ".int.txt', row.names = F, col.names = T, quote = F, sep = '||')";
    end

    if opts.sens
        r(34) = "sens = medsens(out, rho.by = 0.2, effect.type = 'indirect', sims = " + opts.sims + ")";
        r(35) = "write(sens$err.cr.d, file='" + file + ".sens.txt', sep = '||')";
    end

elseif opts.method == "MeMoBootR"
    if isempty(opts.covar) || all(opts.covar == "") 
        cvs = "NULL";
    else
        cvs = "c('" + join(opts.covar, "','") + "')";
    end
    
    if isfield(opts, "mediator2")
        r(9) = "out = MeMoBootR::mediation2('" + opts.outcome + "', '" + opts.treat + ...
            "', '" + opts.mediator + "', '" + opts.mediator2 + "', cvs = " + ...
            cvs + ", data, with_out = T, nboot = " + opts.sims + ...
            ", conf_level = " + opts.conf_level + ")";
        r(10) = "smat = data.frame(Term = c('M1', 'M2', 'M12'), Estimate = out$boot.results$t0, " + ...
            "SE = apply(out$boot.results$t, 2, sd, na.rm=T), " + ...
            "CI = t(sapply(out$boot.ci, function(x) x$normal[2:3])))";
    else
        r(10) = "out = MeMoBootR::mediation1('" + opts.outcome + "', '" + opts.treat + ...
            "', '" + opts.mediator + "', cvs = " + ...
            cvs + ", data, with_out = T, nboot = " + opts.sims + ...
            ", conf_level = " + opts.conf_level + ")";
        r(12) = "smat = data.frame(Term = c('M1', 'M2', 'M12'), Estimate = out$boot.results$t0, " + ...
            "SE = apply(out$boot.results$t, 2, sd, na.rm=T), " + ...
            "CI = t(out$boot.ci$normal[2:3]), P = out$p.value)";
    end
    
    r(13) = "write.table(smat, '" + file + ".txt', quote = F, row.names = F, sep = '||')";

end

MATLAB2Rconnector(file + ".r", code=r, delr=true);

res = readtable(file + ".txt", TextType="string", ...
    VariableNamingRule="preserve", Delimiter="||", TreatAsMissing="NA");

if opts.interaction
    res_int = readtable(file + ".int.txt", TextType="string", ...
        VariableNamingRule="preserve", Delimiter="||", TreatAsMissing="NA");
    res_int.Term = "ACME(1)-ACME(0) interaction";
    res_int = movevars(res_int, "Term", Before=1);
    res = [res; res_int];
    delete(file + ".int.txt")
end

if opts.sens
    res_sens = readmatrix(file + ".sens.txt", Delimiter="||");
    delete(file + ".sens.txt")
    % res = [res; res(end, :)];
    % res.Term(end) = "R";
    % res{end, ["95% CI Lower", "95% CI Upper", "P"]} = nan;
    res.Properties.UserData = res_sens;
end

% tcdf(mdl.Coefficients.tStat(2), 20, "upper")*2
% out.CI = [out.(bcol) - out.SE.*tinv(1-0.05/2, dfe),...
%         out.(bcol) + out.SE.*tinv(1-0.05/2, dfe)];

delete(file + ".mat")
delete(file + ".txt")
end % END
