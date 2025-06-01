function res = fitSurvCox(preds, DV, covariates, opts2, opts)
% called by gwasrunner for survival analysis
% May 2024, University of Gothenburg, Oveis Jamialahmadi.

% Pipelines:
%   - https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
%   - https://stats.stackexchange.com/questions/648473/time-to-event-analysis-with-left-truncation-and-right-censoring-depending-on-the
%
% Summary of possible approaches (using a combination of "timeOrigin",
%   "LT", and "prevalent"):
%   In the following, tt is time since enrollment to diagnosis and can be
%   < 0 for prevalent cases. tt0 is age at diagnosis. In case of
%   time-to-event, user must define tt as time between diagnosis to event
%   , and tt0 as delayed time (W, see "TTE" when "LT" is true). Possible
%   scenarios one can define are:
%       1- age-of-onset (age-of-diagnosis) with left-truncated data
%       (timeOrigin="AOO", LT=true, prevalent=false): for a single
%       event/diagnosis and appropriate for genetic variants. Note that in
%       this case, prevalent cases are removed (prevalent=false).
% 
%       2- similar to above but without left-truncation: adjust for age at
%       enrollment, but (1) is preferred.
% 
%       3- time-to-event with left-truncation (timeOrigin="TTE", LT=true,
%       prevalent=true): when following up from diagnosis to death for an
%       event (origin is diagnosis time). In this case tt0 and tt should be
%       realigned to date of diagnosis (currently phenoParser does not
%       support this, but can be easily defined using two calls to
%       phenoParser).
% 
%       4- similar to above but without left-truncation:
%           - if origin is time from diagnosis, adjust for age at
%           diagnosis, but (3) is preferred.
%           - if origin is time from enrollment (most common and simplest
%           approach), adjust for age at enrollment, but (1) is preferred.
% 
%   IMPORTANT: "TTE" and "LT" together ONLY means time from diagnosis to
%              event. However, "TTE" without "LT" can also be time to event
%              from enrollment. "AOO" and "LT" together, sets "prevalent"
%              to false. Note that, irrespective of choice of "timeOrigin"
%              and "prevalent", all tt and tt0 <= 0 are excluded.
% 
%   "timeOrigin":
%   - "TTE": (default) time-to-event. It uses time-to-event (tt) by
%   subtracting enrollment date and date of event, or time between a
%   diagnosis or event (tt). This is the most common approach used in
%   published papers with no left-truncation ("LT" = false).
%       - If "LT" is true (left-truncated) time-to-event when the interest
%       is studying the time between a diagnosis (e.g. HCC) to death
%       (because of HCC). In this scenario, some individuals may be
%       diagnosed before enrollment. Note that in this case, any data (such
%       as covariates measured later) cannot be used to adjust the model.
%       So, we align with left-truncation to date of diagnosis. We define
%       W, so that if someone is diagnosed before enrollment (W > 0: time
%       before enrollment) and if diagnosed after enrollment (W = 0). Let,
%       Y be the time from diagnosis to death, then we fit this model:
%               
%           coxph(Surv(W, Y, status, type="counting") ~ exposure
% 
%           Consider the following, ▲: diagnosis, O: enrollment, X: death.
%           We align all at time of diagnosis (e.g HCC). W for first
%           individual is 0.
%           O[1950]-------------▲[1990]_______X[2000] {tt0=W=0, tt=Y=10}
%           ▲[1970]---------O[1990]_____X[2005] {tt0=W=20, tt=Y=35}
% 
% 
%   - "AOO": age-of-onset. It uses the age-of-onset (tt0) which is similar to
%   approaches used in GWAS tools such as GATE and SpaCox. While there is
%   no need to adjust for age, this analysis ignores left-truncation and is
%   not appropriate for relatively old participants (See GATE paper).
%   Recommended to use with "LT".
%       - if "LT" is true date of birth is used as time of origin (t=0). 
%         Note that in this case, any time-dependent covariates measured
%         later cannot be used to adjust the model.
% 
%           Consider the following, ●: birth, ▲: event, O: enrollment
%           ●[1990]-----▲[2002]---O[2010] {LT=NA, AOO=12, tt=-8, tt0=12} --> remove
%           ●[2000]-----------O[2010]____▲[2015] {LT=10, AOO=15, tt=5, tt0=15}
% 
%           Then the model is:
%               coxph(Surv(LT, AOO, status, type="counting") ~ exposure)
% 
%   - "prevalent" is set to false when "LT" is true with "AOO". 
% 
%   All above-mentioned scenarios can simply be summarized as this cross
%   tab (ordered as start, stop):
%                            LT
%                   -------------------
%                       T    |    F
%                   -------------------                
%              TTE   tt0,tt  |   tt
%                   -------------------
%              AOO   tt,tt0  |   tt0
%                    ↓
%                  NA,    <0     
%                 ⎨tt0-tt, >0
% 
% Detailed notes:
% 
%   1- Prevalent cases were defined as cases whose time from date of
%   diagnosis until study entry was >6 months
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2695697/
% 
%   2- Note that all other competing risks are by default censored.
%   Example: in case of CLD, viral hepatitis and death of other causes are
%   competing risks. So, by default, cause-specific HR is calculated. It's
%   possible to  also calculate subdistribution hazards as well by choosing
%   "FGR" option.

%   3- To define a time-to-event phenotype for genome-wide studies, one
%   should consider age-of-onset (AAO), since the exposure. in this case is
%   a genetic variant existing at time of birht. In this case, one should
%   not adjust the model for age:
%              
%              https://www.nature.com/articles/s41467-021-22538-w#Sec9:
%              "We did not include any covariates of age or year of birth
%              because these are directly associated with our phenotypes."
% 
%              Other useful refereces:
%              - https://www.nature.com/articles/s41398-022-01782-8#Sec2
%              - https://github.com/WenjianBI/SPACox
%              - https://www.nature.com/articles/s41467-023-41210-z
%              
%   4- On truncation and censoring:
%              section 3.7.3 Alternate time scales
%              Modeling Survival Data: Extending the Cox Model
% 
%              https://erj.ersjournals.com/content/36/3/549.full
%              -In order to remove survivor bias that results from the
%              inclusion of prevalent patients, survival estimates and a
%              Cox proportional hazards model from time to diagnosis were
%              adjusted for the left truncation arising from the delay
%              between diagnosis and study entry.
%               
%              https://academic.oup.com/aje/article/173/9/1078/122651?login=true
%              - Left censoring occurs if a participant is entered into the
%              study when the milestone of interest occurred prior to study
%              entry but the age at that milestone is unknown. Left
%              truncation occurs when individuals who have already passed
%              the milestone at the time of study recruitment are not
%              included in the study.
% 
%              For left-truncation, age should be accounted for and not
%              used as a covariate:
%               
%               https://www.usu.edu/math/jrstevens/biostat/projects2013/pres_LeftTruncation.pdf
%               - "Age is often used as a covariate when it should be used
%               as a left-truncation point.  When age is used as a
%               left-truncation point, it is unnecessary to use it as a
%               covariate in the model"
% 
%              For left-truncated datasets see chapters 3 and 5:
%                https://xsliulab.github.io/Workshop/2021/week3/survival-analysis-book.pdf
%                - "However, in some circumstances one may wish to compare
%                survival times starting from time from diagnosis, and then
%                it is essential to account for the left truncation."
% 
%              on left-truncation (if using AAO), adjusted for age also as
%              covariate: Surv(enrollment_age, AAO, event)
%              - https://www.sciencedirect.com/science/article/pii/S2666247721000580
%              - https://www.nature.com/articles/s41416-019-0465-y
%              - On age-of-onset: https://www.nature.com/articles/s41467-022-32885-x
% 
%             on case ascertainment and left-truncation, not adjusted for age as covariate:
%             - - Main reference: https://academic.oup.com/biomedgerontology/article/53A/5/M337/588259
%             - https://www.nature.com/articles/s41467-020-17374-3#Sec7
%               "All participants diagnosed at dates prior to enrollment in
%               the UK Biobank were considered prevalent at baseline, while
%               participants diagnosed after enrollment were considered
%               incident"
% 
%             When the outcome is death, left-truncation is underpowered,
%             and if not handled, is biased. Because a subject who died
%             before enrollment could not be in the study, ergo not at
%             risk.
%             - https://www.sciencedirect.com/science/article/pii/S2405844024010065
%               "In survival analysis, left truncation bias [46] was
%               corrected by considering only the time between patient
%               sampling and time of death. Left truncation bias occurs
%               when risk of death is measured over a time in which it
%               could not have occurred; by definition a patient who had
%               already died could not have been recruited into a study.
%               The rationale is to exclude survival time that occurred
%               before sampling, because the patient could not, by
%               definition, have died during this period. Analysis of left
%               truncated data can lead to false positive associations
%               [46], however, left truncation bias correction can lead to
%               under-powered analysis because of the loss of information."


arguments
    preds % {mustBeNumeric,mustBeReal}
    DV {mustBeNumeric,mustBeReal,mustBeVector}
    covariates = []
    opts2 {mustBeA(opts2, "struct")} = struct() % from gwasrunner

    opts.plotFlag (1,1) logical = false
    opts.checkDiagnostics (1,1) logical = false
    opts.logRank (1,1) logical = false % not that doesn't support LT. LT overrides this
    opts.logRankFun {mustBeTextScalar, mustBeMember(opts.logRankFun, ["event", "cumhaz", "pct", "null"])} = "null" % "event" plots cumulative events (f(y) = 1-y), "cumhaz" plots the cumulative hazard function (f(y) = -log(y)), and "pct" for survival probability in percentage.
    opts.FGR {mustBeTextScalar, mustBeMember(opts.FGR, ["none", "crr", "fastcrr"])} = "none" % Fine-Gray for competing risks: note that "fastcrr" doesn't support LT

    opts.prevalent (1,1) logical = false % exclude prevalent cases
    opts.timeOrigin {mustBeMember(opts.timeOrigin, ["TTE", "AOO"])} = "TTE"
    opts.LT (1,1) logical = false % no left-truncation
    opts.landmark (1,1) double = 0 % landmark analysis, individuals with a non-zero censoring status (including competing risks) with a tt <= landmark (yr) are excluded. 

end

if numel(fieldnames(opts2))
    fis = string(fieldnames(opts2));
    for k = 1:numel(fis)
        opts.(fis(k)) = opts2.(fis(k));
        opts2 = rmfield(opts2, fis(k));
    end
end


opts.plotFlag = upper(string(opts.plotFlag));
opts.checkDiagnostics = upper(string(opts.checkDiagnostics));

%@14AUG2024: now fixed it.
% if opts.LT, opts.logRank = false; end

% for R prepare a table
opts.hascovars = false;
if ~isempty(covariates), opts.hascovars = true; end
df = gather([preds, covariates, DV]);
assert(size(preds, 1) == numel(opts.qc_eid), "input eids and covariates must have the same size!")
clear preds covariates DV


% add time-to-event
[f1, f2] = ismember(opts.surv_data.eid, opts.qc_eid);
opts.surv_data = opts.surv_data(f1, :);
opts.qc_eid = opts.qc_eid(f2(f1));

if opts.timeOrigin == "TTE"
    method_cols = ["censor", "tt"];
else % AOO
    method_cols = ["censor", "tt0"];
end

% for left-truncation use start,stop times
if opts.LT
    method_cols = union(method_cols, ["tt", "tt0"]);
end

if ~opts.prevalent, method_cols = [method_cols, "base"]; end
df = [df(f2(f1), :), opts.surv_data{:, [method_cols, "eid"]}];

% conver to a table to be transferred to R
opts.surv_preds.preds = matlab.lang.makeValidName(opts.surv_preds.SNP);
opts.surv_preds.Pheno2 = matlab.lang.makeValidName(opts.surv_preds.Pheno);
opts.covartag = opts.covartag(opts.surv_covaridx);
opts.covartag2 = matlab.lang.makeValidName(opts.covartag);
if ~isempty(opts.catCovarTag)
    opts.catCovarTag2 = matlab.lang.makeValidName(opts.catCovarTag); 
end
opts.catDic = dictionary();
df = array2table(df, "VariableNames", [opts.surv_preds.preds', ...
    opts.covartag2', opts.surv_preds.Pheno(1), method_cols, "eid"]);

% if OOA, make it consistent with TTE
if opts.timeOrigin == "AOO"
    if opts.LT % swap tt0 and tt
        new_tt = df.tt0;
        df.tt0 = df.tt;
        df.tt = new_tt;
        clear new_tt
    else
        df = renamevars(df, "tt0", "tt");
    end
end

if ~opts.prevalent % remove prevalent cases
    rmidx = df.base > 0;
    if opts.LT
        rmidx = rmidx | df.tt0 <= 0 | df.tt <= 0;
    else
        rmidx = rmidx | df.tt <= 0;
    end
    df(rmidx, :) = [];
    df.base = [];
end

% define left-truncated age at enrollment age
if opts.timeOrigin == "AOO" && opts.LT
    df(df.tt0 <= 0, :) = [];
    df.tt0 = df.tt - df.tt0;
end

% remove any negative values and throw a warning
% @02JULY2024: changed <=0 to <0 (for TTE and LT)

if opts.LT
    idx = df.tt0 < 0 | df.tt < 0;
else
    idx = df.tt <= 0;
end

if any(idx)
    fprintf("\n\tWARNING: %d negative TTE/AOO were removed. If you beleive no negative values should exist check your trait!\n", sum(idx))
    df(idx, :) = [];
end

% apply landmark
idx = df.tt <= opts.landmark & df.censor ~= 0;
if any(idx)
    error("landmark option should be placed before changes to LT/timescale!")
    % fprintf("\n\tLandmark analysis: %d individuals with a tt <= %.2g yr were excluded.\n", sum(idx), opts.landmark)
    % df(idx, :) = [];
end

df.eid = [];
% summary of predictors
fnan = any(ismissing(df), 2);
event_idx = df{~fnan, opts.surv_preds.Pheno(1)} == 1;
control_idx = df{~fnan, opts.surv_preds.Pheno(1)} ~= 1;
tab = struct;
tab.SNP = opts.surv_preds.SNP;
tab.A2FREQ = mean(df{:, opts.surv_preds.preds})'./opts.denom; % freq calculated here is independent of missing values in covariates
% tab.N = repmat(sum(~fnan), numel(tab.A2FREQ), 1); % N corresponds to effectize sample size used for adjusted model
tab.AF_case = mean(df{event_idx, opts.surv_preds.preds})'./opts.denom;
tab.AF_control = mean(df{control_idx, opts.surv_preds.preds})'./opts.denom;
tab = struct2table(tab);

if isfield(opts, "interactiontag")
    opts.interactiontag2 = matlab.lang.makeValidName(opts.interactiontag);
end

% convert covariates to categorical if num unique < 3
for k = 1:numel(opts.covartag2)
    if startsWith(lower(opts.covartag2(k)), "sex") && numel(unique(df.(opts.covartag2(k)))) < 3
        df.(opts.covartag2(k)) = replace(string(df.(opts.covartag2(k))), ["0", "1"], [opts.sex_0, opts.sex_1]);
        df.(opts.covartag2(k)) = categorical(df.(opts.covartag2(k)));
        opts.catDic{opts.covartag2(k)} = opts.covartag2(k) + string(categories(df.(opts.covartag2(k))));
    elseif ~isempty(opts.catCovarTag) && any(opts.covartag2(k) == opts.catCovarTag2)
        df.(opts.covartag2(k)) = categorical(df.(opts.covartag2(k)));
        opts.catDic{opts.covartag2(k)} = opts.covartag2(k) + string(categories(df.(opts.covartag2(k))));
    end
end
% df.censor = replace(string(df.censor), ["1", "2", "3"], ["event", "death", "competing"]);
% df.censor = categorical(df.censor);

[opts.dir, out1, out2] = fileparts(string(opts.output));
opts.output = out1 + out2;
if opts.dir == ""
    opts.dir = string(pwd);
end
opts.dirr = opts.dir.replace(filesep, "/");

% df.SMS = categorical(df.SMS, [0, 1, 2], ["Control", "Stage 1", "Stage 2"]);

% find limits of the time scale variables for CIF
tcols = colnames(df);
tcols(~tcols.startsWith("tt")) = [];
axlims = table2array(varfun(@(x)[min(x), max(x)], df(:, tcols)));
axlims = [min(axlims(:)), max(axlims(:))];
if opts.timeOrigin == "AOO"
    axlims(2) = axlims(2) + 2; % 2 years to the right
else
    axlims(1) = 0;
    axlims(2) = axlims(2) + 1; % TTE
end


% transfer to R
[df, tmp_file] = table2Rdf(df, "write", true, "dir", opts.dir, "dt", "parquet");

r = string;
r(1) = "setwd('" + opts.dirr + "')";
r(2) = "suppressMessages(lapply(c" + ...
    '("dplyr", "survival", "ggsurvfit", "broom", ' + ...
    '"survminer", "tibble", "tidycmprsk","riskRegression", "pshBAR",' + ...
    '"fastcmprsk", "finalfit"), require, character.only = TRUE))';
r = [r, df'];

% r(numel(r) + 1) = "df = sample_n(df, 1000)";

r(numel(r) + 1) = "plotFlag = " + opts.plotFlag;
r(numel(r) + 1) = "checkDiagnostics = " + opts.checkDiagnostics;
r(numel(r) + 1) = "FGR = '" + opts.FGR + "' # can be 'crr' or 'fastcrr'";

r(numel(r) + 1) = "pruneStab = function(df, tag, flag=0){";
r(numel(r) + 1) = "  df = summary(df, digits = 3)";
r(numel(r) + 1) = "  df = as.data.frame(cbind(df$coef, df$conf.int))";
r(numel(r) + 1) = '  new_cols = c(tag, "SE", "P", "CIL", "CIH")';
r(numel(r) + 1) = "  if (flag == 0)";
r(numel(r) + 1) = '    old_cols = c("exp(coef)", "se(coef)", "Pr(>|z|)", "lower .95", "upper .95")';
r(numel(r) + 1) = "  else if (flag == 1)";
r(numel(r) + 1) = '    old_cols = c("exp(coef)", "se(coef)", "p-value", "2.5%", "97.5%")';
r(numel(r) + 1) = "  else";
r(numel(r) + 1) = '    old_cols = c("exp(coef)", "se(coef)", "Pr(>|z|)", "2.5%", "97.5%")';    
r(numel(r) + 1) = "  df = df %>% rename_with(~ new_cols, .cols = old_cols)";
r(numel(r) + 1) = '  df$CI = paste0("[", sprintf("%.2f", df$CIL), ",", sprintf("%.2f", df$CIH), "]")';
r(numel(r) + 1) = "  df$Term = rownames(df)";
r(numel(r) + 1) = '  df = df[, c("Term", tag, "SE", "P", "CI")]';
r(numel(r) + 1) = "  return(df)";
r(numel(r) + 1) = "}";

r(numel(r) + 2) = "# median follow-up using reverse KM: https://www.bioinfo-scrounger.com/archives/follow_up_time/";
if opts.LT
    r(numel(r) + 1) = "fit <- survfit(Surv(tt0, tt, `" + opts.surv_preds.Pheno(1) + "`==0, type='counting') ~ 1, data = df)";
else
    r(numel(r) + 1) = "fit <- survfit(Surv(tt, `" + opts.surv_preds.Pheno(1) + "`==0) ~ 1, data = df)";
end
r(numel(r) + 1) = "mf = surv_median(fit)";
r(numel(r) + 1) = "write.table(mf, 'mf." + opts.surv_preds.Pheno2(1) + ".txt', quote = F, row.names = F, col.names = T, sep = '||')";

r(numel(r) + 2) = "# KM plot: overall survival probability ";
r(numel(r) + 1) = "if (plotFlag){";

if opts.LT
    r(numel(r) + 1) = "  kmfig = survfit2(Surv(tt0, tt, `" + opts.surv_preds.Pheno(1) + "`, type='counting') ~ 1, data = df) %>% ";
else
    r(numel(r) + 1) = "  kmfig = survfit2(Surv(tt, `" + opts.surv_preds.Pheno(1) + "`) ~ 1, data = df) %>% ";
end

r(numel(r) + 1) = "    ggsurvfit() + ";
r(numel(r) + 1) = "    labs(";
r(numel(r) + 1) = '      x = "Years",';
r(numel(r) + 1) = '      y = "Overall survival probability"';
r(numel(r) + 1) = "    ) + ";
r(numel(r) + 1) = "    add_confidence_interval() +";
r(numel(r) + 1) = "    add_risktable()";
r(numel(r) + 1) = "  grDevices::cairo_pdf('km." + opts.surv_preds.Pheno2(1) + "." + opts.output + ".pdf', onefile=FALSE)";
r(numel(r) + 1) = "  print(kmfig)";
r(numel(r) + 1) = "  dev.off()}";

r(numel(r) + 2) = "# Competing risks";
r(numel(r) + 1) = "# 1- Cause-specific hazard of a given event: this represents the rate per ";
r(numel(r) + 1) = "# unit of time of the event among those not having failed from other events (below)";
r(numel(r) + 1) = "# 2- Subdistribution hazard of a given event: this represents the rate per unit";
r(numel(r) + 1) = "# of time of the event as well as the influence of competing events";

% loop over variants
r(numel(r) + 2) = "preds = colnames(df)[1:" + numel(opts.surv_preds.preds) + "]";

if opts.LT
    r(numel(r) + 1) = "dependent_dss = 'Surv(tt0, tt, `" + opts.surv_preds.Pheno(1) + "`)'";
    r(numel(r) + 1) = "dependent_crr = 'Surv(tt0, tt, censor)'";
    r(numel(r) + 1) = "dependent_crr2 = 'Surv(tt, censor2)'"; % "dependent_crr2 = 'Surv(tt0, tt, censor2)'";
else
    r(numel(r) + 1) = "dependent_dss = 'Surv(tt, `" + opts.surv_preds.Pheno(1) + "`)'";
    r(numel(r) + 1) = "dependent_crr = 'Surv(tt, censor)'";
    r(numel(r) + 1) = "dependent_crr2 = 'Surv(tt, censor2)'";
end

if opts.hascovars
    covs = createRvec(opts.covartag2, 10, true);
    covs(1) = "covs = " + covs(1);
    r = [r, covs'];
else
    r(numel(r) + 1) = "covs = NULL";
end
r(numel(r) + 1) = "res = list()";
r(numel(r) + 1) = "for (k in 1:" + numel(opts.surv_preds.preds) + "){";

% check logRank test
if opts.logRank
    r(numel(r) + 2) = "   if(data.table::uniqueN(df[[preds[k]]], na.rm = T) < 5){";
    if isfield(opts, "interactiontag")
        r(numel(r) + 1) = "     frmlr = as.formula(paste(c(dependent_dss, '~', preds[k], '+', '" + opts.interactiontag2 + "'), collapse = ''))";  
    else
        r(numel(r) + 1) = "     frmlr = as.formula(paste(c(dependent_dss, '~', preds[k]), collapse = ''))";  
    end
    r(numel(r) + 1) = "     fit = surv_fit(frmlr, data = df)";
    r(numel(r) + 1) = "     sfit = summary(fit)";
    r(numel(r) + 1) = '     median_line = ifelse(all(is.na(sfit$table[, "median"])), "none", "hv")';
    if opts.logRankFun ~= "null"
        r(numel(r) + 1) = "     fun = '" + opts.logRankFun + "'";
    else
        r(numel(r) + 1) = "     fun = NULL";
    end
    
    r(numel(r) + 1) = "     ylim = c(min(fit$surv), 1)";
    r(numel(r) + 1) = "     pxy = c(1, mean(ylim))";
    r(numel(r) + 1) = '     if (!is.null(fun) && fun == "event"){ylim = c(0, 1-min(fit$surv))}';
    
    if opts.LT
        r(numel(r) + 1) = "     gg = ggsurvplot(fit, pval = F, pval.coord = pxy,conf.int = TRUE, xlab = 'Years', break.time.by = 10,";
    else
        r(numel(r) + 1) = "     gg = ggsurvplot(fit, pval = T, ylim = ylim, pval.coord = pxy,conf.int = TRUE, xlab = 'Years', break.time.by = 10,";
    end
    r(numel(r) + 1) = "       ggtheme = theme_light(), risk.table = 'absolute', ";
    r(numel(r) + 1) = "       risk.table.y.text.col = T, risk.table.y.text = F,";
    if ~opts.LT
        r(numel(r) + 1) = '       surv.median.line = median_line, linetype = "solid",';
    else
        r(numel(r) + 1) = '       linetype = "solid",';
    end
    r(numel(r) + 1) = "       size = .5, fun = fun, censor=F, fontsize = 3, xlim =c(" + axlims(1) + "," + axlims(2) + ")";
    
    if opts.logRankFun ~= "event"
        r(numel(r) + 1) = "       ,ylim = ylim)";
    else
        r(numel(r) + 1) = "       )";
    end

    if isfield(opts, "interactiontag")
        r(numel(r) + 1) = '     gg = gg$plot +theme_bw() + theme (legend.position = "right")+';
        r(numel(r) + 1) = "       facet_grid(as.formula(paste(c(preds[k], '~', '" + opts.interactiontag2 + "'), collapse = '')))";
    end

    r(numel(r) + 1) = "     grDevices::cairo_pdf(paste(c('LogRankPlot','" + opts.surv_preds.Pheno2(1) + "' , preds[k], '" + opts.output + "', 'pdf'), collapse = '.'), onefile=FALSE)";
    r(numel(r) + 1) = "     print(gg)";
    r(numel(r) + 1) = "     dev.off()";
    
    r(numel(r) + 1) = "   }";
end

% CIF for each predictor: Aalen-Johansen estimator. For ordinary (single
% event) survival this reduces to the Kaplan-Meier estimate (i.e, the plots
% from logrank test are complement to these plots, or alternatively, if set
% logrankfun to "event", they are the same plots).
r(numel(r) + 2) = "# Cumulative incidence for competing risks (group by predictor)";
r(numel(r) + 1) = "   if (plotFlag){";
r(numel(r) + 1) = '      df$censor2 = factor(df$censor, levels = 0:3, labels = c("0", "event", "death", "competing"))';
r(numel(r) + 1) = "      df$censor2 = droplevels(df$censor2)";
if isfield(opts, "interactiontag")
    r(numel(r) + 1) = "     frmlr2 = as.formula(paste(c(dependent_crr2, '~', preds[k], '+', '" + opts.interactiontag2 + "'), collapse = ''))";  
else
    r(numel(r) + 1) = "     frmlr2 = as.formula(paste(c(dependent_crr2, '~', preds[k]), collapse = ''))";  
end
r(numel(r) + 1) = "      tryCatch({";
r(numel(r) + 1) = "        cm = tidycmprsk::cuminc(frmlr2, data = df, subset = complete.cases(df))";
r(numel(r) + 1) = '        fig = ggcuminc(cm, outcome = "event") +';
r(numel(r) + 1) = '          labs(x = "Years") +';
r(numel(r) + 1) = "          add_confidence_interval() +";
r(numel(r) + 1) = "          add_risktable(risktable_stats = 'n.risk', size=3) + add_risktable_strata_symbol(symbol = '\u2022') + scale_ggsurvfit() + scale_x_continuous(n.breaks = 6, limits =c(" + axlims(1) + "," + axlims(2) + "))";
r(numel(r) + 1) = "        grDevices::cairo_pdf(paste(c('CIF','" + opts.surv_preds.Pheno2(1) + "' , preds[k], '" + opts.output + "', 'pdf'), collapse = '.'), onefile=FALSE)";
r(numel(r) + 1) = "        print(fig)";
r(numel(r) + 1) = "        dev.off()";
r(numel(r) + 1) = "      }, error = function(e) {";
r(numel(r) + 1) = '        message("Error encountered: ", e$message)';
r(numel(r) + 1) = '        message("Skipping plot generation.")})';
r(numel(r) + 1) = "   }";

% add formulas
if isfield(opts, "interactiontag")
    r(numel(r) + 1) = "   explanatory = setdiff(covs, '" + opts.interactiontag2 + "')";
    r(numel(r) + 1) = "   explanatory = c(paste0(c(preds[k], '" + opts.interactiontag2 + "'), collapse = '*'), explanatory)";
else
    r(numel(r) + 1) = "   explanatory = c(preds[k], covs)";
end
r(numel(r) + 1) = "   spec = paste0(explanatory, collapse = '+')";
r(numel(r) + 1) = "   frm = as.formula(paste(c(dependent_dss, '~', spec), collapse = ''))";

r(numel(r) + 2) = "   # Cox regression (cause-specific)";
r(numel(r) + 1) = "   mdl0 = df %>% coxphmulti(dependent_dss, explanatory[1]) # univariate";
r(numel(r) + 1) = "   mdl = coxph(frm, data = df, subset = complete.cases(df))";
% r(numel(r) + 1) = '   tab0 = df %>% coxphmulti(dependent_dss, explanatory[1]) %>% fit2df(estimate_name = "HR", metrics = T, condense=F, digits = c(2, 2, 3))';
r(numel(r) + 1) = "   tab0 = mdl0 %>% fit2df(estimate_name = 'HR', metrics = T, condense=F)";
r(numel(r) + 1) = "   tab1 = mdl %>% fit2df(estimate_name = 'aHR', metrics = T, condense=F)";
r(numel(r) + 1) = "   mdl_metrics = rbind(as.character(tab0[[2]]), as.character(tab1[[2]]))";
r(numel(r) + 1) = "   write.table(mdl_metrics, file=paste(c('metrics','" + opts.surv_preds.Pheno2(1) + "' , preds[k], 'txt'), collapse = '.'), quote = F, row.names = F, col.names = F)";
% r(numel(r) + 1) = "   tab0 = tab0[[1]]";
% r(numel(r) + 1) = "   tab1 = tab1[[1]]";
% r(numel(r) + 1) = "   out <- df %>% summary_factorlist(dependent_dss, explanatory, column = TRUE, fit_id = TRUE, cont_cut=1) %>% ff_merge(tab0, ref_symbol = NA) %>% ff_merge(tab1, ref_symbol = NA)";
r(numel(r) + 1) = '   tab0 = pruneStab(mdl0, "HR")';
r(numel(r) + 1) = '   tab1 = pruneStab(mdl, "aHR")';
r(numel(r) + 1) = '   out = merge(tab0, tab1, by="Term", all=T)';

r(numel(r) + 2) = "   # Subdistribution hazard";
r(numel(r) + 1) = '   if (FGR == "crr"){';
r(numel(r) + 1) = "     # Fine and Gray competing risks regression";
r(numel(r) + 1) = "     mdl2 = df[complete.cases(df), ] %>% crrmulti(dependent_crr, explanatory)";
% r(numel(r) + 1) = "     tab2 = mdl2 %>% fit2df(estimate_name = 'sHR', condense=F)";
% r(numel(r) + 1) = "     out <- out %>% ff_merge(tab2, ref_symbol = NA)";
r(numel(r) + 1) = '     tab2 = pruneStab(mdl2, "sHR", flag = 1)';
r(numel(r) + 1) = '     out = merge(out, tab2, by="Term", all=T)';
r(numel(r) + 1) = '   } else if (FGR == "fastcrr") {';
r(numel(r) + 1) = "     frm2 = as.formula(paste(c('Crisk(tt, censor) ~', spec), collapse = ''))";
r(numel(r) + 1) = "     mdl2 = fastCrr(frm2, data=df)";
% r(numel(r) + 1) = "     mdl_fast = summary(mdl_fast, digits = 3)";
% r(numel(r) + 1) = "     mdl_fast = as.data.frame(cbind(mdl_fast$coef, mdl_fast$conf.int))";
% r(numel(r) + 1) = '     tcol1 = out[, c("aHR", "L95.y", "U95.y", "p.y")]';
% r(numel(r) + 1) = "     idx = complete.cases(tcol1)";
% r(numel(r) + 1) = "     tcol2 = data.frame(sHR = mdl_fast$`exp(coef)`, L95.z = mdl_fast$`2.5%`, U95.z = mdl_fast$`97.5%`, p.z = mdl_fast$`Pr(>|z|)`)";
% r(numel(r) + 1) = "     tcol1[idx, ] = tcol2";
% r(numel(r) + 1) = "     colnames(tcol1) = colnames(tcol2)";
% r(numel(r) + 1) = "     out = cbind(out, tcol1)";
r(numel(r) + 1) = '     tab2 = pruneStab(mdl2, "sHR", flag = 2)';
r(numel(r) + 1) = '     tab2$Term = tab1$Term';
r(numel(r) + 1) = '     out = merge(out, tab2, by="Term", all=T)';
r(numel(r) + 1) = "   }";

% r(numel(r) + 2) = "   out <- out %>% select(-fit_id, -index)";
r(numel(r) + 2) = "   out$SNP = preds[k]";
r(numel(r) + 1) = "   res[[k]] = out";

% more http://www.sthda.com/english/wiki/cox-model-assumptions
r(numel(r) + 1) = "   # Assessing proportional hazards";
r(numel(r) + 1) = "   # 1-Schoenfeld residuals to check the proportional hazards assumption";
r(numel(r) + 1) = "   # 2-Martingale residual to assess nonlinearity";
r(numel(r) + 1) = "   # 3-Deviance residual (symmetric transformation of the Martinguale residuals), to examine influential observations";
r(numel(r) + 1) = "   if (checkDiagnostics){";
r(numel(r) + 1) = "     cz = cox.zph(mdl)";
r(numel(r) + 1) = "     czplot = ggcoxzph(cz)";
r(numel(r) + 1) = "     grDevices::cairo_pdf(paste(c('Diagnostics','" + opts.surv_preds.Pheno2(1) + "' , preds[k], '" + opts.output + "', 'pdf'), collapse = '.'))";
r(numel(r) + 1) = "     for(m in 1:length(czplot)){print(czplot[m])}";
r(numel(r) + 1) = "     dev.off()";
r(numel(r) + 1) = "   }";
r(numel(r) + 1) = "}";

% # Time-dependent covariate (should be implemented)
% # An alternative to a landmark analysis is incorporation of a time-dependent covariate
% # This may be more appropriate than landmark analysis when:
% # the value of a covariate is changing over time
% # there is not an obvious landmark time
% # use of a landmark would lead to many exclusions

r(numel(r) + 2) = "# Cumulative incidence for competing risks (overall)";
r(numel(r) + 1) = "if (plotFlag){";
r(numel(r) + 1) = '   df$censor = factor(df$censor, levels = 0:3, labels = c("0", "event", "death", "competing"))';
r(numel(r) + 1) = "   df$censor = droplevels(df$censor)";

% cuminc doesn't support left truncation
% if opts.LT
%     r(numel(r) + 1) = "   cm = tidycmprsk::cuminc(Surv(tt0, tt, censor) ~ 1, data = df, subset = complete.cases(df))";
% else
    % r(numel(r) + 1) = "   cm = tidycmprsk::cuminc(Surv(tt, censor) ~ 1, data = df, subset = complete.cases(df))";
% end
% 
r(numel(r) + 1) = "   tryCatch({";
r(numel(r) + 1) = "     cm = tidycmprsk::cuminc(Surv(tt, censor) ~ 1, data = df, subset = complete.cases(df))";
r(numel(r) + 1) = '     fig = ggcuminc(cm, outcome = setdiff(levels(df$censor), "0")) +';
r(numel(r) + 1) = '       labs(x = "Years") +';
r(numel(r) + 1) = "       add_confidence_interval() +";
r(numel(r) + 1) = "       add_risktable()";
r(numel(r) + 1) = "     grDevices::cairo_pdf('CIF." + opts.surv_preds.Pheno2(1) + "." + opts.output + ".pdf', onefile=FALSE)";
r(numel(r) + 1) = "     print(fig)";
r(numel(r) + 1) = "     dev.off()";
r(numel(r) + 1) = "   }, error = function(e) {";
r(numel(r) + 1) = '     message("Error encountered: ", e$message)';
r(numel(r) + 1) = '     message("Skipping plot generation.")})';
r(numel(r) + 1) = " }";

res_name = getRandomName("res", 6);
r(numel(r) + 2) = "res = bind_rows(res)";
r(numel(r) + 1) = "data.table::fwrite(res, '" + res_name + ".txt', sep = ';', row.names = F, quote = F, col.names = T)";

r = addRfont(r);
MATLAB2Rconnector(fullfile(opts.dir, "survHelper2" + res_name + ".r"), code=r);
delete(tmp_file)

% read results
file = fullfile(opts.dir, res_name + ".txt");
res = readtable(file, "Delimiter", ";", ...
    "NumHeaderLines", 0, "TextType", "string", "VariableNamingRule", "preserve", ...
    "ReadVariableNames", true);
delete(file)

old_cols = ["SE.x", "P.x", "CI.x", "SE.y", "P.y", "CI.y"];
new_cols = [["SE", "P", "CI"] + ".uni", ["SE", "P", "CI"] + ".multi"];
if any(colnames(res) == "sHR")
    old_cols = [old_cols, "SE", "P", "CI"];
    new_cols = [new_cols, ["SE.sub", "P.sub", "CI.sub"]];
end
res = renamevars(res, old_cols, new_cols);

% remove nuisance covars
idx = cellfun(@(x)any(regexpi(x, "cov\d+")), res.Term);
res(idx, :) = [];

% read median follow-up (reverese KM)
file = fullfile(opts.dir, "mf." + opts.surv_preds.Pheno2(1) + ".txt");
mf = readtable(file, "Delimiter", "||", ...
    "NumHeaderLines", 0, "TextType", "string", "VariableNamingRule", "preserve", ...
    "ReadVariableNames", true);
mf = compose("%.3f (%.3f-%.3f)", mf.median, mf.lower, mf.upper);
delete(file)

model_metrics = struct;
% read model metrics
for k = 1:numel(opts.surv_preds.preds)
    file = fullfile(opts.dir, "metrics." + opts.surv_preds.Pheno2(1) + ...
        "." + opts.surv_preds.preds(k) + ".txt"); 
    tmp = readmatrix(file, "OutputType", "string", "Delimiter", ",");
    model_metrics.N{k, 1} = double(tmp(:, 2).extractAfter("Number in model = "))';
    model_metrics.N_events{k, 1} = double(tmp(:, 4).extractAfter("Number of events = "))';
    model_metrics.concordance{k, 1} = tmp(:, 5).extractAfter("Concordance = ")'; 
    lrt = double(tmp(:, 7).extractBetween("Likelihood ratio test = ", "(df"));
    lrt_df = double(tmp(:, 7).extractAfter("(df = "));
    lrt_p = double(tmp(:, 8).extractBetween("p = ", ")"));
    model_metrics.LRT{k, 1} = compose("LRT = %.2f (df = %d, p = %.2g)", lrt, lrt_df, lrt_p)';
    delete(file)
end

model_metrics = struct2table(model_metrics);
model_metrics = convertvars(model_metrics, 1:width(model_metrics), @(x)vertcat(x{:}));
model_metrics = splitvars(model_metrics);
newCols = ["N", "N_events", "concordance", "LRT"] + [" (uni)"; " (multi)"];
model_metrics = renamevars(model_metrics, colnames(model_metrics), newCols(:));
model_metrics.("Median follow-up")(:) = mf;
cols = ["Pheno", "SNP", "CHR", "POS", "A1", "A2"];
for k = 1:numel(cols), model_metrics.(cols(k)) = opts.surv_preds.(cols(k)); end
model_metrics = movevars(model_metrics, cols, Before=1);

% prune and merge tables
res.Pheno(:) = opts.surv_preds.Pheno(1);

% variable names should be back to original names
% if ~isempty(opts.catCovarTag)
%     catDict = opts.catDic.entries;
%     [~, idx] = ismember(catDict.Key, opts.covartag2);
%     catDict.Key2 = opts.covartag(idx);
% 
% end

mapTab = struct2table(opts.surv_preds);
for k = 1:height(mapTab)
    idx = ismember(res.Term, mapTab.preds(k)) | startsWith(res.Term, mapTab.preds(k) + ":");
    res.Term(idx) = replace(res.Term(idx), mapTab.preds(k), mapTab.SNP(k));
    res.SNP(res.SNP == mapTab.preds(k)) = mapTab.SNP(k);
end

% append SNP freq data to model_metric table
appCols = setdiff(colnames(tab), "SNP");
[~, idx] = ismember(model_metrics.SNP, tab.SNP);
model_metrics(:, appCols) = tab(idx, appCols);

% merge all tables
res = outerjoin(res, model_metrics, "Keys", ["SNP", "Pheno"], "MergeKeys", true, "Type", "full");
res = movevars(res, ["Pheno", "SNP", "CHR", "POS", "A1", "A2", "A2FREQ"], Before=1);
res = renamevars(res, ["AF_case", "AF_control"], ["AF.case", "AF.control"]);
res = movevars(res, ["Median follow-up", "AF.case", "AF.control"], Before="concordance (uni)");

end % END