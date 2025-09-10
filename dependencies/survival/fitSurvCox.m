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
    
    % competing-risk model choice
    opts.FGR {mustBeTextScalar, mustBeMember(opts.FGR, ["none", "crr", "fastcrr"])} = "none" % Fine-Gray for competing risks: note that "fastcrr" doesn't support LT
    opts.nBoot (1,1) double = 5  % bootstrap reps for CIF ribbons
    
    % filtering and choice of time-scale
    opts.prevalent (1,1) logical = false % exclude prevalent cases
    opts.timeOrigin {mustBeMember(opts.timeOrigin, ["TTE", "AOO"])} = "TTE"
    opts.LT (1,1) logical = false % no left-truncation
    opts.landmark (1,1) double = 0 % landmark analysis, individuals with a non-zero censoring status (including competing risks) with a tt <= landmark (yr) are excluded. 

    % path for helper R function fitSurvCox.R
    opts.rfun {mustBeFile} = fullfile(fileparts(which("fitSurvCox.m")), "fitSurvCox.R")

end

if numel(fieldnames(opts2))
    fis = string(fieldnames(opts2));
    for k = 1:numel(fis)
        opts.(fis(k)) = opts2.(fis(k));
        opts2 = rmfield(opts2, fis(k));
    end
end

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
method_cols = ["censor", "tt", "tt0"];

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

if opts.timeOrigin == "TTE"
    %@27JUNE2025: see examples above
    idx = df.tt < 0; % df.tt0 < 0 | df.tt < 0;
else
    idx = df.tt <= 0;
end

if any(idx)
    fprintf("\n\tWARNING: %d negative TTE/AOO were removed. If you beleive no negative values should exist check your trait!\n", sum(idx))
    df(idx, :) = [];
end

% % apply landmark
% if opts.timeOrigin ~= "TTE"
%     idx = df.tt <= opts.landmark & df.censor ~= 0;
%     if any(idx)
%         error("landmark option should be placed before changes to LT/timescale!")
%         % fprintf("\n\tLandmark analysis: %d individuals with a tt <= %.2g yr were excluded.\n", sum(idx), opts.landmark)
%         % df(idx, :) = [];
%     end
% end

% time-varying covariate
if isfield(opts, "tv_covar")
    tv_covar = opts.tv_covar;

    % verify required fields
    fnames = string(fieldnames(tv_covar));
    reqnames = [ "eid", "start", "stop", "value", "name"];
    assert(all(ismember(reqnames, fnames)), ...
        "fitSurvCox:Time-varying covariates must have the columns:" + ...
        "eid, start, stop, value, and name, using the same time scale " + ...
        "as the outcome.")

    tv_covar = struct2table(rmfield(tv_covar, setdiff(fnames, reqnames)));
    assert(all(tv_covar.start <= tv_covar.stop), "fitSurvCox:tv_covar:start cannot be > stop time!")
    [~, keep_idx] = ismember(df.eid, tv_covar.eid); 
    keep_idx(keep_idx < 1) = [];
    tv_covar = tv_covar(keep_idx, :);

    if isempty(tv_covar)
        fprintf("\tfitSurvCox:after merging with df no subject remained in tv_covar --> skipped\n")
        df.eid = [];
    end
    
else
    tv_covar = table;
    df.eid = [];
end

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

% transfer to R ===========================================================
% 1- prepare R function flags

ropts = struct;
ropts.plotCIF = false;
ropts.overallKM = false;
if opts.plotFlag
    ropts.plotCIF = true; ropts.overallKM = true;
end

ropts.plotCIF = false;

ropts.diag = false;
if opts.checkDiagnostics, ropts.diag = true; end

% 2- write to parquet (df and/or tv covar) --------------------------------
[df, tmp_file] = table2Rdf(df, ...
    write=true, ...
    dir=opts.dir, ...
    dt="parquet", ...
    name="df");

if ~isempty(tv_covar)
    [~, tv_file] = fileparts(tmp_file);
    tv_file = tv_file + "_tv";
    [df_tv, tv_file] = table2Rdf(tv_covar, ...
        write=true, ...
        dir=opts.dir, ...
        dt="parquet", ...
        fname=tv_file, ...
        name="tv_covar");
    df = [df; df_tv];
    df = unique(df, "stable");
else
    tv_file = [];
end

% 3- prepare other inputs -------------------------------------------------
r = string;
r(1) = "setwd('" + opts.dirr + "')";
r(2) = "source('" + opts.rfun.replace("\", "/") + "')";
r = [r, df'];

if opts.hascovars
    ropts.covars = opts.covartag2;
else
    ropts.covars = nan;
end

% r(numel(r) + 1) = 'outcomeVar = "' + opts.surv_preds.Pheno(1) + '"';
ropts.outcomeVar = opts.surv_preds.Pheno(1); 
ropts.competeVar = "censor";
ropts.preds = opts.surv_preds.preds;

if isfield(opts, "interactiontag")
    ropts.interactionVar = opts.interactiontag2;
else
    ropts.interactionVar = nan;
end

if isempty(tv_covar)
    ropts.tv_df = nan;
else
    ropts.tv_df = {'tv_covar'};
end

sfis = ["timeOrigin", "LT", "prevalent", "landmark", "FGR", "logRank", ...
    "logRankFun", "parallel", "timeDepVars", "nonlinVars", "splineDF", ...
    "timeDepFun"];
for k = 1:numel(sfis)
    if isfield(opts, sfis(k))
        checkVal = opts.(sfis(k));
        if ismember(class(checkVal), ["cell", "char", "string"]) && ...
                (all(ismissing(checkVal)) || all(checkVal == ""))
            continue
        end
        ropts.(sfis(k)) = opts.(sfis(k));

        if sfis(k) == "timeDepFun"
            ropts.timeDepFun = cellstr(ropts.timeDepFun);
        end
    end
end

tmp_rstr = struct2rfun(ropts,par=false, pretty=true);
tmp_rstr(1) = "res = fitSurvCox(df," + tmp_rstr(1);
tmp_rstr(end) = tmp_rstr(end) + ")";
r = [r, tmp_rstr'];

res_name = getRandomName("res", 6);
r(numel(r) + 1) = "data.table::fwrite(res$results, '" + res_name + ".txt', sep = ';', row.names = F, quote = F, col.names = T)";
r(numel(r) + 1) = "write.table(res$followUp, 'mf." + res_name + ".txt', quote = F, row.names = F, col.names = T, sep = '||')";
r(numel(r) + 1) = "write.table(res$metrics, 'metrics." + res_name + ".txt', quote = F, row.names = F, col.names = T, sep = '||')";

MATLAB2Rconnector(fullfile(opts.dir, "survHelper2" + res_name + ".r"), code=r);
delete(tmp_file)
if ~isempty(tv_covar), delete(tv_file); end

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
file = fullfile(opts.dir, "mf." + res_name + ".txt");
mf = readtable(file, "Delimiter", "||", ...
    "NumHeaderLines", 0, "TextType", "string", "VariableNamingRule", "preserve", ...
    "ReadVariableNames", true);
mf = compose("%.3f (%.3f-%.3f)", mf.median, mf.lower, mf.upper);
delete(file)

% read model metrics
file = fullfile(opts.dir, "metrics." + res_name + ".txt");
model_metrics = readtable(file, "Delimiter", "||", ...
    "NumHeaderLines", 0, "TextType", "string", "VariableNamingRule", "preserve", ...
    "ReadVariableNames", true);
delete(file)

model_metrics = unstack(model_metrics, ["N", "N_events", "Concordance", "LRT"], "Model");
old_cols = colnames(model_metrics);
new_cols = old_cols.replace("_multi" + textBoundary("end"), " (multi)");
new_cols = new_cols.replace("_uni" + textBoundary("end"), " (uni)");
model_metrics = renamevars(model_metrics, old_cols, new_cols);

model_metrics.("Median follow-up")(:) = mf;
cols = ["Pheno", "SNP", "CHR", "POS", "A1", "A2"];
for k = 1:numel(cols), model_metrics.(cols(k)) = opts.surv_preds.(cols(k)); end
model_metrics = movevars(model_metrics, cols, Before=1);

% prune and merge tables
res.Pheno(:) = opts.surv_preds.Pheno(1);
res = renamevars(res, "Predictor", "SNP");

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
res = movevars(res, ["Median follow-up", "AF.case", "AF.control"], Before="Concordance (multi)");

end % END