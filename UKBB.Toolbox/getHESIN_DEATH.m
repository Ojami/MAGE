function UKB_STRUCT_ALL = getHESIN_DEATH(interms, opts)
% a subfunction of UKBtrait_toolbox for getting hospitalized (HESIN) or
% death registry codes (ICD9 and ICD10) from directly downloaded data from
% UKB. This data is the latest update for the above-mentioned information
% released during COVID-19 pandemic by UKB, and therefore was not available
% in form of typical databaskets.
%
% Oveis Jamialahmadi, GU, Oct. 2020.
% @04/07/2022: major modifications were made.
% @31/10/2022: instead of returning only 1 code per individual, now merges
%              all codes (, separated).
% @22/02/2024: 'surv' flag was added for survival analysis (if true,
%               returns also date of events)
% @22/02/204: 'ex' option was added for codes for exclusion criteria joined
%              with ",". Note that the 'ex' codes are not removed from
%              returned pheno, so be cautious when using this field (see
%              phenoParser).
% 
% @24/02/2024: 'dlf' option was added and accepts a struct with 'eid' and
%              'rawUKB' fields, containing the id and date each individual
%              lost the follow-up (optional).
% @28MAY2024: 'mortality' flag was added (default is false). If true, only
%             returns death data.

arguments
    interms 
    opts.basketdir {mustBeFolder} = fullfile(fileparts(which('phenoParser.m')), 'UKBFileParser_HESIN_DEATH')
    opts.meaning (1,1) logical = true % if true, additionalDict field is present (deprecated)
    opts.verbose (1,1) logical = true
    opts.surv (1,1) logical = false 
    opts.ex {mustBeTextScalar} = ""
    opts.dlf {mustBeA(opts.dlf, "struct")} % Date lost to follow-up, must be struct (optional)
    opts.mortality (1,1) logical = false
end

UKB_STRUCT_ALL = struct; 
df = interms.fieldID{1}; % input data-fields
keyterms = split(interms.secondaryTerms{1}, ',');
opts.keyterms = keyterms; % original keyterms

if opts.ex ~= "" && ~ismissing(opts.ex)
    opts.exflag = true;
    opts.ex = unique(split(opts.ex, ","));  % exclusion criteria

    % @14MAY2024: merge 'ex' with keyterms to speed-up the matching
    keyterms = union(keyterms, opts.ex);
    keyterms(keyterms == "" | ismissing(keyterms)) = [];
else
    opts.exflag = false;
end

if opts.meaning
    [idx1, idx2] = ismember(keyterms, interms.additionalDict(:, 1));
    if ~all(idx1) % some terms cannot be found
        fprintf('some terms cannot be found: %s\n', join(keyterms(~idx1), ';'))
    end
    keyterms = keyterms(idx1);
    keymeanings = interms.additionalDict(idx2, 2);
end

% load censoring date and take the min between HESIN/DEATH censoring dates.
% These date are fetched directly from UKBB.
% @26FEB2023: each registry should be handled separately: death and HESIN.
if opts.surv
    cen_date = load(fullfile(opts.basketdir, 'censoring_info.mat')).cd;
    % fis = string(fieldnames(cen_date.d));
    % for k = 1:numel(fis) % for England, Scottland, and Wales
    %     cen_date.(fis(k)) = min(cen_date.d.(fis(k)), cen_date.h.(fis(k)));
    % end
end

% step1: check death registry ---------------------------------------------
death_cause_df = find(ismember({'40001', '40002'}, df));
[death, death_ex, dc] = deal(struct);
if ~isempty(death_cause_df) % was cause of death requested?
    if opts.verbose
        tic
        fprintf('parsing death data...')
    end
    
    dc = struct2table(load(fullfile(opts.basketdir, 'death_cause.mat')));
    if all(keyterms == "") % read all causes of death 
        f_cause = true(height(dc), 1);
    else
        f_cause = ismember(dc.cause_icd10, keyterms);
    end

    % if ~opts.exflag % exclusion criteria
    %     f_cause_ex = false(numel(dc.cause_icd10), 1);
    % else
    %     f_cause_ex = ismember(dc.cause_icd10, opts.ex);
    % end

    %@22FEB2024: 'level' is no longer used and both primary/secondary
    %causes are being used
    if any(colnames(dc) == "level")
        f_level = ismember(dc.level, death_cause_df); % primary or secondary is queried 
    else
        f_level = true(numel(dc.eid), 1);
    end

    % subset of dc with 'ex' and 'keyterms' codes
    dc2 = dc(f_cause & f_level, :);
    
    if opts.mortality
        f_cause = true(height(dc2), 1);
    else
        f_cause = ismember(dc2.cause_icd10, opts.keyterms); % original keyterms
    end

    death.eid = dc2.eid(f_cause);
    death.rawUKB = dc2.cause_icd10(f_cause);
    death.source = dc2.source(f_cause);

    if opts.exflag
        f_cause_ex = ismember(dc2.cause_icd10, opts.ex); 
        death_ex.eid = dc2.eid(f_cause_ex);
        death_ex.rawUKB = dc2.cause_icd10(f_cause_ex);
    end

    if opts.surv % fetch also date of death and date of attending
        death.date = dc2.date_of_death(f_cause);
        death.regdate = dc2.regdate(f_cause);
        death.dob = dc2.dob(f_cause);

        if opts.exflag
            death_ex.date = dc2.date_of_death(f_cause_ex);
            death_ex.regdate = dc2.regdate(f_cause_ex);
            death_ex.dob = dc2.dob(f_cause_ex);
            death_ex.source = dc2.source(f_cause_ex);
        end
    end

    clear f_cause f_level dc2
    if ~opts.surv, clear("dc"); end

    % keep all codes joined by ,
    if ~isempty(death.eid)
        death = struct2table(death);

        if opts.surv % keep the code with min date of death
            death = groupfilter(death, "eid", @(x)x == min(x), "date");
            death_tmp = death;
            % death = groupsummary(death, "eid", @(x)join(unique(string(x)), ','), ["rawUKB", "date", "regdate", "source", "dob"]);
            death = groupsummary(death, "eid", @(x)join(unique(x), ','), "rawUKB");
            [~, idx] = ismember(death.eid, death_tmp.eid);
            death(:, ["date", "regdate", "dob", "source"]) = death_tmp(idx, ["date", "regdate", "dob", "source"]);
            clear death_tmp

        else
            death = groupsummary(death, "eid", @(x)join(unique(x), ','), "rawUKB"); 
        end

        death.GroupCount = [];
        death.Properties.VariableNames = regexprep(death.Properties.VariableNames, "^fun1_", "");
        if ~opts.surv, death = table2struct(death, "ToScalar", true); end
    end

    if opts.exflag && ~isempty(death_ex.eid)
        death_ex = struct2table(death_ex);

        if opts.surv % keep the code with min date of death
            death_ex = groupfilter(death_ex, "eid", @(x)x == min(x), "date");
            death_ex_tmp = death_ex;
            death_ex = groupsummary(death_ex, "eid", @(x)join(unique(x), ','), "rawUKB");
            [~, idx] = ismember(death_ex.eid, death_ex_tmp.eid);
            death_ex(:, ["date", "regdate", "dob", "source"]) = death_ex_tmp(idx, ["date", "regdate", "dob", "source"]);
            clear death_ex_tmp
            % death_ex = groupsummary(death_ex, "eid", @(x)join(unique(string(x)), ','), ["rawUKB", "date", "regdate", "source", "dob"]);
            
        else
            death_ex = groupsummary(death_ex, "eid", @(x)join(unique(x), ','), "rawUKB");
        end
        death_ex.GroupCount = [];
        death_ex.Properties.VariableNames = regexprep(death_ex.Properties.VariableNames, "^fun1_", "");
        if ~opts.surv, death_ex = table2struct(death_ex, "ToScalar", true); end

    end
    
%     % keep only one icd code per each eid (first record)
%     [death.eid, uniqueIdx] = unique(death.eid, 'stable');
%     death.rawUKB = death.rawUKB(uniqueIdx);
%     clear uniqueIdx

    if ~opts.surv
        death.source = [];
    end

    if opts.verbose
        fprintf('\b\b Done.\n')
        toc
    end
end

% step2: check HESIN diagnosis data ---------------------------------------
hesin_df = find(ismember({'41270', '41202', '41204',...
    '41271', '41203', '41205'}, df));
[hesin, hesin_ex] = deal(struct);
if ~opts.mortality && ~isempty(hesin_df)
    if opts.verbose
        tic
        fprintf('parsing HESIN data...')
    end
    
    icd10flag = false; icd9flag = false;
    if any(ismember(hesin_df, 1:3)) % ICD10
        icd10flag = true;
    end
    if any(ismember(hesin_df, 4:6))
        icd9flag = true;
    end
    
    % @22FEB2024: level is no longer used
    if icd9flag && icd10flag
        vars = {'eid', 'diag_icd9', 'diag_icd10'};
    elseif icd9flag
        vars = {'eid', 'diag_icd9'};
    else
        vars = {'eid', 'diag_icd10'};
    end

    if opts.surv
        vars = union(vars, {'epistart', 'regdate', 'source', 'dob'});
    end
    
    hd = load(fullfile(opts.basketdir, 'hesin_diag.mat'), vars{:});
    
    % @22FEB2024: level is no longer used
    % if isfield(hd, "level")
    %     if all(hesin_df == 2) || all(hesin_df == 5)% 41202 or 41203: level 1
    %         f_level = hd.level == 1;
    %     elseif all(hesin_df == 3) || all(hesin_df == 6)% 41204 or 41205: level 2
    %         f_level = hd.level == 2;
    %     else
    %         f_level = true(numel(hd.level), 1);
    %     end
    %     hd = rmfield(hd, 'level');
    % else
    %     f_level = true(numel(hd.eid), 1);
    % end
    
    % hd.eid(f_level) = [];
    if icd9flag
        % hd.diag_icd9(f_level) = [];
        if all(keyterms == "")
            f_cause9 = false(numel(hd.eid), 1);
        else
            f_cause9 = ismember(hd.diag_icd9, keyterms);
        end
    else
        f_cause9 = false(numel(hd.eid), 1);
    end
    
    if icd10flag
        % hd.diag_icd10(f_level) = [];
        if all(keyterms == "")
            f_cause10 = false(numel(hd.eid), 1);
        else
            f_cause10 = ismember(hd.diag_icd10, keyterms);
        end
    else
        f_cause10 = false(numel(hd.eid), 1);
    end
    
    % @22FEB2024: level is no longer used
    % idx = f_cause9 & f_cause10 & f_level;
    idx = f_cause9 | f_cause10;

    hd = struct2table(hd);
    hd(~idx, :) = [];

    if icd9flag
        if all(opts.keyterms == "") && ~opts.surv
            f_cause9 = true(numel(hd.eid), 1);
        else
            f_cause9 = ismember(hd.diag_icd9, opts.keyterms);
        end
    else
        f_cause9 = true(numel(hd.eid), 1);
    end
    
    if icd10flag
        if all(opts.keyterms == "") && ~opts.surv
            f_cause10 = true(numel(hd.eid), 1);
        else
            f_cause10 = ismember(hd.diag_icd10, opts.keyterms);
        end
    else
        f_cause10 = true(numel(hd.eid), 1);
    end
    idx = f_cause9 & f_cause10;

    if opts.exflag % exclusion criteria
        if icd9flag
            f_cause9_ex = ismember(hd.diag_icd9, opts.ex);
        else
            f_cause9_ex = true(numel(hd.eid), 1);
        end
        
        if icd10flag
            f_cause10_ex = ismember(hd.diag_icd10, opts.ex);
        else
            f_cause10_ex = true(numel(hd.eid), 1);
        end

        idx_ex = f_cause9_ex & f_cause10_ex;
    end


    hesin.eid = hd.eid(idx);
    hesin.rawUKB = hd.diag_icd10(idx);

    if opts.exflag
        hesin_ex.eid = hd.eid(idx_ex);
        hesin_ex.rawUKB = hd.diag_icd10(idx_ex);
    end

    if opts.surv % fetch also date of death and date of attending
        hesin.date = hd.epistart(idx);
        hesin.regdate = hd.regdate(idx);
        hesin.dob = hd.dob(idx);
        hesin.source = hd.source(idx);

        if opts.exflag
            hesin_ex.date = hd.epistart(idx_ex);
            hesin_ex.regdate = hd.regdate(idx_ex);
            hesin_ex.dob = hd.dob(idx_ex);
            hesin_ex.source = hd.source(idx_ex);
        end
    end

    clear hd f_cause f_level f_cause10 f_cause10 idx idx_ex
    
    % keep all codes joined by ,
    if ~isempty(hesin.eid)
        hesin = struct2table(hesin);

        if opts.surv % keep the code with min date of death
            hesin = groupfilter(hesin, "eid", @(x)x == min(x), "date");
            hesin_tmp = hesin;
            hesin = groupsummary(hesin, "eid", @(x)join(unique(x), ','), "rawUKB");
            [~, idx] = ismember(hesin.eid, hesin_tmp.eid);
            hesin(:, ["date", "regdate", "dob", "source"]) = hesin_tmp(idx, ["date", "regdate", "dob", "source"]);
            clear hesin_tmp
            
        else
            hesin = groupsummary(hesin, "eid", @(x)join(unique(x), ','), "rawUKB");
            
        end
        hesin.GroupCount = [];
        hesin.Properties.VariableNames = regexprep(hesin.Properties.VariableNames, "^fun1_", "");
        if ~opts.surv, hesin = table2struct(hesin, "ToScalar", true); end
    end

    if opts.exflag && ~isempty(hesin_ex.eid)
        hesin_ex = struct2table(hesin_ex);

        if opts.surv % keep the code with min date of death
            hesin_ex = groupfilter(hesin_ex, "eid", @(x)x == min(x), "date");
            hesin_ex_tmp = hesin_ex;
            hesin_ex = groupsummary(hesin_ex, "eid", @(x)join(unique(x), ','), "rawUKB");
            [~, idx] = ismember(hesin_ex.eid, hesin_ex_tmp.eid);
            hesin_ex(:, ["date", "regdate", "dob", "source"]) = hesin_ex_tmp(idx, ["date", "regdate", "dob", "source"]);
            clear hesin_ex_tmp
            
        else
            hesin_ex = groupsummary(hesin_ex, "eid", @(x)join(unique(x), ','), "rawUKB");
            
        end
        hesin_ex.GroupCount = [];
        hesin_ex.Properties.VariableNames = regexprep(hesin_ex.Properties.VariableNames, "^fun1_", "");
        if ~opts.surv, hesin_ex = table2struct(hesin_ex, "ToScalar", true); end
    end
    
%     % keep only one icd code per each eid (first record)
%     [hesin.eid, uniqueIdx] = unique(hesin.eid, 'stable');
%     hesin.rawUKB = hesin.rawUKB(uniqueIdx);

    if opts.verbose
        fprintf('\b\b Done.\n')
        toc
    end
end

if opts.surv 

    % merge all together but keep also the dates for later analysis
    % time to death: 
    %   1- those who are dead: date of death - attending date
    %   2- those who lost follow-up: date lost follow up - attending date
    %   3- for the rest: censoring date - attending date
    % censoring:
    %   1- those developed outcome: code 1
    %   2- those who died before outcome: code 2
    %   3- those with competing risk: code 3, either
    %       1)before (death or diagnosis) outcome 
    %       2)no overlap with outcome (still should be censored)
    %   4- those after censoring date/lost follow-up: code 0

    if ~isstruct(death) || (opts.exflag && ~isstruct(death_ex))
        % we don't need to check if anyone left the study before death,
        % otherwise the info woudn't be appeared in the death registry.
        % Moreover, if anything, it increases the power by extending the
        % censoring date (date of death)
        % @28MAY2024: we now consider censoring for those who died after
        % leaving the study: to be consistent with other traits (e.g. HESIN
        % vs overall mortality)
        
        if isstruct(death)
            death = death_ex;
            death.date(:) = datetime("today", Format="dd/MM/uuuu"); % to be censored, since there are no cases
            ex_code_force = 2;  % only exclusiong criteria
        else
            ex_code_force = 1;
        end
        % death = convertvars(death, ["date", "regdate"], @(x)datetime(x, Format="dd/MM/uuuu"));

        % check exclusion criteria and merge it with death table
        % exclusion criteria (overall, more fine-tuned definition is
        % 'cencor')
        death.ex(:) = 0;
        if ~isstruct(death_ex)
            % death_ex = convertvars(death_ex, ["date", "regdate"], @(x)datetime(x, Format="dd/MM/uuuu"));
            [f1, f2] = ismember(death.eid, death_ex.eid); f2(f2<1) = [];
            death.ex(f1) = 1; % having both exclusion and main codes
            death.excode(:) = "";
            death.excode(f1) = death_ex.rawUKB(f2);
            
            % for censoring: check if exclusion (or competing) codes
            % happened before the death because of outcome
            idx_ex_before_outcome = find(f1);
            f12 = death_ex.date(f2) < death.date(f1);
            idx_ex_before_outcome = idx_ex_before_outcome(f12);
            death.ex(idx_ex_before_outcome) = -1; % will be changed to 3 later
            death.date(idx_ex_before_outcome) = death_ex.date(f2(f12)); % replace date with date of exclusion/competing risk

            death_ex(f2, :) = []; % overlapping
            death_ex.ex(:) = 2; % purely exclusion criteria (or competing risk)
            death_ex.excode = death_ex.rawUKB;
            death = [death; death_ex];
        else
            death.excode(:) = "";
        end
        
        % check source of data: England, Scottland and Wales. Also, check
        % if censoring date is < date of death (if so, censor)
        % death.source = double(death.source);
        idx0 = death.source == 0 & (cen_date.d.HES < death.date); %E/W
        idx1 = death.source == 1 & (cen_date.d.SMR < death.date); % Scottland 
        death.censor(:) = 1;
        death.censor(idx0 | idx1) = 0; % censored because died after the last day of censoring date
        death.date(idx0) = cen_date.d.HES; % replace censoring date for time-to-even calculation
        death.date(idx1) = cen_date.d.SMR; % same as above
        death.source = [];

        % those who died of exclusion criteria (e.g. competing risk)
        % before the date of death because of the outcome of interest
        % should be censored (or properly handled later using Fine-Gray
        % analysis)
        idx_ex_before_outcome = death.ex < 0;
        death.censor(idx_ex_before_outcome) = 3;
        death.ex(idx_ex_before_outcome) = ex_code_force;

        % remaining dead: died due to competing risk and free of event,
        % still should be censored
        death.censor(death.ex == 2 & death.censor ~= 0) = 3; 
        
    end
    
    % now analyse the HESIN data
    if ~isstruct(hesin) || (opts.exflag && ~isstruct(hesin_ex))

        if isstruct(hesin)
            hesin = hesin_ex;
            hesin.date(:) = datetime("today", Format="dd/MM/uuuu"); % to be censored, since there are no cases
            ex_code_force = 2; % only exclusion criteria
        else
            ex_code_force = 1;
        end

        % hesin = convertvars(hesin, ["date", "regdate"], @(x)datetime(x, Format="dd/MM/uuuu"));

        % check exclusion criteria and merge it with hesin table
        % exclusion criteria (overall, more fine-tuned definition is
        % 'cencor')
        hesin.ex(:) = 0;
        if ~isstruct(hesin_ex)
            % hesin_ex = convertvars(hesin_ex, ["date", "regdate"], @(x)datetime(x, Format="dd/MM/uuuu"));

            [f1, f2] = ismember(hesin.eid, hesin_ex.eid); f2(f2<1) = [];
            hesin.ex(f1) = 1; % having both exclusion and main codes
            hesin.excode(:) = "";
            hesin.excode(f1) = hesin_ex.rawUKB(f2);
            
            % for censoring: check if exclusion (or competing) codes
            % happened before the outcome diagnosis
            idx_ex_before_outcome = find(f1);
            f12 = hesin_ex.date(f2) < hesin.date(f1);
            idx_ex_before_outcome = idx_ex_before_outcome(f12);
            hesin.ex(idx_ex_before_outcome) = -1; % will be changed to 3 later
            hesin.date(idx_ex_before_outcome) = hesin_ex.date(f2(f12)); % replace date with date of exclusion/competing risk

            hesin_ex(f2, :) = []; % overlapping
            hesin_ex.ex(:) = 2; % purely exclusion criteria (or competing risk)
            hesin_ex.excode = hesin_ex.rawUKB;
            hesin = [hesin; hesin_ex];
        else
            hesin.excode(:) = "";
        end
        
        % check source of data: England, Scottland and Wales. Also, check
        % if censoring date is < date of diagnosis (if so, censor)
        % hesin.source = double(hesin.source);
        idx0 = hesin.source == 0 & (cen_date.h.HES < hesin.date); % England
        idx1 = hesin.source == 1 & (cen_date.h.SMR < hesin.date); % Scottland 
        idx2 = hesin.source == 2 & (cen_date.h.PEDW < hesin.date); % Wales 
        hesin.censor(:) = 1;
        hesin.censor(idx0 | idx1 | idx2) = 0; % censored because diagnoses after the last day of censoring date
        hesin.date(idx0) = cen_date.h.HES; % replace censoring date for time-to-even calculation
        hesin.date(idx1) = cen_date.h.SMR; % same as above
        hesin.date(idx2) = cen_date.h.PEDW; % same as above
        hesin.source = [];

        % those who diagnosed with exclusion criteria (e.g. competing risk)
        % before the date of diagnosis with the outcome of interest
        % should be censored (or properly handled later using Fine-Gray
        % analysis)
        idx_ex_before_outcome = hesin.ex < 0;
        hesin.censor(idx_ex_before_outcome) = 3;
        hesin.ex(idx_ex_before_outcome) = ex_code_force;

        % remaining alive: diagnoses with competing risk and free of event,
        % still should be censored
        hesin.censor(hesin.ex == 2 & hesin.censor ~= 0) = 3; 
        
        % check if they left the study before the diagnosis (either ex or main outcome): code 0
        if isfield(opts, "dlf")
            [f1, f2] = ismember(hesin.eid, opts.dlf.eid); f2(f2<1) = [];
            if any(f1)
                idx_lf_before_outcome = find(f1); % lost follow-up before outcome
                idx = opts.dlf.rawUKB(f2) < hesin.date(f1);
                idx2 = idx_lf_before_outcome(idx);
                hesin.censor(idx2) = 0;
                hesin.date(idx2) = opts.dlf.rawUKB(f2(idx)); % use the date lost-follow up for these individuals
            end
        end

    end

    % final step: merge HESIN and DEATH
    %   1- those who are dead: date of death - attending date
    %   2- those who lost follow-up: date lost follow up - attending date
    %   3- for the rest: censoring date - attending date
    % censoring:
    %   1- those developed outcome: code 1
    %   2- those who died before outcome (death of other causes, i.e. anything except main and ex criteria): code 2
    %   3- those with competing risk before (death or diagnosis) outcome: code 3
    %   4- those after censoring date/lost follow-up before the outcome: code 0

    % those who died before outcome (death of other causes, i.e. anything
    % except main and ex criteria): code 2
    % dc = struct2table(dc);
    if opts.mortality
        dc = death; clear death
        dc.died(:) = true;
        dc.df(:) = "DEATH";
    elseif ~isstruct(dc) % if isstruct, then death data are not querried
         rm_idx = ismember(dc.eid, death.eid);
        if ~isstruct(hesin)
            rm_idx = ismember(dc.eid, hesin.eid) | rm_idx;
        end
       
        dc(rm_idx, :) = [];
    
        % unique ids for death: date of death cannot be different per each
        % icd_10!
        dc = renamevars(dc, "date_of_death", "date");
        dc = groupsummary(dc, "eid", @(x)unique(x), ["date", "regdate", "dob", "source"]);
        dc.GroupCount = [];
        dc.Properties.VariableNames = regexprep(dc.Properties.VariableNames, "^fun1_", "");
        
        % add ex and excode cols
        dc.ex(:) = 2;
        dc.excode(:) = "";
        dc.rawUKB(:) = "";
        dc = movevars(dc, "rawUKB", After="eid");
        idx = true(height(dc), 1);
        dc.excode(idx) = missing; dc.rawUKB(idx) = missing;
    
        % apply censoring for dc table
        idx0 = dc.source == 0 & (cen_date.d.HES < dc.date); %E/W
        idx1 = dc.source == 1 & (cen_date.d.SMR < dc.date); % Scottland 
        dc.censor(:) = 2; % died of any other cause
        dc.date(idx0) = cen_date.d.HES; % replace censoring date for time-to-even calculation
        dc.date(idx1) = cen_date.d.SMR; % same as above
        dc.censor(idx0 | idx1) = 0;
        dc.source = [];
    
        % for calculating time-to-even first reconcile overlapping individuals
        % between death and HESIN registeries.
        if ~isstruct(hesin)
            [f1, f2] = ismember(hesin.eid, death.eid); f2(f2<1) = [];
            hesin.died(:) = false;
        
            if any(f1)       
                assert(~any(death.date(f2) < hesin.date(f1)), "No one can die before HESIN diagnosis!") % no one can die before diagnosis!
                death(f2, :) = []; % we use HESIN since this happened before death for these subjects
                hesin.died(f1) = true; % also they died!
            end
        end
        
        dc = [death; dc];
        dc.died(:) = true;
        
        dc.df(:) = "DEATH";

        if ~isstruct(hesin)
            hesin.df(:) = "HESIN";
            dc = [hesin; dc];
            clear hesin
        end

    else % not death data is requested: only HESIN
        hesin.died(:) = false;
        hesin.df(:) = "HESIN";
        dc = hesin;
        clear hesin
    end
    dc.excode(dc.excode == "") = missing;
    dc.rawUKB(dc.rawUKB == "") = missing;
    % dc.df(:) = "HESIN";
    % dc.df(dc.died) = "DEATH";

    % check if they left the study before the diagnosis (either ex or main outcome): code 0
    if isfield(opts, "dlf")
        [f1, f2] = ismember(dc.eid, opts.dlf.eid); f2(f2<1) = [];
        if any(f1)
            idx_lf_before_outcome = find(f1); % lost follow-up before outcome
            idx = opts.dlf.rawUKB(f2) < dc.date(f1);
            idx2 = idx_lf_before_outcome(idx);
            dc.censor(idx2) = 0;
            dc.date(idx2) = opts.dlf.rawUKB(f2(idx)); % use the date lost-follow up for these individuals
        end
    end
    dc.tt = years(dc.date - dc.regdate);
    dc.tt0 = years(dc.date - dc.dob);
    UKB_STRUCT_ALL = dc;
    clear dc

    return
end

% step 3-merge hesin and death structs
if ~isempty(death.eid) && ~isempty(hesin.eid)
    [idx1, idx2] = ismember(hesin.eid, death.eid);
    idx2(idx2 < 1) = [];

    % add non-overlapping samples
    idx3 = setdiff(1:numel(death.eid), idx2);
    UKB_STRUCT_ALL.eid = [hesin.eid(~idx1); death.eid(idx3)];
    if ~isempty(idx3)
        UKB_STRUCT_ALL.rawUKB = [hesin.rawUKB(~idx1); death.rawUKB(idx3)];
    else
        UKB_STRUCT_ALL.rawUKB = hesin.rawUKB(~idx1);
    end
    UKB_STRUCT_ALL.source = ...
        [repmat("HESIN", sum(~idx1), 1);...
        repmat("DEATH", numel(idx3), 1)];
    
    % join overlapping samples
    if any(idx1)
        UKB_STRUCT_ALL.eid = [UKB_STRUCT_ALL.eid; hesin.eid(idx1)];
        tmpcodes = hesin.rawUKB(idx1) + "," + death.rawUKB(idx2);
        tmpcodes = arrayfun(@(x)join(unique(split(x, ",")), ","), tmpcodes);
        UKB_STRUCT_ALL.rawUKB = [UKB_STRUCT_ALL.rawUKB; tmpcodes];
        UKB_STRUCT_ALL.source = [UKB_STRUCT_ALL.source; 
            repmat("HESIN&DEATH", sum(idx1), 1)];
    end
    
elseif ~isempty(death.eid)  % only death 
    UKB_STRUCT_ALL = death;
    UKB_STRUCT_ALL.source = repmat("DEATH", numel(death.eid), 1);
else % only hesin
    UKB_STRUCT_ALL = hesin;
    UKB_STRUCT_ALL.source = repmat("HESIN", numel(hesin.eid), 1);
end

if ~isempty(UKB_STRUCT_ALL.eid) && ~isempty(duplicates(UKB_STRUCT_ALL.eid))
    error('getHESIN_DEATH: non-unique eids!')
end

UKB_STRUCT_ALL.numericFlag = false;
if all(keyterms ~= "")
    [~, matchedIDX] = ismember(UKB_STRUCT_ALL.rawUKB, keyterms);
    if opts.meaning
        UKB_STRUCT_ALL.termMeaning = keymeanings(matchedIDX);
        UKB_STRUCT_ALL.tag = interms.tags;
    end
end

end %END