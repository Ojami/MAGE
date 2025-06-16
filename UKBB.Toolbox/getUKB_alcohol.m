function UKB_STRUCT_ALL = getUKB_alcohol(inst)

if nargin < 1
    inst = 0; % baseline
end

% Jan. 2020
%% Notes-------------------------------------------------------------------
%           Alcohol             Unit   g
% Red/white wine or champagne   1.5    12
% Beer/cider                    2      16
% Spirits                       1      8
% Fortified wine                1.5    12
% Other (e.g. alcopops)         1.5    12

% Based on paper: DOI: https://doi.org/10.1017/S0033291719000357
%           Alcohol             Unit   g
% Red/white wine or champagne   1.67   13.36
% Beer/cider                    2.3    18.4
% Spirits                       1      8
% Fortified wine                2.5    20
% Other (e.g. alcopops)         -      -

% Note: in several papers, 'other alcoholic drinks' were not considered or
% set to 0 (doi: 10.1093/ije/dyz064,
% https://www.nature.com/articles/s41467-019-12424-x), except this one
% (https://academic.oup.com/cercor/article/29/10/4169/5232542)

% @25NOV2024: major update
%--------------------------------------------------------------------------
%%
% Step 0: read Alcohol_intake_frequency.mat
if ~isfile("aif.mat")
    aif = phenoParser("query", "", "df", "1558;1568;1578;1588;1598;1608;5364;4407;4418;4429;4440;4451;4462", "ukbrap", true, ...
        "save", false, "surv", false, "all", true);
    save("aif.mat", "aif")
else
    aif = load("aif.mat").aif;
end

idx = aif.tag.contains("Alcohol intake frequency. | Instance " + inst);
UKB_STRUCT_ALL = table2struct(aif(idx, :));

% add sex
sex_df = load("Sex.mat").UKB_STRUCT_ALL;
[~, idx] = ismember(UKB_STRUCT_ALL.eid, sex_df.eid);
UKB_STRUCT_ALL.sex = sex_df.termMeaning(idx);

idx_weekly = aif.tag.lower.contains("average weekly") & aif.tag.lower.contains("instance " + inst);
idx_monthly = aif.tag.lower.contains("average monthly") & aif.tag.lower.contains("instance " + inst);
assert(sum(idx_weekly) == sum(idx_monthly), "number of monthly and weekly report files do not match!")
weekly_reports = aif(idx_weekly, :);
monthly_reports = aif(idx_monthly, :);

UKB_STRUCT_ALL = convunint2g(UKB_STRUCT_ALL, weekly_reports, 'week');
UKB_STRUCT_ALL.g_week = UKB_STRUCT_ALL.g_week./7; % Convert to g/day

UKB_STRUCT_ALL = convunint2g(UKB_STRUCT_ALL, monthly_reports, 'month');
UKB_STRUCT_ALL.g_month = UKB_STRUCT_ALL.g_month./30.4375; % Convert to g/day

%==================== Merge weekly/monthly data ===========================
% Merge monthly and weekly data (give precedence to weekly intake data)
% a. If data in weekly reports is -1 or -3, replace it with corresponding
%    monthly data, if not NaN, -1 or -3.
% b. If data in weekly reports is missing (NaN), replace it with
%    corresponding monthly data.
UKB_STRUCT_ALL.g_day = inf(size(UKB_STRUCT_ALL.g_week));
for ii = 1:size(UKB_STRUCT_ALL.g_week, 1)
    for jj = 1:size(UKB_STRUCT_ALL.g_week, 2)
        if isnan(UKB_STRUCT_ALL.g_week(ii, jj))
            UKB_STRUCT_ALL.g_day(ii, jj) = UKB_STRUCT_ALL.g_month(ii, jj);
        elseif UKB_STRUCT_ALL.g_week(ii, jj) < 0 % -1 or -3
            if ~isnan(UKB_STRUCT_ALL.g_month(ii, jj)) && UKB_STRUCT_ALL.g_month(ii, jj) >= 0
                UKB_STRUCT_ALL.g_day(ii, jj) = UKB_STRUCT_ALL.g_month(ii, jj);
            else
                UKB_STRUCT_ALL.g_day(ii, jj) = UKB_STRUCT_ALL.g_week(ii, jj);
            end
        else
            UKB_STRUCT_ALL.g_day(ii, jj) = UKB_STRUCT_ALL.g_week(ii, jj);
        end
    end
end

f_anynan_day = sum(isnan(UKB_STRUCT_ALL.g_day), 2);

f_negative = sum(UKB_STRUCT_ALL.g_day < 0, 2);
f_anynan_day = f_anynan_day + f_negative; % Some rows have only -1/-3/NaN, but not values, should be set to NaN
f_negative_all = f_anynan_day == size(UKB_STRUCT_ALL.g_day, 2); % Some rows have only -1/-3/NaN, but not values, should be set to NaN
f_negative_zero = f_negative < size(UKB_STRUCT_ALL.g_day, 2); % Some values are available (at least one),so -1 or -3 should be set to 0
f_negative = UKB_STRUCT_ALL.g_day < 0;

UKB_STRUCT_ALL.g_day(f_negative & f_negative_zero) = 0;
UKB_STRUCT_ALL.g_day(f_negative & f_negative_all) = NaN;

% All the columns for each individual should be missing. Meaning, even if
% one type of alcohol has a available value, others will be set to 0, and
% will not be removed.
f_nan_day = all(isnan(UKB_STRUCT_ALL.g_day), 2); 

UKB_STRUCT_ALL.g_day(isnan(UKB_STRUCT_ALL.g_day)) = 0; % Set NaN values to 0: only for summation, NaN cells are stored in f_nan_day
UKB_STRUCT_ALL.g_day = sum(UKB_STRUCT_ALL.g_day, 2); % Sum over all alcoholic drinks
%==========================================================================
%===================== Set "Never" or "Special occasions" to non-missing ==
never_idx = UKB_STRUCT_ALL.termMeaning == "Never";
occasional_idx = UKB_STRUCT_ALL.termMeaning == "Special occasions only";
f_nan_day(never_idx | occasional_idx) = false;
% =========================================================================
% ======================== Check unjustified 0 drinkers ===================
% Individuals who already responded as one of these categories of
% overall alcohol intake frequency: "Once or twice a week", "Daily or
% almost daily", "Three or four times a week", but reported 0 g/day of
% alcohol intake, will be excluded due to inconsistency.

strange_zero_idx1 = (UKB_STRUCT_ALL.termMeaning == "Once or twice a week")...
    & ~f_nan_day & (UKB_STRUCT_ALL.g_day == 0);
strange_zero_idx2 = (UKB_STRUCT_ALL.termMeaning == "Daily or almost daily")...
    & ~f_nan_day & (UKB_STRUCT_ALL.g_day == 0);
strange_zero_idx3 = (UKB_STRUCT_ALL.termMeaning == "Three or four times a week")...
    & ~f_nan_day & (UKB_STRUCT_ALL.g_day == 0);
strange_zero_idx = strange_zero_idx1 | strange_zero_idx2 | strange_zero_idx3;
fprintf('%d samples were removed due to unjustified 0 g/day for once or twice a week\n', ...
    sum(strange_zero_idx1))
fprintf('%d samples were removed due to unjustified 0 g/day for daily or almost daily\n', ...
    sum(strange_zero_idx2))
fprintf('%d samples were removed due to unjustified 0 g/day for three or four times a week\n', ...
    sum(strange_zero_idx3))

UKB_STRUCT_ALL.eid(strange_zero_idx) = [];
UKB_STRUCT_ALL.rawUKB(strange_zero_idx) = [];
UKB_STRUCT_ALL.g_day(strange_zero_idx) = [];
UKB_STRUCT_ALL.g_week(strange_zero_idx, :) = [];
UKB_STRUCT_ALL.g_month(strange_zero_idx, :) = [];
UKB_STRUCT_ALL.sex(strange_zero_idx) = [];
UKB_STRUCT_ALL.termMeaning(strange_zero_idx) = [];
f_nan_day(strange_zero_idx) = [];
% f_nan_day(strange_zero_idx) = true;
%==========================================================================
%============================= Outlier removal ============================
% Remove outliers based on this paper:
% https://www.nature.com/articles/mp2017153 Individuals with an alcohol
% consumption quantity deviating >5 s.d. from their sex-specific mean will
% be excluded.

% First exclude all f_nan_day indices because they've been "wrongly" set to 0.
TEMP_STRUCT = UKB_STRUCT_ALL;
TEMP_STRUCT.eid(f_nan_day) = [];
TEMP_STRUCT.rawUKB(f_nan_day) = [];
TEMP_STRUCT.termMeaning(f_nan_day) = [];
TEMP_STRUCT.sex(f_nan_day) = [];
TEMP_STRUCT.g_day(f_nan_day) = [];
TEMP_STRUCT.g_month(f_nan_day) = [];
TEMP_STRUCT.g_week(f_nan_day) = [];

sex_cat = unique(TEMP_STRUCT.sex);
if numel(sex_cat) > 2
    error('Wrong sex categories')
end
f_sex1 = find(TEMP_STRUCT.sex == sex_cat(1));
f_sex2 = find(TEMP_STRUCT.sex == sex_cat(2));
alc_sex1 = TEMP_STRUCT.g_day(f_sex1);

% outliers_sex1 = isoutlier(alc_sex1, 'quartiles');
outliers_sex1 = (alc_sex1 > (mean(alc_sex1) + 5*std(alc_sex1)));
% outliers_sex1 = alc_sex1 > 20;

alc_sex2 = TEMP_STRUCT.g_day(f_sex2);

% outliers_sex2 = isoutlier(alc_sex2, 'quartiles');
outliers_sex2 = (alc_sex2 > (mean(alc_sex2) + 5*std(alc_sex2)));
% outliers_sex2 = alc_sex2 > 30;

fprintf('%d outliers were found in %s\n', sum(outliers_sex1), sex_cat(1))
fprintf('%d outliers were found in %s\n', sum(outliers_sex2), sex_cat(2))
alc_outliers = union(f_sex1(outliers_sex1), f_sex2(outliers_sex2));
alc_outliers = TEMP_STRUCT.eid(alc_outliers);
clear TEMP_STRUCT
outlier_idx = ismember(UKB_STRUCT_ALL.eid, alc_outliers);
UKB_STRUCT_ALL.eid(outlier_idx) = [];
UKB_STRUCT_ALL.rawUKB(outlier_idx) = [];
UKB_STRUCT_ALL.termMeaning(outlier_idx) = [];
UKB_STRUCT_ALL.sex(outlier_idx) = [];
UKB_STRUCT_ALL.g_month(outlier_idx, :) = [];
UKB_STRUCT_ALL.g_week(outlier_idx, :) = [];
UKB_STRUCT_ALL.g_day(outlier_idx, :) = [];
f_nan_day(outlier_idx, :) = [];
%==========================================================================
%===================== Set median values to unknown data ==================
% For participants who had unknown grams/day of alcohol but who
% reported one of these categories of overall alcohol intake frequency:
% "Once or twice a week", "One to three times a month", "Daily or almost
% daily", "Three or four times a week", the median value (g/day) from their
% category was assigned (doi: 10.1093/ije/dyz064)
once_twice_week_idx = (UKB_STRUCT_ALL.termMeaning == "Once or twice a week") & f_nan_day;
daily_idx = (UKB_STRUCT_ALL.termMeaning == "Daily or almost daily") & f_nan_day;
three_four_week_idx = (UKB_STRUCT_ALL.termMeaning == "Three or four times a week") & f_nan_day;
one_three_month_idx = (UKB_STRUCT_ALL.termMeaning == "One to three times a month") & f_nan_day;

if any(once_twice_week_idx(:))
    fprintf('%d missing values for once or twice a week were replaced by median of known values\n', sum(once_twice_week_idx(:)))
    this_category_idx = (UKB_STRUCT_ALL.termMeaning == "Once or twice a week") & ~f_nan_day;
    UKB_STRUCT_ALL.g_day(once_twice_week_idx) = median(UKB_STRUCT_ALL.g_day(this_category_idx));
end
if any(daily_idx(:))
    fprintf('%d missing values for daily or almost daily were replaced by median of known values\n', sum(daily_idx(:)))
    this_category_idx = (UKB_STRUCT_ALL.termMeaning == "Daily or almost daily") & ~f_nan_day;
    UKB_STRUCT_ALL.g_day(daily_idx) = median(UKB_STRUCT_ALL.g_day(this_category_idx));
end
if any(three_four_week_idx(:))
    fprintf('%d missing values for three or four times a week were replaced by median of known values\n', sum(three_four_week_idx(:)))
    this_category_idx = (UKB_STRUCT_ALL.termMeaning == "Three or four times a week") & ~f_nan_day;
    UKB_STRUCT_ALL.g_day(three_four_week_idx) = median(UKB_STRUCT_ALL.g_day(this_category_idx));
end
if any(one_three_month_idx(:))
    fprintf('%d missing values for one to three times a month were replaced by median of known values\n', sum(one_three_month_idx(:)))
    this_category_idx = (UKB_STRUCT_ALL.termMeaning == "One to three times a month") & ~f_nan_day;
    UKB_STRUCT_ALL.g_day(one_three_month_idx) = median(UKB_STRUCT_ALL.g_day(this_category_idx));
end

%==========================================================================
f_nan_day = f_nan_day & (UKB_STRUCT_ALL.g_day == 0);
UKB_STRUCT_ALL.g_day(f_nan_day) = [];
UKB_STRUCT_ALL.eid(f_nan_day) = [];
UKB_STRUCT_ALL.termMeaning(f_nan_day) = [];

UKB_STRUCT_ALL.rawUKB = UKB_STRUCT_ALL.g_day;
UKB_STRUCT_ALL = rmfield(UKB_STRUCT_ALL, {'g_week', 'sex', 'g_month', 'g_day'});
UKB_STRUCT_ALL.tag = "Alcohol consumption (g/day)";
if inst == 2, UKB_STRUCT_ALL.tag = UKB_STRUCT_ALL.tag + " MRI"; end
UKB_STRUCT_ALL.numericFlag = true;
% UKB_STRUCT_ALL.termMeaning = '';

save(matlab.lang.makeValidName(UKB_STRUCT_ALL.tag) + ".mat", 'UKB_STRUCT_ALL')
end

%% subfunctions -----------------------------------------------------------
function alcoh_freq = convunint2g(alcoh_freq, report_data, timeFlag)

if strcmp(timeFlag, 'week')
    alcoh_freq.g_week = nan(numel(alcoh_freq.eid), height(report_data));
else
    alcoh_freq.g_month = nan(numel(alcoh_freq.eid), height(report_data));
end

for k = 1:height(report_data)
    xukb = table2struct(report_data(k, :));
    % xukb.rawUKB = double(xukb.rawUKB);
    % f_calc = ismissing(xukb.termMeaning);
    % if ~all(f_calc == (xukb.rawUKB > -1))
    %     error('Someting went wrong!')
    % end
    f_calc = xukb.rawUKB >= 0; % -1	Do not know, -3	Prefer not to answer
    
    tag = xukb.tag.lower;
    if tag.contains("beer plus cider")
        xukb.rawUKB(f_calc) = xukb.rawUKB(f_calc).*16;
    elseif tag.contains("champagne plus white wine")
        xukb.rawUKB(f_calc) = xukb.rawUKB(f_calc).*12;
    elseif tag.contains("fortified wine")
        xukb.rawUKB(f_calc) = xukb.rawUKB(f_calc).*12;
    elseif tag.contains("red wine")
        xukb.rawUKB(f_calc) = xukb.rawUKB(f_calc).*12;
    elseif tag.contains("spirits")
        xukb.rawUKB(f_calc) = xukb.rawUKB(f_calc).*8;
    elseif tag.contains("other alcoholic drinks")
        xukb.rawUKB(f_calc) = xukb.rawUKB(f_calc).*12;
    else
        error("unknown tag")
    end
    [nidx1, nidx2] = ismember(alcoh_freq.eid, xukb.eid);
    nidx2(nidx2 < 1) = [];
    if strcmp(timeFlag, 'week')
        alcoh_freq.g_week(nidx1, k) = xukb.rawUKB(nidx2);
    else
        alcoh_freq.g_month(nidx1, k) = xukb.rawUKB(nidx2);
    end
end
end