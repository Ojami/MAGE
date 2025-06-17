function getAgeAtInstance(instance)
% uses dfs 34 (date of birth) and 52 (month of birth) to infer age at
% follow-ups
% was performed using df 53 (date of attending assessment center) array 2
% first imaging visit
% Oveis Jamialahmadi. University of Gothenburg.

if nargin < 1, instance = 2; end

if instance == 3
    desc = "first repeat imaging visit (2019+)";
    tag = " at i3";
elseif instance == 2
    desc = "imaging visit (2014+)";
    tag = " at MRI";
elseif instance == 1
    desc = "First repeat assessment visit (2012-13)";
    tag = " at i1";
elseif instance == 0
    desc = "Initial assessment visit (2006-2010)";
    tag = "";
else
    error("values can be either: 0, 1, 2 or 3!")
end

% year of birth
yob = table2struct(phenoParser(query="", ...
    df="34", ...
    save=false, ...
    surv=false, ...
    ukbrap=true));

% month of birth
mob = table2struct(phenoParser(query="", ...
    df="52", ...
    save=false, ...
    surv=false, ...
    ukbrap=true));

% date of attending
doa = table2struct(phenoParser(query="", ...
    df="53", ...
    instance=instance, ...
    save=false, ...
    surv=false, ...
    ukbrap=true, ...
    all=true)); 

% date of birth
[f1, f2] = ismember(yob.eid, mob.eid); f2(f2<1) = [];
dob = table(yob.eid(f1), datetime(yob.rawUKB(f1), ...
    double(mob.rawUKB(f2)), ones(sum(f1), 1)), ...
    'VariableNames', {'eid', 'dob'});

% keep only those attended first imaging visit
[f1, f2] = ismember(dob.eid, doa.eid); f2(f2<1) = [];
dob(~f1, :) = [];
dob.doa = doa.rawUKB(f2);

dob.agemri = years(dob.doa - dob.dob);

% write to a new UKB struct
UKB_STRUCT_ALL.eid = dob.eid;
UKB_STRUCT_ALL.rawUKB = dob.agemri;
UKB_STRUCT_ALL.numericFlag = true;
UKB_STRUCT_ALL.termMeaning = '';
UKB_STRUCT_ALL.tag = "Age" + tag;
UKB_STRUCT_ALL.info.date = mob.info.date;
UKB_STRUCT_ALL.info.df = "34,52,53";
UKB_STRUCT_ALL.info.dfraw = "x34_0_0,x52_0_0," + doa.info.dfraw;
UKB_STRUCT_ALL.info.basket = mob.info.basket;
UKB_STRUCT_ALL.info.desciption = "difference between year/month of birth and date of attending assessment center for " + desc + " instance " + instance;

pth = fullfile(fileparts(which("phenoParser.m")), "UKB_PHENO"); % default path to save phenotypes
out_name = fullfile(pth, matlab.lang.makeValidName(UKB_STRUCT_ALL.tag) + ".mat");
save(out_name, "UKB_STRUCT_ALL")

fprintf("saved to %s\n", out_name)


end % END