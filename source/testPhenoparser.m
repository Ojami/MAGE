function testPhenoparser
clc

% a tutorial on how to use phenoParser function to define time-to-event
% traits.

% trait: CLD + HCC

% exclusion: 
% 1-non-cancer self reported CLD: df 20002: instance 0 (baseline)
ex1 = [1136, 1155, 1156, 1157, 1158, 1159, 1578, 1579, 1580, 1581, 1582, 1506, 1507, 1604, 1507, 1141];
ex1 = join(string(ex1), ";");

%2-cancer self reported blood cancer: df 20001: instance 0 (baseline)
ex2 = [1074, 1070, 1058, 1056, 1047, 1048, 1050, 1051, 1052, 1053, 1055];
ex2 = join(string(ex2), ";");


%3-blood cancer from cancer registry df 40006 (coding 19): all instances 
ex3 = "C" + (81:96) + (0:9)';
ex3 = join(ex3(:), ";");

%4-viral hep. from coding 19: all instance
ex4 = ["B18", "B180", "B181", "B182","B188", "B189", "B19", "B190","B199"];
ex4 = join(ex4, ";");

% inclusion (coding 19): i.e. in HESIN, death, and cancer reported
inc1 = ["C220", "K70", "K700", "K701","K702","K703","K704", "K709", "K721","K729","K730","K731", "K732","K738","K739", "K740", "K741", "K742", "K746", "K760", "K766", "K767", "K768", "K769", "I850", "I859"];
inc1 = join(inc1, ";");

% example 0: HCC: C220
ex = struct;

% example 1
ex = struct;
ex.df20002_0 = ex1;
ex.df20001_0 = ex2;
ex.cd19 = ex3 + ";" + ex4; 
t = phenoParser("query", inc1, "coding", "19", "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'remove', ex, 'desc', "surv", "save", false, "ukbrap", true);

% example 2
ex5 = ex3 + ";" + ex4;
ex = struct("df20001_0", ex2, "cd19", ex5);

t = phenoParser("query", [ex1, inc1], "coding", ["", "19"], "df", ["20002", ""], "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'remove', ex, 'desc', "surv");

% example 3
ex5 = ex3 + ";" + ex4 + ";" + inc1;
ex = struct("df20001_0", ex2, "cd19", ex5);
t = phenoParser("query", ex1, "coding", "", "df", "20002", "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'remove', ex, 'desc', "surv");

% example 4
ex = struct("df20001_0", ex2);
t = phenoParser("query", ex1, "coding", "", "df", "20002", "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'remove', ex, 'desc', "surv");

% example 5
ex = struct("df20001_0", ex2, "df20002_0", ex1);
t = phenoParser("query", inc1, "coding", "19", "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'remove', ex, 'desc', "surv");

% example 6: no exclusion criteria
t = phenoParser("query", [ex1, inc1], "coding", ["", "19"], "df", ["20002", ""], "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'desc', "surv");

% % example 7: no exclusion criteria
t = phenoParser("query", [ex1, ex2], "df", ["20002", "20001"], "all", true, ...
    'phenodir', pwd, "tag", "CLD_surv", ...
    'desc', "surv");

% % example 8
inc1 = ["K700", "K701","K702","K703","K704", "K709", "K721","K729","K730","K731", "K732","K738","K739", "K740", "K741", "K742", "K746", "K760", "K766", "K767", "K768", "K769", "I850", "I859"];
inc1 = join(inc1, ";");
ex = struct;
ex.cd19 = ex4; 
t = phenoParser("query", inc1, "coding", "19", "all", true, ...
    "tag", "Chronic liver disease", ...
    'remove', ex, "save", true, "surv", true);

% example 9: overall mortality
t = phenoParser("query", "mortality", "all", true, ...
    "tag", "Mortality2", "save", true, "surv", true);

s.snp = ["rs12238997", "rs2642438", "rs371890604", "rs148120343"]; s.chr = ["1", "1", "1", "1"];


gwasrunner(s, "trait", ["Mortality.surv", "HCC", "ALT", "ChronicLiverDisease.surv"], "output", ...
    fullfile(pwd, "myfirst", "test"), "covar", ["Age", "Sex", "PNPLA3_rs738409", "BMI"], ...
    "catCovar", "PNPLA3_rs738409")


end % END