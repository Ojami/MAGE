function UKB_HESIN_DEATH_2mat(opts)
% parses hesin and death data for phenoParser function. These files can be
% dowlonaded from UKBB application/project page, under data tables.
% 
% Oveis Jamialahmadi, University of Gothenburg, March 2020.
% 
% @20FEB2024: major update to accept time-to-event data (see phenoParser
%             also)
% @21FEB2024: merges date of attending the assessment center df53 with
%             heasin and death data.
% @03JUNE2024: for genetic studies, age-of-onset (AOO) should be used.
%              Therefore we add also date of birth.
% @25JULY2024: Use UKB RAP to fetch parquet files for hesin/death data. For
%              more instructions see fetch_hesin.ipynb in /MAGE/dxtoolkit
% @20SEP2024: dx_extract_dataset can fetch HESIN and DEATH tables with
%             functionality implemented in phenoParser. So, it's
%             recommended NOT to run this function separately. 'dir'
%             argument was added for output directory: default is
%             UKB_HESIN_DEATH_2mat.m directory.

arguments
    opts.hesin {mustBeFile}
    opts.hesin_diag {mustBeFile}
    opts.death {mustBeFile}
    opts.death_cause {mustBeFile}
    opts.workers (1,1) double = 0 % to invoke parallel pool
    opts.instance (1,1) double = 0 % 0 is baseline, 2 is imaging visit: 'regdate' depends on this.

    %@20SEP2024
    opts.dir {mustBeTextScalar} = fileparts(which("UKB_HESIN_DEATH_2mat.m"))
end

fclose('all');


% both 'hesin' and 'hesin_diag' must be provided
if ~all(isfield(opts, ["hesin", "hesin_diag", "death", "death_cause"]))
    error("all files must be provided: hesin, hesin_diag, death and death_cause!")
end

if ~isfolder(opts.dir), mkdir(opts.dir); end 

if isempty(gcp('nocreate')) && opts.workers > 0
    parpool("Processes", opts.workers);
    % parpool("Threads", opts.workers);
end

if opts.workers > 1
    opts.workers = true;
else
    opts.workers = false;
end

% read date attending the assessment center (for prospective analysis)
da = table2struct(phenoParser("query", "", "df", "53", "save", false, ...
    "instance", opts.instance, "surv", false, "all", true));

% read data of birht based on year and month of birth
yob = table2struct(phenoParser("query", "", "df", "34", "save", false, ...
    "surv", false, "all", true));
mob = table2struct(phenoParser("query", "", "df", "52", "save", false, ...
    "surv", false, "all", true));
[f1, f2] = ismember(yob.eid, mob.eid); f2(f2<1) = [];
dob = table(yob.eid(f1), datetime(yob.rawUKB(f1), double(mob.rawUKB(f2)), ones(sum(f1), 1), ...
    Format="dd/MM/uuuu"), VariableNames=["eid", "dob"]);


% find the censoring date according to UKBB recommendations
% https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates
% ref: https://link.springer.com/article/10.1186/s12916-023-02772-3
% fetchMore = true;
% [death_found, hesin_found] = deal(false);
regions = ["England", "Scot", "Wales"];
tags = ["HES", "SMR", "PEDW"];
% ct = 1;
cd = struct; % censored date struct

%@24SEP2024: use pandas: readtable throws some unexpected errors;
pd = py.importlib.import_module("pandas");
url = "https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates";
tabs = cell(pd.read_html(url));

for k = 1:numel(tabs)
    tab = table(tabs{k});
    tab = convertvars(tab, 1:width(tab), @string);
    tab.(width(tab)) = erase(tab.(width(tab)), "*");

    if tab.(1)(1).lower == "death"

        for j = 1:numel(regions)
            idx = tab.(1).contains(regions(j));
            if any(idx)
                tab.(width(tab))(idx) = tab.(width(tab))(idx).strip;
                cd.d.(tags(j)) = datetime(tab.(width(tab))(idx).strip("both", char(160)).replace(whitespacePattern, "-"), Format="dd/MM/uuuu");
            end
        end

    elseif tab.(1)(1).lower.startsWith("hospital admission")

        for j = 1:numel(regions)
            idx = tab.(1).contains(regions(j));
            if any(idx)
                tab.(width(tab))(idx) = tab.(width(tab))(idx).strip;
                cd.h.(tags(j)) = datetime(tab.(width(tab))(idx).strip("both", char(160)).replace(whitespacePattern, "-"), Format="dd/MM/uuuu");
            end
        end

    end
end
        
% while fetchMore
%     try
%         % tab = readtable("https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates", ...
%         %     FileType="html",ReadVariableNames=true, ...
%         %     TableSelector="//TABLE[" + ct + "]", ...
%         %     VariableNamingRule="preserve");
% 
%         cols = colnames(tab);
%         idx = find(cols.lower.contains(["censoring", "date"]));
%         tab(ismissing(tab.(idx)), :) = [];
%         tab.(idx) = erase(tab.(idx), "*");
%         tab.(idx) = datetime(tab.(idx).strip.replace(whitespacePattern, "-"), Format="dd/MM/uuuu");
% 
%         if any(cols.lower == "death")
%             dcol = cols(cols.lower == "death");
%             tab = tab(:, [cols(idx), dcol]);
% 
%             for k = 1:numel(regions)
%                 f = tab.(2).contains(regions(k));
%                 if any(f)
%                     cd.d.(tags(k)) = tab.(1)(f);
%                 end
%             end
%             death_found = true;
% 
%         elseif any(cols.lower.contains("hospital admissions"))
%             dcol = cols(cols.lower.contains("hospital admissions"));
%             tab = tab(:, [cols(idx), dcol]);
% 
%             for k = 1:numel(regions)
%                 f = tab.(2).contains(regions(k));
%                 if any(f)
%                     cd.h.(tags(k)) = tab.(1)(f);
%                 end
%             end
%             hesin_found = true;
%         end
% 
%         ct = ct + 1;
% 
%         if death_found && hesin_found
%             fetchMore = false;
%         end
% 
%     catch % no more tables
%         fetchMore = false;
%     end
% end


% read death/death_cause files --------------------------------------------
if opts.death.endsWith(".parquet")
    ds = parquetDatastore(opts.death, "OutputType", "table", "VariableNamingRule", "preserve");
else
    ds = tabularTextDatastore(opts.death, "TextType", "string", ...
        "DatetimeType", "text", "NumHeaderLines", 0, ...
        "SelectedVariableNames", ["eid", "ins_index", "date_of_death", "dsource"], ...
        "VariableNamingRule", "preserve");
end

if opts.workers
    ds = gather(tall(ds));
else
    ds = readall(ds);
end

if opts.death.endsWith(".parquet")
    ds.eid = double(ds.eid);
    ds.date_of_death = string(ds.date_of_death);
end

if opts.death_cause.endsWith(".parquet")
    dc = parquetDatastore(opts.death_cause, "OutputType", "table", "VariableNamingRule", "preserve");
else
    dc = tabularTextDatastore(opts.death_cause, "TextType", "string", ...
        "DatetimeType", "text", "NumHeaderLines", 0, ...
        "SelectedVariableNames", ["eid", "ins_index", "cause_icd10"], ...
        "VariableNamingRule", "preserve");
end

if opts.workers
    dc = gather(tall(dc));
else
    dc = readall(dc);
end

if opts.death_cause.endsWith(".parquet")
    dc.eid = double(dc.eid);
end


ds = outerjoin(ds, dc, "Keys", ["eid", "ins_index"], "MergeKeys", true, "Type", "left");
clear("dc")
ds.ins_index = [];

% add date of attending the assessment center
[f1, f2] = ismember(ds.eid, da.eid);
if ~any(f1)
    fprintf("WARNING: there are some eids in DEATH table which did not register in UKBB!");
    fprintf("\tThose participants will be removed from the final HESIN table")
end
f2(f2<1) = [];
ds(~f1, :) = [];
ds.regdate = da.rawUKB(f2);

[~, f2] = ismember(ds.eid, dob.eid); f2(f2<1) = [];
ds.dob = dob.dob(f2);
ds.date_of_death = datetime(ds.date_of_death, Format="dd/MM/uuuu");

% add censoring info (E/W = 0, Scottland = 1, rest = -1)
e_w_idx = ds.dsource == "E/W"; % England and Wales are together in DEATH tables
s_idx = ds.dsource == "SCOT";
ds.source(:) = int8(-1);
ds.source(e_w_idx) = int8(0);
ds.source(s_idx) = int8(1);
ds.dsource = [];

% parquetwrite(regexprep(opts.death_cause, '.txt$', '.parquet'), ds)
ds = table2struct(ds, 'ToScalar', true);  
[~, death_name] = fileparts(opts.death_cause);
save(fullfile(opts.dir, death_name + ".mat"), '-struct', 'ds');
clear("ds")


% read hesin/hesin_diag files ---------------------------------------------
if opts.hesin.endsWith(".parquet")
    hs = parquetDatastore(opts.hesin, "OutputType", "table", "VariableNamingRule", "preserve");
else
    hs = tabularTextDatastore(opts.hesin, "TextType", "string", ...
        "DatetimeType", "text", "NumHeaderLines", 0, ...
        "SelectedVariableNames", ["eid", "ins_index", "dsource", "epistart", "admidate"], ...
        "VariableNamingRule", "preserve");
end

if opts.workers
    hs = tall(hs);
else
    hs = readall(hs);
end

if opts.hesin.endsWith(".parquet")
    hs.eid = double(hs.eid);
    hs.epistart = string(hs.epistart);
    hs.admidate = string(hs.admidate);
end

% @21FEB2024: text is lighter 
% replace missing epistart with admidate (if available)
if isdatetime(hs.epistart) % parquet
    idx = ismissing(hs.epistart);
else
    idx = hs.epistart == "" | ismissing(hs.epistart);
end
hs.epistart(idx) = hs.admidate(idx);
hs.admidate = [];

% convert epistart and date
hs.epistart = datetime(hs.epistart, Format="dd/MM/uuuu");

if opts.workers, hs = gather(hs); end

if opts.hesin_diag.endsWith(".parquet")
    hd = parquetDatastore(opts.hesin_diag, "OutputType", "table", "VariableNamingRule", "preserve");
else
    hd = tabularTextDatastore(opts.hesin_diag, "TextType", "string", ...
        "DatetimeType", "text", "NumHeaderLines", 0, ...
        "SelectedVariableNames", ["eid", "ins_index", "diag_icd10", "diag_icd9"], ...
        "VariableNamingRule", "preserve");
    hd.SelectedFormats = ["%f", "%f", "%q", "%q"];
end

if opts.workers
    hd = gather(tall(hd));
else
    hd = readall(hd);
end

if opts.hesin_diag.endsWith(".parquet")
    hd.eid = double(hd.eid);
end

hs = outerjoin(hd, hs, "Keys", ["eid", "ins_index"], "MergeKeys", true, "Type", "left");
clear("hd")
hs.ins_index = [];

% add date of attending the assessment center
[f1, f2] = ismember(hs.eid, da.eid);
if ~any(f1)
    fprintf("WARNING: there are some eids in HESIN table which did not register in UKBB!");
    fprintf("\tThose participants will be removed from the final HESIN table")
end
f2(f2<1) = [];
hs(~f1, :) = [];
hs.regdate = da.rawUKB(f2);

[~, f2] = ismember(hs.eid, dob.eid); f2(f2<1) = [];
hs.dob = dob.dob(f2);

% add cencoring date (England/HES = 0, Scottland/SMR = 1, Wales/PEDW = 2, rest = -1)
hs.source(:) = int8(-1);
for k = 1:numel(tags)
    idx = hs.dsource == tags(k);
    hs.source(idx) = (k -1);
end
hs.dsource = [];

% parquetwrite(regexprep(opts.hesin_diag, '.txt$', '.parquet'), hs)
hs = table2struct(hs, 'ToScalar', true);  
[~, hesin_name] = fileparts(opts.hesin_diag);
save(fullfile(opts.dir, hesin_name + ".mat"), '-struct', 'hs');
clear("hs")

% save censoring info
save(fullfile(opts.dir, "censoring_info.mat"), "cd")

% construct the variable mapper -------------------------------------------
variableMapper = cell(13, 2);
variableMapper(:, 2) = {'x41202_'; 'x41204_'; 'x41270_'; 'x40001_'; ...
    'x40002_'; 'x41271_'; 'x41203_'; 'x41205_'; ...
    'x41262_'; 'x41280_'; 'x40000_'; 'x41281_'; 'x41263_'; };
datev = datetime('now', 'Format', 'dd-MMM-uuuu');
save(fullfile(opts.dir, 'variableMapper.mat'), 'variableMapper', 'datev')

toc
end %END