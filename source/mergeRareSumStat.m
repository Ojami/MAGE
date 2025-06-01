function mergeRareSumStat(opts)
% merges seperate summari stat files for rare variant analyses from
% gwasrunner. Files can either be generated for single or multiple genes.
% They can contain different phenotypes (phenome-wide), or/and cover
% several schemes (e.g. LoF, missense,...). Files must be named as
% gene.scheme.xlsx (e.g. GPAM.LoF.xlsx). 
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, Jan 2022.

arguments
    opts.files {mustBeFile}
    opts.wd {mustBeFolder} = pwd
    opts.output {mustBeTextScalar} = "jointStats.xlsx";
%     opts.multigene (1,1) logical = true % wether summary stat files are for single (false) or multiple genes (true) 
end

if ~endsWith(opts.output, [".xlsx", ".xls"])
    opts.output = opts.output + ".xslx";
end

if ~isfield(opts, "files") || any(ismissing(opts.files) | (opts.files == ""))
    opts.files = string({dir(fullfile(opts.wd, '*.xlsx')).name}.');
    opts.files(opts.files == opts.output + ".xlsx" | startsWith(opts.files, "~$")) = [];
end

if isempty(opts.files)
    error('mergeRareSumStat:noFile', 'no summary stat file was found!')
end

% check if files belong to multigene scheme
tmpgenes = unique(extractBefore(opts.files, '.'));
if numel(tmpgenes) > 1
    opts.multigene = true;
else
    opts.multigene = false;
end

if opts.multigene % summary stats for multiple genes (each file can contain different traits)
    % files must be named as gene.scheme.xlsx
    % rows of output table contain genes with each trait stored in a separate sheet
    % if there are different approaches, each approach will be written to a separate file
    opts.genes = extractBefore(opts.files, '.');
    opts.tags = unique(extractBetween(opts.files, opts.genes + ".", ".xlsx")); % different approaches
    opts.ugenes = unique(opts.genes, 'stable');

    for i = 1:numel(opts.tags) % 1 file per each tag (scheme)
        inopt = opts;
        idx = startsWith(opts.files, opts.ugenes + ("."|"_"|"-") + opts.tags(i));
        inopt.files = opts.files(idx);
        inopt.output = erase(inopt.output, (".xlsx"|".xls") + letterBoundary('end'));
        inopt.output = inopt.output + "." + opts.tags(i) + ".xlsx"; % pick a new name
        
        % check traits
        inopt.sheets = sheetnames(fullfile(inopt.wd, inopt.files(1))); % all files should have the same shape
        tmp = readtable(fullfile(inopt.wd, inopt.files(1)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', inopt.sheets(contains(lower(inopt.sheets), "skat")));
        inopt.pheno = tmp.Pheno(~ismissing(tmp.Pheno));
        
        inopt = extractSheets(inopt); % extract data from files 

        for j = 1:numel(inopt.pheno)
            % each trait is written to a spearate sheet
            phenoOpt = inopt;
            phenoOpt.skat = phenoOpt.skat(:, j);
            phenoOpt.collapse = phenoOpt.collapse(:, j);
            phenoOpt.sheetname = inopt.pheno(j);
            if isfield(phenoOpt, 'overall')
                phenoOpt.overall = phenoOpt.overall(:, j);
            end
            writeInnerStats(phenoOpt)
        end

    end

else
    writeInnerStats(opts)
end

end % END

%% ========================================================================
function writeInnerStats(opts)
sz = numel(opts.files);
sheets = sheetnames(fullfile(opts.wd, opts.files(1))); % all files are assumed to be for the same gene

if ~opts.multigene % each gene has a different info sheet
    res.info = readtable(fullfile(opts.wd, opts.files(1)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "info")));
end
res.overall = table(nan(sz, 1), nan(sz, 1), nan(sz, 1), nan(sz, 1), nan(sz, 1),...
    'VariableNames', {'N', 'N (case)', 'MAC.case', 'MAC.control', 'β|OR'});
methods = [];
for i = 1:sz
    
    if ~opts.multigene
        sheets = sheetnames(fullfile(opts.wd, opts.files(i)));
        opts.skat{i, 1} = readtable(fullfile(opts.wd, opts.files(i)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "skat")));
    end
    methods = unique([methods; opts.skat{i}.Method]);

    if any(contains(lower(sheets), "overall"))
        if ~opts.multigene
            overall = readtable(fullfile(opts.wd, opts.files(i)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "overall")));
        else
            overall = opts.overall{i};
        end
        effCol = contains(overall.Properties.VariableNames, {'OR', 'β'});
        res.overall.("β|OR")(i) = overall{1, effCol};
    end

    if ~opts.multigene
        collapse = readtable(fullfile(opts.wd, opts.files(i)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "collapse")));
    else
        collapse = opts.collapse{i};
    end

    idx = contains(lower(collapse.Properties.VariableNames), 'model');
    idx = contains(lower(collapse{:, idx}), 'dominant');
    
    if ~any(contains(lower(sheets), "overall")) % use beta/OR in collapse sheet: this beta/OR comes from a GLM and not Firth (which is found in overall sheet)
        effCol = contains(collapse.Properties.VariableNames, {'OR', 'β'});
        res.overall.("β|OR")(i) = collapse{idx, effCol};
    end

    res.overall.N(i) = collapse.N(idx); % for dominant model
    res.overall.("N (case)")(i) = collapse.("N (case)")(idx);
    res.overall.("MAC.case")(i) = collapse.("MAC.case")(idx);
    res.overall.("MAC.control")(i) = collapse.("MAC.control")(idx);
end

res.skat = splitvars(table("", nan(1, numel(methods) + 1)));
if opts.multigene
    res.skat.Properties.VariableNames = ["Gene", methods.', "N.markers"];
else
    res.skat.Properties.VariableNames = ["Pheno", methods.', "N.markers"];
end

res.skat = repmat(res.skat, sz, 1);
for i = 1:sz
    [~, idx] = ismember(opts.skat{i}.Method, methods);
    if ~opts.multigene
        res.skat.Pheno(i) = opts.skat{i}.Pheno(1);
    end
    res.skat.("N.markers")(i) = opts.skat{i}.("N.marker")(1);
    for j = 1:numel(idx)
        res.skat.(methods(idx(j)))(i) = opts.skat{i}.P(j);
    end
end

if opts.multigene
    res.skat.Gene = extractBefore(opts.files, ".");
end

if ~opts.multigene 
    writetable(res.info, opts.output, 'Sheet', 'info', 'WriteVariableNames', true, 'AutoFitWidth', true);
    writetable([res.skat, res.overall], opts.output, 'Sheet', 'SKAT', 'WriteVariableNames', true, 'AutoFitWidth', true);
end
    if opts.sheetname.strlength > 31 % sheet's name cannot exceed 31 characters
        opts.sheetname = opts.sheetname.char;
        opts.sheetname = string(opts.sheetname(1:31));
    end
    writetable([res.skat, res.overall], opts.output, 'Sheet', opts.sheetname, 'WriteVariableNames', true, 'AutoFitWidth', true);
end

%% ------------------------------------------------------------------------
function opts = extractSheets(opts)
% reads sheets from files to be used later in writeInnerStats function
% (only for multigene flag)
sz = numel(opts.files);
for i = 1:sz
    sheets = sheetnames(fullfile(opts.wd, opts.files(i)));
    opts.skat{i, 1} = readtable(fullfile(opts.wd, opts.files(i)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "skat")));
    opts.skat{i}(all(ismissing(opts.skat{i}), 2), :) = []; % remove blank rows between traits
    if any(contains(lower(sheets), "overall"))
        opts.overall{i, 1} = readtable(fullfile(opts.wd, opts.files(i)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "overall")));
        opts.overall{i}(all(ismissing(opts.overall{i}), 2), :) = [];
    end
    opts.collapse{i, 1} = readtable(fullfile(opts.wd, opts.files(i)), 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Sheet', sheets(contains(lower(sheets), "collapse")));
    opts.collapse{i}(all(ismissing(opts.collapse{i}), 2), :) = [];
end

% check the presence of multiple phenotypes
if numel(opts.pheno) > 1
    idx_skat = [find(ismember(opts.skat{1}.Pheno, opts.pheno)); height(opts.skat{1})+1];
    if any(contains(lower(sheets), "overall"))
        idx_overall = [find(ismember(opts.overall{1}.Pheno, opts.pheno)); height(opts.overall{1})+1];
    end
    idx_collapse = [find(ismember(opts.collapse{1}.Pheno, opts.pheno)); height(opts.collapse{1})+1];
    
    [skat, overall, collapse] = deal(cell(sz, numel(opts.pheno)));
    for i = 1:sz
        for j = 1:numel(opts.pheno)
            skat{i, j} = opts.skat{i}(idx_skat(j):idx_skat(j+1)-1, :);
            if any(contains(lower(sheets), "overall"))
                overall{i, j} = opts.overall{i}(idx_overall(j):idx_overall(j+1)-1, :);
            end
            collapse{i, j} = opts.collapse{i}(idx_collapse(j):idx_collapse(j+1)-1, :);
        end
    end
    opts.skat = skat; opts.collapse = collapse;
    if any(contains(lower(sheets), "overall"))
        opts.overall = overall;
    end
end

end