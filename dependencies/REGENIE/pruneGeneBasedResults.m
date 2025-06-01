function [res, leg] = pruneGeneBasedResults(tab, opts)
% prunes gene-based summary stats from REGENIE. This is a subfunction of
% 'gwasrunner' function.

arguments
    tab {mustBeA(tab, 'table')}
    opts.output {mustBeTextScalar} % output file
end

idx = ismember(colnames(tab), ["GENPOS", "CHISQ", "LOG10P"]);
tab(:, idx) = [];

if any(colnames(tab) == "ALLELE1")
    tab.A1 = tab.ALLELE1 + "." + tab.TEST;
    tab(:, ["ALLELE0", "ALLELE1", "TEST"]) = [];
end

% check if Pheno column is present
if ~any(colnames(tab).lower == "pheno")
    tab.Pheno = repmat("NA", height(tab), 1);
    tab = movevars(tab, "Pheno", 'Before', 1);
end

% check ADD/DOM/REC test
adc = ["ADD", "DOM", "REC"];
adcIdx = arrayfun(@(x)contains(tab.A1, x), "." + adc, uni=false);
adc = adc(cellfun(@any, adcIdx));

tests = extract(tab.A1, "." + adc + wildcardPattern + textBoundary("end"));
tab.Mask = strrep(tab.A1, tests, "");
pheno = unique(tab.Pheno);
mask = unique(tab.Mask);
test = adc + ["-BURDEN-MINP", "-BURDEN-SBAT", "-BURDEN-ACAT", "",...
    "-ACATO", "-ACATV", "-SKAT", "-SKATO",...
    "-SKATO-ACAT", "-ACATO-FULL"];% for reformatting the table

% maf groups (+singleton and Joint)
% maf = "0.01"; % maf cutoff to filter
maf = unique(extractAfter(mask, "."));
maf(ismissing(maf)) = [];
joint_mask = find(mask.lower == "joint", 1);
if ~isempty(joint_mask), maf = union(maf, mask(joint_mask)); end

mask(~endsWith(mask, maf)) = [];
selcols = string(union(setdiff(tab.Properties.VariableNames, ["A1", "P"], "stable"), "GO", "stable"));
sortedCols = ["Pheno", "ID", "CHR", "Mask", "BETA", "SE", test, "N"];
for j = 1:numel(maf)
    out = ({}); ct = 1;
    tmpout = tab;  
    tmpout(~endsWith(tmpout.Mask, maf(j)), :) = [];

    for i = 1:numel(pheno)
        tmp1 = tmpout;
        tmp1(tmp1.Pheno ~= pheno(i), :) = [];
        if isempty(tmp1), continue; end
        for k = 1:numel(mask)
            tmp = tmp1;
            tmp(tmp.Mask ~= mask(k), :) = [];
            if isempty(tmp), continue; end
            
            uID = unique(tmp.ID); 
            for l = 1:numel(uID) % loop over genes
                tmp2 = tmp(tmp.ID == uID(l), :);
                tmp2 = sortrows(tmp2, "SE");
                A1 = strrep(tmp2.A1, tmp2.Mask + ".", "");
                P = tmp2.P;
                tmp3 = nan(1, numel(test));
                [f1, f2] = ismember(test, A1);
                tmp3(1, f1) = P(f2(f1));
                
                cols = string(tmp2.Properties.VariableNames);
                if ~any(cols.lower.startsWith("go"))
                    selcols(selcols == "GO") = [];
                end
                tmp2 = tmp2(1, selcols);
                tmp2{:, test} = tmp3;
                out{ct} = tmp2; ct = ct + 1;
            end
        end 
    end
    
    out(cellfun(@isempty, out)) = [];
    out = vertcat(out{:});
    out = rmmissing(out, 2, "MinNumMissing", height(out));
    cols = string(out.Properties.VariableNames);
    beta_idx = contains(cols.lower, ["or|", "β|"]) | ismember(cols.lower, ["β", "or"]);
    scols = sortedCols;
    if any(beta_idx)
        scols(scols == "BETA") = cols(beta_idx);
    else
        scols(scols == "BETA") = [];
    end
    [f1,f2] = ismember(scols, cols);
    f3 = ~ismember(cols, scols);
    scols = [cols(f2(f1)), setdiff(cols(f3), scols, "stable")];

    res.("x"+j) = out(:, scols);
end

res.dict = dictionary("x"+(1:numel(maf)), maf');

% add legend
% find different masks and tests
ucats = unique(tab.A1);
lastDot = regexp(ucats, '[.]','end');
leg.test = strings(numel(ucats), 1);
leg.mask = extractBefore(ucats, ".");
for i = 1:numel(lastDot)
    pos = lastDot{i}(end);
    tmp = ucats(i).char;
    leg.test(i) = string(tmp(pos+ 1:end));
end
leg.maf = replace(ucats, leg.mask + ".", "");
leg.maf = extractBefore(leg.maf, "." + leg.test);

fis = string(fieldnames(leg));
for i = 1:numel(fis)
    leg.(fis(i))(ismissing(leg.(fis(i))) | leg.(fis(i)) == "") = [];
    leg.(fis(i)) = unique(leg.(fis(i)));
end

leg.mask = erase(leg.mask, "mask_");
leg.mask(leg.mask.lower == "joint") = [];
tab = array2table(strings(numel(leg.mask), 5), "VariableNames", ...
    ["cadd", "revel", "polyphen", "sift", "clinvar"]);
cols = string(tab.Properties.VariableNames);
for i = 1:numel(cols)
    switch cols(i)
        case "cadd", tag = "[c](\d+)";
        case "revel", tag = "[r](\d+)";
        case "polyphen", tag = "[^\w]p(?<n>\w{1})([^\w]|$)";
        case "sift", tag = "[^\w]s(?<n>\w{1})([^\w]|$)";
        case "clinvar", tag = "(cvar)";
    end
    
    if any(ismember(cols(i), ["polyphen", "sift"]))
        idx = regexp(leg.mask, tag, 'names');
    else
        idx = regexp(leg.mask, tag, 'tokens');
    end

    if isempty(idx), continue; end % a custom mask
    rmidx = cellfun(@isempty, idx);
    idx(rmidx) = []; idx = vertcat(idx{:});
    if isstruct(idx), idx = struct2array(idx); end
    tab.(cols(i))(~rmidx) = string(idx);
    if cols(i) == "revel"
        tab.(cols(i)) = string(double(tab.(cols(i)))/10);
    elseif any(ismember(cols(i), ["polyphen", "sift"]))
        tab.(cols(i)) = replace(tab.(cols(i)), ["r", "s"], ["Relaxed", "Strict"]);
    elseif cols(i) == "clinvar"
        tab.(cols(i)) = replace(tab.(cols(i)), "cvar", "ClinVar");
    end
    
    if any(ismember(cols(i), ["cadd", "revel"]))
        tab.(cols(i))(~rmidx) = "≥" + tab.(cols(i))(~rmidx);
    elseif any(ismember(cols(i), ["polyphen", "sift"]))
        tab.(cols(i))(~rmidx) = "=" + tab.(cols(i))(~rmidx);
    end
end

tab = standardizeMissing(tab, "");
tab.clinvar(~ismissing(tab.clinvar)) = "";
cols = cols.upper;
cols = replace(cols, ["POLYPHEN", "CLINVAR"], ["PolyPhen", "ClinVar"]);
for i = 1:height(tab)
    rule = "";
    if contains(leg.mask(i).lower, "lof")
        rule = "LoF or ";
    end

    if contains(leg.mask(i).lower, "missense")
        rule = rule + "missnese or ";
    end
    
    if contains(leg.mask(i), "&")
        glue = " and ";
    else
        glue = " or ";
    end
    
    idx = ~ismissing(tab{i, :});
    glued = join(cols(idx) + tab{i, idx}, glue);
    glued(ismissing(glued)) = [];

    if ~ismissing(glued)
        if any(contains(glued, glue))
            rule = rule + "(" + glued + ")";
        else
            rule = rule + glued;
        end
    else
        rule = regexprep(rule, " or $", "");
    end

    tab.term(i) = rule;
end
leg.mask = table(leg.mask, tab.term, 'VariableNames', {'Mask', 'Definition'});

leg.test = table(["ADD", "SKAT", "SKATO", "SKATO-ACAT", "ACATV", "ACATO", ...
    "ACATO-FULL"]', ...
    ["Burden test"
    "Variance component test"
    "Omnibus test combining features of SKAT and Burden"
    "Same as SKATO but using Cauchy combination method to maximize power across SKATO models"
    "Test using Cauchy combination method to combine single-variant p-values"
    "Omnibus test combining features of ACATV, SKAT and Burden"
    "Same as ACATO but using the larger set of SKATO models used in the SKATO test"], ...
    'VariableNames', {'Test', 'Desciption'});

sheetrange1 = height(leg.test) + 3;
leg.mask.Definition(leg.mask.Definition == "") = "NA";
sheetrange2 = sheetrange1 + height(leg.mask) + 2;
leg.maf = ["Alternate AF bins: ", join(leg.maf, ",")];

leg.dict = dictionary(["test", "mask", "maf"], "A" + [1, sheetrange1, sheetrange2]);

% write to output
if ~isfield(opts, 'output')
    return
else
    opts.output = regexprep(opts.output, [".xlsx$", ".xls$"], "");
    pth = fileparts(opts.output);
    if pth == ""
        pth = pwd;
    end
    if ~isfolder(pth), mkdir(pth); end
end

for i = 1:res.dict.numEntries
    writetable(res.("x"+i), opts.output + ".xlsx", ...
        Sheet=res.dict("x"+i), AutoFitWidth=true, ...
        WriteMode="overwritesheet");
end

% write legend
legkeys = leg.dict.keys;
for i = 1:leg.dict.numEntries
    try
        writetable(leg.(legkeys(i)), opts.output + ".xlsx", ...
            Sheet="Legend", AutoFitWidth=true, Range=leg.dict(legkeys(i)));
    catch
        writematrix(leg.(legkeys(i)), opts.output + ".xlsx", ...
            Sheet="Legend", Range=leg.dict(legkeys(i)));
    end
end

end