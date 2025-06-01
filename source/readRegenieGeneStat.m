function readRegenieGeneStat(file)

arguments
    file {mustBeFile}
end

% read masks
fid = fopen(file);
masks = fgetl(fid);
fclose(fid);
masks = extractBetween(masks, textBoundary("start") + "##MASKS=<", ...
    ">" + textBoundary("end"));
masks = regexp(string(masks), '(.*?)[=](.*?)(?:;|$)', 'tokens');
masks = vertcat(masks{:});
masks = erase(masks, [textBoundary('start') + '"', '"' + textBoundary('end')]);
masks = array2table(masks, 'VariableNames', {'mask', 'category'});

% read summary stat file
ss = readtable(file, 'TextType', 'string', 'VariableNamingRule', ...
    'preserve', 'CommentStyle', '##', 'ReadVariableNames', true);
if height(ss) > 1 % otherwise, nothing to read
    
    ss.Properties.UserData = masks;
    ss = convertvars(ss, ["BETA", "SE", "CHISQ", "LOG10P", "N", "A1FREQ", ...
        "GENPOS", "CHROM"], @double);
    ss.ID = extractBefore(ss.ID, ("."|textBoundary("end")));
    ss.ALLELE1(ss.ALLELE1 == "NA") = "Joint"; % joint test for burden masks
    
    % group by different tests, AAF bins and masks
    maskcats = unique(ss.ALLELE1);
    res = ({}); ct = 1; schema = string;
    for i = 1:numel(maskcats)
        tab = ss(ss.ALLELE1 == maskcats(i), :);
        test = unique(tab.TEST);
        for j = 1:numel(test)
            idx = tab.TEST == test(j);
            res{ct, 1} = tab(idx, :);
            schema(ct, 1) = maskcats(i) + "." + test(j);
            res{ct}.Properties.Description = maskcats(i) + "." + test(j);
            res{ct}.P = 10.^-res{ct}.LOG10P;
            res{ct} = sortrows(res{ct}, 'P', 'ascend');
            fprintf('%d-%s, min P = %.3g\n', ct, maskcats(i) + "." + test(j), res{ct}.P(1))
            ct = ct + 1;
        end
    end

else % nothing to read
    res{1} = table;
    schema = missing;
end

save(regexprep(file, ".txt$", ".mat"), 'res', 'schema')

end