function res = delongtest(df1, df2, opts)
% a wrapper for DeLong's test for two ROC curves using pROC package.
% Oveis Jamialahmadi, University of Gothenburg, March 2023

arguments
    df1 {mustBeA(df1, 'table')} % first df
    df2 {mustBeA(df2, 'table')} % second df
    opts.test {mustBeMember(opts.test, ["bootstrap", "delong"])} = "delong"
end

% last column is response variable
patt = getRandomName("df", 5);
out = fullfile(pwd, patt + ".txt");
out2 = fullfile(pwd, patt + "2.txt");
out1 = fullfile(pwd, patt + "1.txt");

writetable(df1, out1)
writetable(df2, out2)

c1 = colnames(df1);
c2 = colnames(df2);
formula1 = c1(end) + "~" + join(c1(1:end-1), "+");
formula2 = c2(end) + "~" + join(c2(1:end-1), "+");

r = string;
r(1) = "library(pROC)";
r(2) = 'df1 = read.table("' + out1.replace(filesep, "/") + '", sep = ",", header = T)';
r(3) = 'df2 = read.table("' + out2.replace(filesep, "/") + '", sep = ",", header = T)';
r(4) = "mdl1 = glm(" + formula1 + ", data = df1, family = binomial())";
r(5) = "mdl2 = glm(" + formula2 + ", data = df2, family = binomial())";
r(6) = 'prob1 = predict(mdl1, df1, type = "response")';
r(7) = 'prob2 = predict(mdl2, df2, type = "response")';
r(8) = "roc1 = roc(df1$" + c1(end) + " ~ prob1)";
r(9) = "roc2 = roc(df2$" + c2(end) + " ~ prob2)";
r(10) = 'out = roc.test(roc1, roc2, method = "' + opts.test + '")';
r(11) = "out = c(out$estimate, out$statistic, out$p.value)";
r(12) = "names(out)[length(out)] = 'P'";
r(13) = 'write.table(t(out), "' + out.replace(filesep, "/") + ...
    '", quote = F, row.names = F, col.names = T, sep = "||")';
MATLAB2Rconnector(patt, code=r, delr=true);
res = readtable(out, Delimiter="||", ...
    VariableNamingRule="preserve", ...
    ReadVariableNames=true, NumHeaderLines=0);

delete(out1)
delete(out2)
delete(out)

end