function [pval, tstat] = lrtest(nestedModel, complexModel, verbose)
% calculates LRT test statistics and p-value with 1 df
% Oveis Jamialahmadi, 15APR2024. University of Gothenburg.

if nargin < 3
    verbose = true;
end

df_nested = nestedModel.NumPredictors + 1;
df_complex = complexModel.NumPredictors + 1;
df = abs(df_nested - df_complex);
A = nestedModel.LogLikelihood;
B = complexModel.LogLikelihood;
tstat = abs(-2*(A - B));
pval = chi2pval(tstat, df);

% prepare output similar to lrtest R function
r = string;
r(1) = "Model 1: " + string(nestedModel.Formula);
r(2) = "Model 2: " + string(complexModel.Formula);
r(3) = sprintf("#DF\tLogLik\tDf\tChisq\tPr(>Chisq)");
r(4) = sprintf("%d\t%.1f", df_nested, A);                       
r(5) = sprintf("%d\t%.1f\t%d\t%.2f\t%.2g", df_complex, B, df, tstat, pval);
r = r';
if verbose
    for k = 1:numel(r)
        disp(r(k))
    end
end

end % END