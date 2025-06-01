function [p, X2] = chisqtest(O, df)
% INPUT:
%       O    : observed frequencies.
% OUTPUT:
%       p : chi-square p-value.
%       X2: chi-square statistics.
% Oveis Jamialahmadi, Sept. 2020. GU

if nargin < 2
    df = (size(O, 1) - 1)*(size(O, 2) - 1);
end
E = sum(O, 2)*sum(O, 1)/sum(O(:));
X2 = sum(sum((O-E).^2./E));
p = chi2pval(X2, df);

end