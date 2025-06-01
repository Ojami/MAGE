function [p, x2] = hwe(O)
% calculates Hardy-Weinberg equilibrium test
% INPUT: 
%       O: observed frequencies in the order of [hom.1 count, het count, 
%       hom.2 count]. 
% OUTPUT:
%       p: HWE p-value
%       x2: chi-square stat.
% 
% Oveis Jamialahmadi, Sep. 2020, GU

if nargin < 1
    help hwe
    return
elseif numel(O) ~= 3
    fprintf('ERROR: observed genotype frequencies must be a 3X1 numeric vector!\n')
    return
end

N = sum(O);
q = (2*O(1) + O(2)) / (N*2);
E = [N*q^2, 2*N*q*(1-q), N*(1-q)^2];
x2 = sum(sum((O-E).^2./E));

p = chi2pval(x2, 1);
end %END