function chi2 = pval2chi(p, nu)
% p-value 2 chi2 when p-value is close to 0 (inverse survival function). 
% this is the same as chi2inv(1-x, 1) but more accurate for extreme
% values
% this is equivalent to chi2.isf: https://docs.scipy.org/doc/scipy-0.7.x/reference/generated/scipy.stats.chi2.html
 
chi2 =  gammaincinv(p, nu/2, 'upper') .* 2;
end