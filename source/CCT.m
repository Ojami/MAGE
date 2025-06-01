function pval = CCT(pvals, weights, log10p, ignore0s, ignore1s)
    % CCT - Cauchy distribution test
    % A direct re-implementation of https://github.com/seanjosephjurgens/UKBB_200KWES_CVD
    %
    %   pval = CCT(pvals, weights, log10p, ignore0s, ignore1s)
    %
    %   pvals: vector of p-values.
    %   weights: vector of weights (default: equal weights).
    %   log10p: if true, pvals are given as -log10(p) (default: false).
    %   ignore0s: if true, p-values exactly equal to 0 are ignored (default: false).
    %   ignore1s: if true, p-values exactly equal to 1 are ignored (default: false).
    
    % Set defaults if needed
    if nargin < 5, ignore1s = false; end
    if nargin < 4, ignore0s = false; end
    if nargin < 3, log10p = false; end
    if nargin < 2, weights = []; end
    
    % If p-values are provided in -log10 scale, convert them.
    if log10p
        pvals = 10.^(-pvals);
    end
    
    % Remove NaN values from pvals (and corresponding weights, if provided)
    nanIdx = isnan(pvals);
    if any(nanIdx)
        pvals = pvals(~nanIdx);
        if ~isempty(weights)
            weights = weights(~nanIdx);
        end
    end
    if any(isnan(pvals))
        error('ERROR: could not evaluate p-value entries properly, NaNs remain in the data.');
    end
    
    % Check if all p-values are between 0 and 1.
    if any(pvals < 0 | pvals > 1)
        error('All p-values must be between 0 and 1!');
    end
    
    % Identify p-values exactly 0 or 1.
    hasZero = any(pvals == 0);
    hasOne  = any(pvals == 1);
    
    if hasZero && hasOne && ~(ignore0s || ignore1s)
        error('Cannot have both 0 and 1 p-values!');
    end
    
    % Handle p-values that are 0.
    if hasZero
        if ignore0s
            idx = (pvals == 0);
            pvals(idx) = [];
            if ~isempty(weights)
                weights(idx) = [];
            end
        else
            pval = 0;
            return;
        end
    end
    
    % Handle p-values that are 1.
    if hasOne
        if ignore1s
            warning('There are p-values that are exactly 1! Ignoring those...');
            idx = (pvals == 1);
            pvals(idx) = [];
            if ~isempty(weights)
                weights(idx) = [];
            end
        else
            warning('There are p-values that are exactly 1!');
            pval = 1;
            return;
        end
    end
    
    % Validate and set weights.
    if isempty(weights)
        weights = ones(size(pvals)) / numel(pvals);
    elseif numel(weights) ~= numel(pvals)
        error('The length of weights should be the same as that of the p-values!');
    elseif any(weights < 0)
        error('All the weights must be positive!');
    else
        weights = weights / sum(weights);
    end
    
    % Compute the CCT statistic.
    smallIdx = pvals < 1e-16;
    if ~any(smallIdx)
        cct_stat = sum(weights .* tan((0.5 - pvals) * pi));
    else
        cct_stat = sum((weights(smallIdx) ./ pvals(smallIdx))) / pi + ...
                   sum(weights(~smallIdx) .* tan((0.5 - pvals(~smallIdx)) * pi));
    end
    
    % Compute the p-value based on the CCT statistic.
    if cct_stat > 1e15
        pval = (1 / cct_stat) / pi;
    else
        % Standard Cauchy CDF: 0.5 + (1/pi)*atan(x)
        pval = 1 - (0.5 + (1/pi) * atan(cct_stat));
    end
    
    % If log10p flag was set, convert the result.
    if log10p
        pval = -log10(pval);
    end
end
