function str = smartFloat(x, n)
% floatToSmartString  Returns a string representing x with n decimal places.
%   If |x| is in [10^(-n), 10^(n+1)), it uses fixed-point notation (%.nf).
%   Otherwise, it uses scientific notation (%.nE).
%
%   USAGE:
%       str = floatToSmartString(x, n)
%
%   INPUTS:
%       x : numeric scalar
%       n : integer specifying the number of decimal places
%
%   OUTPUT:
%       str : string with appropriately formatted numeric value

    % Define magnitude boundaries
    lowerBound = 10^(-n);
    upperBound = 10^(n+1);

    % Build format specifiers for fixed and scientific at n decimals
    formatFixed = compose("%%.%df", n);  % e.g. "%.2f"
    formatSci   = compose("%%.%dE", n);  % e.g. "%.2E"
    
    str = strings(numel(x), 1);

    idx = (abs(x) >= lowerBound) & (abs(x) < upperBound);
    str(idx) = compose(formatFixed, x(idx));
    str(~idx) = compose(formatSci, x(~idx));

    % for k = 1:numel(x)
    %     % Decide which format to use based on x’s absolute value
    %     if (abs(x(k)) >= lowerBound) && (abs(x(k)) < upperBound)
    %         str(k) = compose(formatFixed, x(k));
    %     else
    %         str(k) = compose(formatSci, x(k));
    %     end
    % end

end % END
