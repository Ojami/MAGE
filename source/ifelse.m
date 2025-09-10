function out = ifelse(cond, yesVal, noVal)
%IFELSE  Element-wise conditional assignment (R-style)
%
%   out = ifelse(cond, yesVal, noVal)
%
%   - cond: logical array
%   - yesVal: value(s) used where cond is true
%   - noVal: value(s) used where cond is false
%   - All inputs must be scalar or same size as cond

    % Validate condition
    if ~islogical(cond)
        error('Condition must be a logical array.');
    end

    % Get output size
    sz = size(cond);

    % Handle scalar expansion
    if isscalar(yesVal)
        yesVal = repmat(yesVal, sz);
    end

    if isscalar(noVal)
        noVal = repmat(noVal, sz);
    end

    % Check size consistency
    if ~isequal(size(yesVal), sz)
        error('yesVal must be scalar or same size as condition.');
    end

    if ~isequal(size(noVal), sz)
        error('noVal must be scalar or same size as condition.');
    end

    % Combine values
    out = yesVal;
    out(~cond) = noVal(~cond);
end
