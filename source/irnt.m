function valout = irnt(valin)
% perform Inverse Normal Transformation (IRNT)
% Oveis Jamialahmadi, Sep. 2020. GU.

if nargin < 1
    error('no input?')
elseif ~isrow(valin) && ~isvector(valin)
    error('row/column vector only!')
end

idx = ~isnan(valin);
valout = nan(size(valin));

% check nomrality 
if ~kstest(valin(idx))
    warning('input data are already normal!')
end
valout(idx) = norminv((tiedrank(valin(idx))-0.5)./sum(idx));
end % END