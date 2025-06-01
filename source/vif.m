function result = vif(model)
% calculates variance inflation factor (VIF). VIF > 5 can be considered as
% collinearity cutoff.
% an implementation of https://github.com/cran/car/blob/master/R/vif.R

result = NaN;
v = model.CoefficientCovariance;
assign = model.CoefficientNames;

if strcmp(model.CoefficientNames(1), '(Intercept)')
    v = v(2:end, 2:end);
    assign(1) = [];
else
    fprintf('WARNING: No intercept: vifs may not be sensible.\n')
end

terms = model.PredictorNames;
if numel(terms) < 2
    fprintf('ERROR: model contains fewer than 2 terms\n')
    return
end

R = corrcov(v);
detR = det(R);

result  = splitvars(table(zeros(numel(terms), 3)));
result.Properties.RowNames = terms;
result.Properties.VariableNames = {'GVIF', 'Df', 'GVIF^(1/(2*Df)'};

for i = 1:numel(terms)
    subs = contains(assign, terms(i));
    result.(1)(i) = det(R(subs, subs)) * det(R(~subs, ~subs)) / detR;
    result.(2)(i) = sum(subs);
end

if all(result.(2) == 1)
    result.(3) = result.(1);
else
    result.(3) = result.(1).^(1./(2.*result.(2)));
end

end %END

%% 
% v = (0);
% if istable(D)
%     cols = D.Properties.VariableNames;
%     D = D.Variables;
% else
%     cols = [];
% end
% 
% for i = 1:size(D,2)
%     model = fitlm(D(:,setdiff(1:size(D,2),i)),D(:,i));
%     v(i) = model.Rsquared.Ordinary;
% end
% v = 1./(1-v);
% if ~isempty(cols)
%     v = splitvars(table([string(cols)', v']));
%     v.Properties.VariableNames = {'varName', 'VIF'};
% end