function out = susie(preds, response, opts)
% runs susie_auto to identify putative causal predictors.
% @09NOV2022

arguments
    preds {mustBeNonempty,mustBeA(preds, ["table", "double"])} % predictors
    response {mustBeNonempty, mustBeA(response, ["table", "double"])} % response
    opts.covars {mustBeNonempty,mustBeA(opts.covars, ["table", "double"])} % covariates to be adjusted for (can be left empty)
    opts.single (1,1) logical = false % convert preds to single (if preds is a large double matrix)
end

% tt = tic;
pred_names = "X" + (1:size(preds, 2));
if istable(preds)
    pred_names = string(preds.Properties.VariableNames);
    preds = preds{:, :};
end

if opts.single, preds = single(preds); end

% response_name = "Y";
if istable(response)
%     response_name = string(response.Properties.VariableNames(1));
    response = response{:, 1};
end

nan_idx = any(isnan(preds), 2) | any(isnan(response), 2);

if isfield(opts, 'covars')
    if istable(opts.covars)
        opts.covars = opts.covars{:, :};
    end
    nan_idx = nan_idx | any(isnan(opts.covars), 2);
    opts.covars(nan_idx, :) = [];
end
preds(nan_idx, :) = [];
response(nan_idx, :) = [];

rcode = string;
rcode(10, 1) = "X <- R.matlab::readMat('susie.raw.mat')";
rcode(11, 1) = "Y <- X$response";
rcode(13, 1) = "X <- X$preds";

if isfield(opts, 'covars')
    covars = opts.covars;
    save('susie.raw.mat', 'preds', 'response', 'covars');
    rcode(1, 1) = "remove.covariate.effects <- function (X, Z, y) {";
    rcode(2, 1) = "  # include the intercept term";
    rcode(3, 1) = "  if (any(Z[,1]!=1)) Z = cbind(1, Z)";
    rcode(4, 1) = "  A   <- Matrix::forceSymmetric(crossprod(Z))";
    rcode(5, 1) = "  SZy <- as.vector(solve(A,c(y %*% Z)))";
    rcode(6, 1) = "  SZX <- as.matrix(solve(A,t(Z) %*% X))";
    rcode(7, 1) = "  y <- y - c(Z %*% SZy)";
    rcode(8, 1) = "  X <- X - Z %*% SZX";
    rcode(9, 1) = "  return(list(X = X,y = y,SZy = SZy,SZX = SZX))}";
    rcode(12, 1) = "Z <- X$covars";
    rcode(16, 1) = "out = remove.covariate.effects(X, Z, Y[,1])";
    rcode(18, 1) = "res <- susie_auto(out$X, out$y)"; % automatic version by a 3-step strategy
else
    rcode(18, 1) = "res <- susie_auto(X, Y)"; % automatic version by a 3-step strategy
    save('susie.raw.mat', 'preds', 'response');
end

rcode(17, 1) = "library(susieR)";
rcode(19, 1) = "beta <- coef(res)";
rcode(20, 1) = "pip <- susie_get_pip(res)";
rcode(21, 1) = "write.table(cbind(beta, c(0,pip)), 'susie.beta.txt', quote = F, row.names = F, col.names = T)";
writematrix(rcode, "susie_featureselection.r", 'QuoteStrings', false, 'FileType','text')
MATLAB2Rconnector("susie_featureselection.r", 'delr', true);
delete('susie.raw.mat');
out = readtable('susie.beta.txt');
out.Properties.VariableNames = {'beta', 'pip'};
out.name = ["Intercept"; pred_names'];
delete('susie.beta.txt')
% fprintf('\n\b\b done (%.1f sec)\n', toc(tt)) 

end % END