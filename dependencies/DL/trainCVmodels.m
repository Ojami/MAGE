function [beta, stat, data] = trainCVmodels(data, opts)
% performs nested cross validation. After splitting the original dataset
% into training/test using hold-out CV (default is 20/80 for test/train),
% nested CV is performed on the training set as follows: outer loop splits
% original training set into new training/validation set, and each training
% set is used to find tuned hyperparemters (if any) in the inner loop CV.
% 
% References: 
%       - the final model is trained on whole dataset ('data' table), with
%         the hyperparameters tuned on training set and performance
%         reported on test set:
%         https://machinelearningmastery.com/train-final-machine-learning-model/
%         https://datascience.stackexchange.com/questions/33008/is-it-always-better-to-use-the-whole-dataset-to-train-the-final-model
%         https://stats.stackexchange.com/questions/184095/should-final-production-ready-model-be-trained-on-complete-data-or-just-on-tra
%         https://sebastianraschka.com/blog/2016/model-evaluation-selection-part3.html
% 
%       - nested cross-validation:
%         https://stats.stackexchange.com/questions/553658/cross-validation-and-hyperparameter-tuning-workflow
%         https://sebastianraschka.com/blog/2018/model-evaluation-selection-part4.html
%         https://ploomber.io/blog/nested-cv/
%         https://www.philipscurve.com/post/2020-10-24-nested-cross-validation/#fnref3
%         (which supports to tune hyperparameters after finding the best model strategy)

% 
% OUTPUTs
%   beta    if ml == "penalized": 
%           a cell array a tables with each table correponding to methods
%           in each row  of stat.best, respectively. Note that beta for
%           each method is a table with three columns: feature names,
%           train, and full. "train" shows the coefficients training set,
%           while "full" shows the coefficients on full dataset (train +
%           test) with the hyperparameters tuned on the training set. 
% 
%           if ml == "overall":
%           best model from the set of machine-learning methods set in
%           'tool', trained on the training set a,d full dataset (train +
%           test)
% 
%   stat    a struct with fields:
%               -mse: performance of each model on whole training set in a
%                      nested CV. AICs in this table cannot be trusted
%                      since they correspond to each fold of CV (outer CV).
%               -aic: performance of each model on whole training set. Used
%                     to find the best model based on lowest AIC. note that
%                     AIC (on whole training set) and CV MSE are both
%                     asymptotically equivalent.
%               -best: performance on test set. This should be formally
%                      reported. Note that either of methods (AIC- or
%                      MSE-based) should be kept in the end (more common:
%                      MSE)
%               -shap: if 'shap' is set to true SHAP values will be
%                      calculated.
% 
%   data    table of original data with proper transformation/normalization
%           on target/features. This can be useful when original features
%           should undergo transformation and full model is used for
%           deployment or other purposes (e.g. ICE or PDP plots).
% 
% Example:
% [beta,stat] = trainCVmodels(data(:, 2:end), "ml",...
%     "penalized", "kfold", 5, "featureNames", colnames(data(:,2:end-1)),...
%     "ignore_aic",true,"parallel",true, "targetTransform","none", ...
%     "tool", [ "elnet", "ridge", "lasso", "plsr"], "target", "PDFF",...
%     "shap", true);
% 
% @15JUNE2023: 'CategoricalPredictors' option was added. This is only
%              supported when ml is 'overall'. If left empty, automatically
%              checks the categorical fatures.
%              Support for visualizing SHAP values and ICE plots has been
%              added.
% 
% @10AUGUST2023: 'shap_method' was added and can be selected as either
%               "interventional" or "conditional" (AKA observational). In
%               case of correlated features, "conditional" is preferred.
%               Notes from MATLAB doc: This approach uses more realistic
%               samples than an interventional algorithm and does not
%               require the feature independence assumption. However, a
%               conditional algorithm is computationally more expensive,
%               does not support ordered categorical predictors, and cannot
%               handle NaNs in continuous features. Also, the algorithm
%               might assign a nonzero Shapley value to a dummy feature,
%               which does not contribute to the prediction, if the dummy
%               feature is correlated with an important feature. See more
%               details in this post:
%               https://pacmedhealth.medium.com/explainability-for-tree-based-models-which-shap-approximation-is-best-6df78bc5d086
% 
% @12SEP2024: SHAP calculation method has been updated (faster in R2024
%             onwards). Also plots summary (swarm) plots with additional
%             colormaps. 
% @13SEP2024: 'toolset' argument was added to use nested CV for a subset of
%             model used in "auto_sep". The methods used in 'toolset' are a
%             subset of models defined in 'auto_sep'. This is useful when
%             use wants to compare only a subset of methods. 'toolset'
%             overrides 'tool' argument.

arguments
    data {mustBeA(data, "table")}
    opts.target {mustBeTextScalar} % by default the last variable is used as the target variable
    opts.kfold (1,1) double = 10 % inner CV
    opts.outerkfold (1,1) double % outer CV, by default is the same as 'kfold'
    opts.holdout (1,1) double {mustBeInRange(opts.holdout, 0.001, 0.999)} = 0.2 % test set size

    % to standardize input features. Note that 'standardization' is always
    % done for penalized regressions, and has a different purpose here. If
    % set to true (default) input features are standardized (mean 0, sd 1,
    % or centered and scaled), and output coefficients are standardized. If
    % false, still the returned coefficients will be returned to their
    % original scale (not recommended if one wants to compare the features'
    % coefficients directly).
    opts.standardize (1,1) logical = true
    opts.targetTransform {mustBeTextScalar, mustBeMember(opts.targetTransform, ["irnt", "none", "log", "log+1", "rescale", "zscore"])} = "none" % transformation on the target response
    opts.parallel (1,1) logical = true 
    opts.elnet {mustBeA(opts.elnet, "struct")} % a struct with lambda and alpha fields to define gridsize for elnet, ridge and lasso models
    opts.featureNames {mustBeText, mustBeVector} % feature names, if empty, input dataset variable names will be used
    opts.ignore_aic (1,1) logical = false % don't compare models based on AIC (on whole training set)
    opts.ml {mustBeMember(opts.ml, ["penalized", "overall"])} = "penalized" % overall: a wider range of ML methods can be used (set with 'tool')
    opts.tool {mustBeMember(opts.tool, ["auto", "nnet_auto", "lm", ...
        "lm_interaction", "gpr_exponential","gpr_rationalquadratic",...
        "gpr_squaredexponential", "gpr_matern52", "tree_fine",...
        "tree_medium", "tree_coarse", "ensemble_boosted", ...
        "ensemble_bagged", "ensemble_auto", "svm_linear", "svm_quadratic", ...
        "svm_cubic", "svm_finegaussian", "svm_mediumgaussian", ...
        "svm_coarsegaussian", "svm_auto", "gam", "auto_sep", ...
        "nnet_narrow", "nnet_medium", "nnet_wide", "nnet_bi", "nnet_tri", ...
        "gam_auto", "susie", "elnet", "ridge", "lasso", "plsr"])} = "auto_sep"
    opts.toolset {mustBeMember(opts.toolset, ["net", "ensemble", ...
        "svm", "kernel", "linear", "tree"])}

    opts.shap (1,1) logical = false % to return SHAP values
    opts.shap_subsample (1,1) double = 1e3 % max number of sumsamples to be drawn randomly for computational tractability.
    opts.shap_method {mustBeMember(opts.shap_method, ["interventional", "conditional"])} = "interventional"

    opts.covariates {mustBeA(opts.covariates, 'table')} % covariates to be adjusted for their effects. They'll be regressed out.
    opts.verbose (1,1) logical = true
    opts.CategoricalPredictors {mustBeText, mustBeVector}

    % plotting options
    opts.plotshap (1,1) logical = false % only if 'shap' is true
    opts.plotice (1,1) logical = false % only if ml is 'overall'
    opts.NumObservationsToSample (1,1) double = 1000 % for 'plotice'
    opts.markersize (1,1) double = 15 % marker size for scatter plot in 'plotshap' 
    opts.WindowState {mustBeTextScalar, mustBeMember(opts.WindowState, ["maximized", "normal"])} = "maximized"
    opts.WindowPosition (1, 4) double = [1, 41, 1280, 720] % windows position in pixels, overrides WindowState
    opts.TileSpacing {mustBeTextScalar, mustBeMember(opts.TileSpacing, ["tight", "loose", "compact"])} = "tight"
    opts.grid (1,1) logical = true
    opts.box (1,1) logical = true
    opts.colormap {mustBeMember(opts.colormap, ["jet", "turbo", ...
        "parula", "hot", "winter", "copper", "bone", "hsv", ...
        "summer", "shap", "magma", "plasma", "inferno", "cividis", ...
        "viridis"])} = "plasma" % colormap for 'plotshap' (MATLAB or 'get_colormap' supported colormaps)
    opts.fontsize (1,1) double = 11
    opts.fontname {mustBeTextScalar, mustBeMember(opts.fontname, ["Bauhaus", "Futura", "Bodoni", "Garamond", "Corbel", "Corbel Light", "Arial", "Tahoma"])} = "Arial"
    opts.figname {mustBeTextScalar} = ""
    opts.format {mustBeMember(opts.format, ["jpg", "png", "tif", "eps", "pdf"])} = "png" % format of output plots
    opts.resolution (1, 1) double = 300 % resolution of output plot
    opts.save (1,1) logical = true
    opts.savefig (1,1) logical = false % Saves .fig files as well
    opts.rng (1,1) double = 0 % default MATLAB random generator seed

end

rng(opts.rng, "twister")

% reset plotting arguments 
if ~opts.shap, opts.plotshap = false; end
if opts.ml ~= "overall", opts.plotice = false; end

if ~isfield(opts, 'outerkfold'), opts.outerkfold = opts.kfold; end
if ~isfield(opts, 'target')
    opts.target = string(data.Properties.VariableNames(end));
end
data = movevars(data, opts.target, "After", width(data));
data = rmmissing(data);

if ~isfield(opts, 'featureNames')
    opts.featureNames = string(data.Properties.VariableNames(1:end-1));
end

if ~isfield(opts, 'CategoricalPredictors')
    opts.CategoricalPredictors = varfun(@iscategorical, data(:, 1:end-1));
    opts.CategoricalPredictors = opts.CategoricalPredictors{:, :}; % opts.featureNames(opts.CategoricalPredictors{:, :});
else
    opts.CategoricalPredictors = ismember(opts.featureNames, opts.CategoricalPredictors);
    % opts.CategoricalPredictors = false(numel(opts.featureNames), 1);
end

if opts.verbose
    disp('====================== summary ======================')
    fprintf('dataset has %d rows and %d features\n', height(data), width(data) - 1)
    fprintf("target variable: %s\n", string(opts.target)) 
    if opts.shap, fprintf("SHAP will be calculated for (max) %d samples\n", opts.shap_subsample); end
    
    if opts.ml == "penalized"
        fprintf("%d inner CV will be performed.\n", opts.kfold)
        fprintf("%d outer CV will be performed.\n", opts.outerkfold)
    else
        fprintf("%d CV will be performed.\n", opts.kfold)
        fprintf("\tthis number is the same for inner/outer CV\n")
    end

    if any(opts.CategoricalPredictors)
        fprintf("%d categorical features have been found.\n", sum(opts.CategoricalPredictors))
    end

end

if isrow(opts.featureNames), opts.featureNames = opts.featureNames'; end

opts.statset = statset('UseParallel', opts.parallel);

% define training/test sets
cv = cvpartition(height(data), "Holdout", opts.holdout);

x_test = data(cv.test, :); % for final model's performance assessment
x_train = data(cv.training, :); % for training
% clear data

% % for visualization pruposes
% opts.raw_x_train = x_train; 
% opts.raw_x_test = x_test;

if opts.verbose
    fprintf('created training (%d rows) and test (%d rows) sets\n', ...
        height(x_train), height(x_test))
end

if (opts.targetTransform ~= "none")
    if opts.targetTransform == "log+1"
        data{:, end} = log(data{:, end} + 1); % for full model fit
        x_train{:, end} = log(x_train{:, end} + 1);
        x_test{:, end} = log(x_test{:, end} + 1);
    else
        data{:, end} = feval(opts.targetTransform, data{:, end});
        x_train{:, end} = feval(opts.targetTransform, x_train{:, end});
        x_test{:, end} = feval(opts.targetTransform, x_test{:, end});
    end

    if opts.verbose
        fprintf('target variable has been transformed: %s\n', opts.targetTransform)
    end
end

% adjust for covariates (regress out from predictors/response)
if isfield(opts, "covariates")
    if opts.verbose
        fprintf('regression out %d covariates.\n', width(opts.covariates))
    end

    if any(opts.CategoricalPredictors)
        error(sprintf("\tcategorical features were found."+ ...
            "\n\ttry to include covariates as features in the input dataset.\n"))
    end
    
    [data(:, 1:end-1), data(:, end)] = removeCovarEffect(data(:, 1:end-1), data(:, end), opts.covariates);
    [x_train(:, 1:end-1), x_train(:, end)] = removeCovarEffect(x_train(:, 1:end-1), x_train(:, end), opts.covariates(cv.training, :));
    [x_test(:, 1:end-1), x_test(:, end)] = removeCovarEffect(x_test(:, 1:end-1), x_test(:, end), opts.covariates(cv.test, :));
end

if opts.standardize
    if opts.verbose
        disp('standardizing (mean 0, variance 1) features.')
    end
    
    if any(opts.CategoricalPredictors)
        fprintf("\tskipping categorical features.\n")
    end

    cidx = [~opts.CategoricalPredictors, false];
    data(:, cidx) = normalize(data(:,cidx));
    x_train(:, cidx) = normalize(x_train(:, cidx));
    x_test(:, cidx) = normalize(x_test(:, cidx));
end


if opts.verbose
    disp('=====================================================')
end

% initialize elnet params
if ~isfield(opts, 'elnet')
    opts.elnet.lambda = logspace(-5, 2, 300);
    opts.elnet.alpha = linspace(1e-6, 0.9999, 50); % for elastic net
end

%% ------------------------------------------------------------------------
if any(ismember(opts.ml, ["penalized", "all"]))
    % for penalized regression models including PLSR (not penalized), LASSO,
    % Ridge, Elastic Net and SuSiE

    x_train = convertvars(x_train, 1:width(x_train), @double);
    x_test = convertvars(x_test, 1:width(x_test), @double);
    
    if numel(opts.tool) > 1 % otherwise no need for model selection
        % outer loop CV
        cvp = cvpartition(size(x_train, 1), 'KFold', opts.outerkfold);
        
        [beta_in, stat_in] = deal(cell(numel(opts.outerkfold), 1));
        for i = 1:opts.outerkfold
            tt = tic;
            fprintf("outer fold:%d of %d\n", i, opts.outerkfold)
            x_train_in = x_train(cvp.training(i), :);
            x_test_in = x_train(cvp.test(i), :);
            [beta_in{i}, stat_in{i}, yhat_in] = trainInnerCV(x_train_in, x_test_in, opts);
        
            if i == 1 % initiate yhat cell
                yhat = nan(size(x_train, 1), size(yhat_in, 2));
            end
            yhat(cvp.test(i), :) = yhat_in;
            fprintf("fold time: %.1f\n\n", toc(tt))
        end
        
        % check CV MSE to find the best model
        % Note that x_train width is num features + 1, so AIC from this stat cannot be used
        stat.mse = evaluatePredictions(x_train.(opts.target), yhat, repmat(width(x_train), size(yhat, 2), 1), stat_in{1}.method);
        stat.mse = sortrows(stat.mse, "mse"); 
        fprintf("best model (lowest MSE in nested CV): %s\n\n", stat.mse.method(1))
        
        if opts.ignore_aic
            opts.tool = stat.mse.method(1).lower;
        else
            % check also based on lowest AIC: no need for nested CV.
            [~, stat.aic] = trainInnerCV(x_train, x_train, opts);
            stat.aic = sortrows(stat.aic, "AICc"); 
            fprintf("best model (lowest AICc): %s\n\n", stat.aic.method(1))
            opts.tool = union(stat.mse.method(1).lower, stat.aic.method(1).lower);
        end
    end
    
    % now train on the training set using the best method
    [beta_train, stat.best, ~, hp] = trainInnerCV(x_train, x_test, opts);
    stat.best.outcome = repmat(opts.target, numel(stat.best.method), 1); 
    stat.best = movevars(stat.best, ["outcome", "method"], "Before", 1);
    [stat.best.train_N, stat.best.test_N] = deal(nan(height(stat.best), 1));
    stat.best.train_N(1) = size(x_train, 1);
    stat.best.test_N(1) = size(x_test, 1);

    % whole dataset (training + test set)
    opts.hp = hp; % use hyperparameters tuned on whole training set
    beta_full = trainInnerCV(data, x_test, opts);
    [beta, shap, shap_idx] = deal(cell(numel(beta_train), 1));
    for i = 1:numel(beta_train)
        beta{i, 1} = table(["Intercept"; opts.featureNames], beta_train{i},...
            beta_full{i}, VariableNames=["feature", "train", "full"]);

        % SHAP values
        if opts.shap
            fprintf('\n\tcalculating SHAP values...')
            stime = tic;
            fnc = @(x) [ones(size(x,1),1),x]*beta_train{i};
            shap_idx{i, 1} = randsample(1:height(x_train), min(opts.shap_subsample, height(x_train)))';
            shap_train = x_train(shap_idx{i}, :); % only a subset for better performance
            shap_train = shap_train{:, 1:end-1};
            explainer = shapley(fnc, shap_train, ...
            UseParallel=opts.parallel, Method=opts.shap_method, ...
            CategoricalPredictors=opts.CategoricalPredictors, ...
            QueryPoints=shap_train);
            % shap{i, 1} = zeros(size(shap_train, 1), size(shap_train, 2));
            % for k = 1:size(shap_train, 1)
            %     explainer = fit(explainer, shap_train(k, :), UseParallel=opts.parallel);
            %     shap{i, 1}(k, :) = explainer.ShapleyValues.ShapleyValue;
            % end
            % beta{i, 1}.shap_mean_abs = [0, mean(abs(shap{i}), 1)]';
            % tcols = colnames(x_train);
            % shap{i} = array2table(shap{i}, VariableNames=tcols(1:end-1)+".shap");

            beta{i, 1}.shap_mean_abs = explainer.MeanAbsoluteShapley;
            shap{i} = array2table(explainer.ShapleyValues.ShapleyValue', VariableNames=explainer.ShapleyValues.Predictor+".shap");
            if isfield(opts, "raw_x_train")
                % use data in original scale
                shap{i} = [shap{i}, opts.raw_x_train(shap_idx{i}, :)];
            else
                shap{i} = [shap{i}, x_train(shap_idx{i}, :)];
            end
            fprintf("\b\b done(%.1f sec)\n", toc(stime))
        end

    end
    
else 
    
    %@13SEP2024: 'toolset' overrides 'tool' and instead runs a nestedCV
    %(recommended)
    if isfield(opts, "toolset")
        opts.tool = "auto_sep";
    end

    for k = 1:numel(opts.tool)
        fprintf('\t%d-%s', k, opts.tool(k))
        
        mltime = tic;
        [tmp.model{k, 1}, tmp.RMSEvalidation(k, 1),...
            tmp.RMSEtest(k, 1), tmp.yhat{k, 1}] = ...
            trainModel(x_train, x_test, data, opts.tool(k), opts);
        tt = toc(mltime);
        fprintf(' (%.1f sec)\n', tt)
        
    end

    yhat = horzcat(tmp.yhat{:}); % CV yhat on x_train

    % note that this may not be applicable to all approaches (e.g.
    % intercept included in numP or the way AIC is calculated)
    stat.mse = evaluatePredictions(x_train.(opts.target), yhat, ...
    repmat(width(x_train), size(yhat, 2), 1), opts.tool);

    [stat.mse, idx] = sortrows(stat.mse, "mse"); 
    if height(stat.mse) > 1
        fprintf("best model (lowest MSE in nested CV): %s\n\n", stat.mse.method(1))
    end
    
    beta = tmp.model{idx(1)}; % report only the best model
    % clear tmp
    
    % check the performance on the best model using x_test
    yhat = beta.predictFcn(x_test);
    stat.best = evaluatePredictions(x_test.(opts.target), yhat, ...
        repmat(width(x_train), size(yhat, 2), 1), stat.mse.method(idx(1)));
    stat.best.train_N(1) = size(x_train, 1);
    stat.best.test_N(1) = size(x_test, 1);

    % SHAP values
    if opts.shap
        fprintf('\tcalculating SHAP values...')
        stime = tic;
        shap_idx = randsample(1:height(x_train), min(opts.shap_subsample, height(x_train)))';
        shap_train = x_train(shap_idx, :); % only a subset for better performance
        shap_train = shap_train(:, 1:end-1);
        warning("off", "all")
        explainer = shapley(beta.model, shap_train, ...
            UseParallel=opts.parallel, Method=opts.shap_method, ...
            CategoricalPredictors=opts.CategoricalPredictors, ...
            QueryPoints=shap_train);
        % shap = zeros(size(shap_train, 1), size(shap_train, 2));
        % for k = 1:size(shap_train, 1)
        %     explainer = fit(explainer, shap_train(k, :), UseParallel=opts.parallel);
        %     shap(k, :) = explainer.ShapleyValues.ShapleyValue;
        % end
        warning("on", "all")
        % beta.shap_mean_abs = mean(abs(shap), 1)';
        beta.shap_mean_abs = explainer.MeanAbsoluteShapley;
        shap = array2table(explainer.ShapleyValues.ShapleyValue', VariableNames=explainer.ShapleyValues.Predictor+".shap");
        if isfield(opts, "raw_x_train")
            % use data in original scale
            shap = [shap, opts.raw_x_train(shap_idx, :)];
        else
            shap= [shap, x_train(shap_idx, :)];
        end
        fprintf("\b\b done(%.1f sec)\n", toc(stime))
    else
        shap = []; shap_idx = [];
    end

end

stat.shap = shap;
stat.shap_idx = shap_idx;

% do the plotting: note that tarined model on the traning dataset is used
% for visualization purposes.
if opts.plotshap
    if iscell(shap), shap = shap{1}; end
    plotSHAP(shap, explainer, opts)
end

if opts.plotice
    opts.NumObservationsToSample = min(opts.NumObservationsToSample, height(x_train));
    plotICE(beta.model, opts)
end


end % END

%% subfunctions ===========================================================
function plotICE(model, opts)
% full model by default is used.

if isa(model, "RegressionLinear")
    return
end

opts.figname = regexprep(opts.figname, "." + opts.format + "$", "");

if isfield(opts, "WindowPosition")
    fig = figure(Position=opts.WindowPosition, Units="pixels");
elseif isfield(opts, "WindowState")
    fig = figure("WindowState", opts.WindowState);
else
    fig = figure;
end

ti = tiledlayout(fig, "flow", "TileSpacing", opts.TileSpacing, "Padding", "compact");
cols = opts.featureNames;
for k = 1:numel(cols)
    ax = nexttile(ti);
    plotPartialDependence(model, cols(k), "Conditional", "centered", ...
        "UseParallel", opts.parallel, ...
        "NumObservationsToSample", opts.NumObservationsToSample, ...
        "Parent", ax)
    
    ax.Title.String = "";
    ax.FontName = opts.fontname;
    ax.FontSize = opts.fontsize;
    if opts.grid
        grid(ax, "on")
    end

    if opts.box
        ax.Box = "on";
    end
end

if opts.save
    exportgraphics(fig, opts.figname + ".ICE." + opts.format, "Resolution", opts.resolution)
end

if opts.savefig
    savefig(fig, opts.figname + ".ICE.fig")
end


end

%% ------------------------------------------------------------------------
function plotSHAP(shap, explainer, opts)

if opts.colormap == "shap" % deprecated
    numLevels = 256;
    ccmap = zeros(numLevels, 3);
    
    % Determine the color transition points
    tpoint = round(numLevels / 2);
    
    % Assign shades of blue for negative values
    ccmap(1:tpoint, 1) = linspace(0, 0.5, tpoint)';
    ccmap(1:tpoint, 3) = linspace(1, 0.8, tpoint)';
    
    % Assign shades of red for positive values
    ccmap(tpoint+1:end, 1) = linspace(0.5, 1, numLevels-tpoint)';
    ccmap(tpoint+1:end, 3) = linspace(0.8, 0, numLevels-tpoint)';
else
    ccmap = get_colormap(opts.colormap, 256);
end

%@12SEP2024: Summary Plot (needs MATLAB > R2024)
if isfield(opts, "WindowPosition")
    fig3 = figure(Position=opts.WindowPosition, Units="pixels");
elseif isfield(opts, "WindowState")
    fig3 = figure("WindowState", opts.WindowState);
else
    fig3 = figure;
end
ax3 = axes(fig3);
swarmchart(ax3, explainer);
axis(ax3, "square");
colormap(ccmap);
ax3.FontName = opts.fontname;
ax3.FontSize = opts.fontsize + 2;
ax3.YAxis.Label.String = "";
ax3.YAxis.FontSize = opts.fontsize + 4;
cbar = findobj(fig3.Children, "Type", "ColorBar");
cbar.Label.String = "";
cbar.FontName = opts.fontname;
ax3.FontSize = opts.fontsize + 2;
if opts.save
    exportgraphics(fig3, opts.figname + ".summarySHAP." + opts.format, "Resolution", opts.resolution)
end

if opts.savefig
    savefig(fig3, opts.figname + ".summarySHAP.fig")
end


opts.figname = regexprep(opts.figname, "." + opts.format + "$", "");

cols =  colnames(shap);
target = cols(end);
cols(end) = [];
cols = unique(regexprep(cols, ".shap$", ""));

if isfield(opts, "WindowPosition")
    fig = figure(Position=opts.WindowPosition, Units="pixels");
elseif isfield(opts, "WindowState")
    fig = figure("WindowState", opts.WindowState);
else
    fig = figure;
end

ti = tiledlayout(fig, "flow", "TileSpacing", opts.TileSpacing, "Padding", "compact");
for k = 1:numel(cols)
    ax = nexttile(ti);
    sval = shap.(cols(k) + ".shap");
    scatter(ax, shap.(cols(k)), shap.(target), opts.markersize, sval, "filled")
    xlabel(ax, cols(k));
    ylabel(ax, target);
    title(ax, sprintf("mean SHAP = %.3f",  mean(abs(shap.(cols(k) + ".shap")))))
    colormap(ax, ccmap)
    cb = colorbar(ax);
    cb.Label.String = "SHAP";
    ax.FontName = opts.fontname;
    ax.FontSize = opts.fontsize;
    if opts.grid
        grid(ax, "on")
    end

    if opts.box
        ax.Box = "on";
    end
end

if opts.save
    exportgraphics(fig, opts.figname + ".SHAP." + opts.format, "Resolution", opts.resolution)
end

if opts.savefig
    savefig(fig, opts.figname + ".SHAP.fig")
end

% create importance SHAP plot
mean_abs_shap = varfun(@(x)mean(abs(x)), shap(:, cols + ".shap"));
mean_abs_shap.Properties.VariableNames = regexprep(colnames(mean_abs_shap), "^fun_", "", "ignorecase");
mean_abs_shap.Properties.VariableNames = regexprep(colnames(mean_abs_shap), ".shap$", "");
mean_abs_shap = rows2vars(mean_abs_shap);
mean_abs_shap.Properties.VariableNames = ["feature", "shap"];
mean_abs_shap.feature = string(mean_abs_shap.feature);
mean_abs_shap = sortrows(mean_abs_shap, "shap", "ascend");
if isfield(opts, "WindowPosition")
    fig2 = figure(Position=opts.WindowPosition, Units="pixels");
elseif isfield(opts, "WindowState")
    fig2 = figure("WindowState", opts.WindowState);
else
    fig2 = figure;
end
ax2 = axes(fig2);
ax2h = barh(ax2, mean_abs_shap.shap);
ax2h.FaceColor = ccmap(1, :);
ax2.YTick = 1:height(mean_abs_shap);
ax2.YTickLabel = mean_abs_shap.feature;
ax2.XLabel.String = "mean(|SHAP value|)";
ax2.FontName = opts.fontname;
ax2.FontSize = opts.fontsize + 2;
if opts.grid
    grid(ax2, "on")
end

if opts.box
    ax2.Box = "on";
end
axis(ax2, "square")
text(mean_abs_shap.shap, 1:height(mean_abs_shap), ...
    compose("%.3f", mean_abs_shap.shap), ...
    Color=ccmap(1, :), FontSize=opts.fontsize);

if opts.save
    exportgraphics(fig2, opts.figname + ".importanceSHAP." + opts.format, "Resolution", opts.resolution)
end

if opts.savefig
    savefig(fig2, opts.figname + ".importanceSHAP.fig")
end

end

%% ------------------------------------------------------------------------
function [beta, stat, yhat, hp] = trainInnerCV(x_train, x_test, opts)

if ~isfield(opts, "tool")
    opts.tool = "all";
end

hp = struct;
lambda = opts.elnet.lambda;
alpha = opts.elnet.alpha; 

% convert all features to numeric (for ridge)
x_train = convertvars(x_train, 1:width(x_train), @double);
x_test = convertvars(x_test, 1:width(x_test), @double);

% convert to matrix form
y_train = x_train.(opts.target);
x_train.(opts.target) = [];
x_train = x_train{:, :};

beta = ({});

% Ridge
if any(ismember(opts.tool, ["all", "ridge"]))
    fprintf('\tridge...')
    tt = tic;
    if isfield(opts, 'hp') % use already tuned hyperparameters (recommended only for final model trained on whole dataset)
        tunedLambda = opts.hp.ridge.lambda;
    else
        cvp = cvpartition(numel(y_train), 'KFold', opts.kfold);
        yhat = nan(numel(y_train), numel(lambda));
        for fold = 1:opts.kfold
            br = ridge(y_train(cvp.training(fold)), x_train(cvp.training(fold), :), lambda, 0);
            yhat(cvp.test(fold), :) = br(1, :) + x_train(cvp.test(fold), :)*br(2:end, :);
        end
        MSE = sum((yhat - y_train).^2)./size(yhat, 1);
        [~, minMSEidx] = min(MSE);
        tunedLambda = lambda(minMSEidx);
    end

    beta{1} = ridge(y_train, x_train, tunedLambda, 0);
    hp.ridge.lambda = tunedLambda;
    fprintf('\b\b done (%.1f sec)\n', toc(tt)) 
end

% lasso
if any(ismember(opts.tool, ["all", "lasso"]))
    fprintf('\tLASSO...')
    tt = tic;
    if isfield(opts, 'hp') % use already tuned hyperparameters (recommended only for final model trained on whole dataset)
        tunedLambda = opts.hp.lasso.lambda;
    else
        tunedLambda = lambda;
    end

    [blasso, fitinfo] = lasso(x_train, y_train, 'CV', opts.kfold, 'Options', opts.statset, ...
        'PredictorNames', opts.featureNames, 'Lambda', tunedLambda, 'Standardize', true);
    beta{2} = [fitinfo.Intercept(fitinfo.IndexMinMSE); blasso(:, fitinfo.IndexMinMSE)];
    hp.lasso.lambda = fitinfo.LambdaMinMSE;
    fprintf('\b\b done (%.1f sec)\n', toc(tt)) 
end

% Elastic net
if any(ismember(opts.tool, ["all", "elnet"]))
    fprintf('\tElastic Net...')
    tt = tic;
    
    if isfield(opts, 'hp') % use already tuned hyperparameters (recommended only for final model trained on whole dataset)
        tunedLambda = opts.hp.elnet.lambda;
        alphaTuned = opts.hp.elnet.alpha;
    else
        minMSE = inf;
        alphaTuned = nan;
        for j = 1:numel(alpha)
            [~, fitinfoa] = lasso(x_train, y_train, 'CV', opts.kfold, 'Options', opts.statset,...
                'PredictorNames', opts.featureNames, 'Alpha', alpha(j), ...
                'Lambda', lambda, 'Standardize', true);
            if fitinfoa.MSE(fitinfoa.IndexMinMSE) < minMSE
                minMSE = fitinfoa.MSE(fitinfoa.IndexMinMSE);
                alphaTuned = alpha(j);
            end
        end
        tunedLambda = lambda;
    end
    
    [belnet, fitinfoa] = lasso(x_train, y_train, 'CV', opts.kfold, 'Options', opts.statset,...
            'PredictorNames', opts.featureNames, 'Alpha', alphaTuned, ...
            'Lambda', tunedLambda, 'Standardize', true);
    beta{3} = [fitinfoa.Intercept(fitinfoa.IndexMinMSE); belnet(:, fitinfoa.IndexMinMSE)];
    hp.elnet.alpha = alphaTuned;
    hp.elnet.lambda = fitinfoa.LambdaMinMSE;
    fprintf('\b\b done (%.1f sec)\n', toc(tt)) 
end

if any(ismember(opts.tool, ["all", "susie"]))
    fprintf('\tSuSiE and L0-learn...')
    tt = tic;
    save('susie.raw.mat', 'x_train', 'y_train');
    rcode = string;
    rcode(1, 1) = "y <- R.matlab::readMat('susie.raw.mat')$y.train";
    rcode(2, 1) = "X <- R.matlab::readMat('susie.raw.mat')$x.train";
    rcode(3, 1) = "library(susieR)";
    % res1 <- susie(X,y,L=100, max_iter = 300)
    rcode(4, 1) = "res <- susie_auto(X, y)"; % automatic version by a 3-step strategy
    rcode(5, 1) = "beta <- coef(res)";
    rcode(6, 1) = "pip <- susie_get_pip(res)";
    rcode(7, 1) = "write.table(cbind(beta, c(0,pip)), 'susie.beta.txt', quote = F, row.names = F, col.names = T)";
    % rcode(8, 1) = "library(L0Learn)";
    % rcode(9, 1) = "cvfit <- L0Learn.cvfit(X, y, nFolds=" + opts.kfold + ", seed=1, penalty='L0L2', nGamma=10, gammaMin=1e-6, gammaMax=.9)";
    % rcode(10, 1) = "mingammaIdx <- which.min(lapply(cvfit$cvMeans, min))";
    % rcode(11, 1) = "minLambdaIndex <- which.min(cvfit$cvMeans[[mingammaIdx]])";
    % rcode(12, 1) = "optimalLambda = cvfit$fit$lambda[[mingammaIdx]][minLambdaIndex]";
    % rcode(13, 1) = "beta<- coef(cvfit, lambda=optimalLambda, gamma=cvfit$fit$gamma[mingammaIdx])";
    % rcode(14, 1) = "write.table(as.numeric(beta), 'l0learn.beta.txt', quote = F, row.names = F, col.names = F)";
    writematrix(rcode, "susie_featureselection.r", 'QuoteStrings', false, 'FileType','text')
    MATLAB2Rconnector("susie_featureselection.r", 'delr', true);
    delete('susie.raw.mat');
    susie = readtable('susie.beta.txt');
    susie.Properties.VariableNames = {'beta', 'pip'};
    susie.name = ["Intercept"; opts.featureNames];
    % l0learn =  table(susie.name, readmatrix('l0learn.beta.txt'), ...
    %     'VariableNames', {'name', 'beta'});
    delete('susie.beta.txt')
    % delete('l0learn.beta.txt')
    beta{4} = susie.beta;
    % beta{5} = l0learn.beta;
    fprintf('\b\b done (%.1f sec)\n', toc(tt)) 
end

if any(ismember(opts.tool, ["all", "plsr"]))
    % Partial least-squares (PLS) regression
    % Note that PLS loadings/coefficients do not tell anything directly
    % about the importance of features, but to what extend each contributes
    % to specific PCs. See Google Notes. Hence, I stick with betas here
    % since they collectively show the contribution of top K (here 10) PCs.
    fprintf('\tPLSR...')
    tt = tic;
    [~,~,~,~, beta{5}] = plsregress(x_train, y_train, min(10, size(x_train, 2))); % with 10 PCs
    fprintf('\b\b done (%.1f sec)\n', toc(tt))  
end

idx = cellfun(@isempty, beta);
beta(idx) = [];

% check performance on test set 
yhat = cell(numel(beta), 1);
numP = zeros(numel(beta), 1);
for j = 1:numel(beta)
    yhat{j} = beta{j}(1) + x_test{:, 1:end-1}*beta{j}(2:end);
    numP(j, 1) = sum(beta{j} ~= 0);
end

methods = ["Ridge"; "LASSO"; "ElNet"; "SuSiE"; "PLSR"];
yhat = horzcat(yhat{:});
stat = evaluatePredictions(x_test.(opts.target), yhat, numP, methods(~idx));

end

%% ------------------------------------------------------------------------
function stat = evaluatePredictions(y, yhat, numP, methods)
% AIC (and other IC) can be used to compare penalized models accounting for
% both sample size and number of parameters.
% https://stats.stackexchange.com/questions/442121/aic-and-its-degrees-of-freedom-for-linear-regression-models
% https://stats.stackexchange.com/questions/479603/adjusting-the-number-of-parameters-for-aic-bic-calculation-in-case-of-correlat
% 
% However, note that AIC (on whole training set) and CV MSE are both
% asymptotically equivalent, so there is no point in having mean of AIC
% from CV, and MSE does that. So, it's better to use AIC on whole training
% set to compare models (keep in mind the differences between assumptions
% each make). So, AIC is asymptotically equivalent to leave-one-out CV
% error [see e.g. http://www.petrkeil.com/?p=836 ], so using AIC as a
% computationally efficient proxy for CV is reasonable.

% https://stats.stackexchange.com/questions/319666/aic-with-test-data-is-it-possible
% https://stats.stackexchange.com/questions/573694/how-does-lmstepaic-work-in-caret-when-using-cross-validation
% https://stats.stackexchange.com/questions/377527/stepwise-aic-does-there-exist-controversy-surrounding-this-topic

% Note: numc is number of features + 1 (intercept) for regression models

stat = struct;

if nargin < 4
    methods = "method." + (1:numel(numP))';
else
    if isrow(methods), methods = methods'; end
end


n = size(yhat, 1);
for i = 1:size(yhat, 2)
    sst = sum((y-mean(y)).^2);
    sse = sum((y-yhat(:, i)).^2);
    stat.mse(i, 1) = sse/n; % fitlm MSE is sse/dfe where dfe is error dfe: n - numc
    stat.r2adj(i, 1) = 1 - ((n-1)/(n-numP(i)))*(sse/sst);
    stat.r2(i, 1) = corr(y, yhat(:, i)).^2;
    stat.AIC(i, 1) = n*log(sse/n)+2*numP(i); % see here https://stats.stackexchange.com/questions/43733/what-is-the-difference-between-aic-and-extractaic-in-r
    stat.AICc(i, 1) = stat.AIC(i) + (2*numP(i)*(numP(i) + 1))/(n - numP(i) - 1); % for small samples (taken from MATLAB aicbic)
end

stat.n = repmat(n, size(yhat, 2), 1);
if isrow(numP), numP = numP'; end
stat.numP = numP;
stat.method = methods;
stat = struct2table(stat);

% stat.method = ["Ridge"; "LASSO"; "ElNet"; "SuSiE"; "L0Learn"; "PLSR"];
end % END

%% ------------------------------------------------------------------------
function [mdl, vRMSE, tRMSE, yhat] = trainModel(tab, testSet, ds, ml, inopts)
% a modified version of an automatic generated code from RegressionLearner
% app. 

KFolds = inopts.kfold;
response = inopts.target;

tab = movevars(tab, response, "After", width(tab)); % last column is response
testSet = movevars(testSet, response, "After", width(tab));
ds = movevars(ds, response, "After", width(ds));
features = tab.Properties.VariableNames(1:end-1);

% Bayesian optimizer options
bopts = struct('UseParallel', inopts.parallel, 'Verbose', 0, 'ShowPlots', false, ...
        'AcquisitionFunctionName','expected-improvement-plus');

%@13SEP2024: proper nested CV when subset of "auto" tagged algorithms are
%selected.

% Train the model for prediction
if startsWith(ml, "auto") % performs hyperparameter optimization for all ml types except gam

    % nested CV: outer loop splits tab into training/validation sets, and
    % each training set is sent to fitrauto for hyperparameter tuning using
    % ASHA optimizer using another CV as inner loop.
    learners = array2table( ...
        ["ensemble"    "CompactRegressionEnsemble" "fitrensemble"    
        "gp"          "CompactRegressionGP" "fitrgp"          
        "kernel"      "RegressionKernel" "fitrkernel"             
        "linear"      "RegressionLinear" "fitrlinear"              
        "net"         "CompactRegressionNeuralNetwork" "fitrnet"
        "svm"         "CompactRegressionSVM" "fitrsvm"
        "tree"        "CompactRegressionTree" "fitrtree"], ...
        "VariableNames", ["name", "class", "func"]);
    

    %@13SEP2024: subset of auto models
    if isfield(inopts, "toolset")
        learners(~ismember(learners.name, inopts.toolset), :) = [];
    end
    

    % learners = learners(4, :); % DEBUG
    % ml = "auto_sep";
    cvp = cvpartition(size(tab, 1), 'KFold', KFolds);

    if ml == "auto"
        bopts.Optimizer = "asha";
        bopts = rmfield(bopts, "AcquisitionFunctionName");
    elseif ml == "auto_sep" % don't use futrauto
        bopts.AcquisitionFunctionName = "expected-improvement-plus";
        bopts.Optimizer = "bayesopt";
        learners.RMSE = nan(height(learners), KFolds);
    end
    
    inRMSEv = table(nan(KFolds, 1), strings(KFolds, 1), ...
        VariableNames=["RMSE", "class"]);
    pbar = progressGen(KFolds);
    for i = 1:KFolds
        progressGen(pbar, i);
        tabtrain = tab(cvp.training(i), :);
        tabtest = tab(cvp.test(i), :);
        
        bopts.CVPartition = cvpartition(size(tabtrain.(response), 1), 'KFold', KFolds);

        if ml == "auto"
            tmpMdl = fitrauto(tabtrain, response, ...
               "OptimizeHyperparameters","all", "Learners", "all", ...
               "HyperparameterOptimizationOptions", bopts, ...
               "CategoricalPredictors", inopts.CategoricalPredictors);
            inRMSEv.RMSE(i) = sqrt(tmpMdl.loss(tabtest, response, "LossFun", 'mse'));
            inRMSEv.class(i) = string(class(tmpMdl));
        elseif ml == "auto_sep"
            % alternatively, loop over learners spearately
            for j = 1:height(learners)
                tmpMdl = feval(learners.func(j), tabtrain, response, ...
                    "OptimizeHyperparameters","all", ...
                    "HyperparameterOptimizationOptions", bopts, ...
                    "CategoricalPredictors", inopts.CategoricalPredictors);
                learners.RMSE(j, i) = sqrt(tmpMdl.loss(tabtest, response, "LossFun", 'mse'));
            end
        end

    end

    if ml == "auto"
        idx = contains(inRMSEv.class, ".");
        inRMSEv.class(idx) = extract(inRMSEv.class(idx), "." + lettersPattern + textBoundary("end"));
        inRMSEv.class = erase(inRMSEv.class, ".");
        inRMSEv = sortrows(inRMSEv, "RMSE", "ascend");
        bestLearner = learners(learners.class == inRMSEv.class(1), :);
    elseif ml == "auto_sep"
        learners.RMSE = mean(learners.RMSE, 2);
        learners = sortrows(learners, "RMSE", "ascend");
        bestLearner = learners(1, :);
    end
    
    % train the whole traning dataset on the best model from nested CV   
    fprintf(", best learner: %s", bestLearner.name)
    bopts.CVPartition = cvpartition(size(tab.(response), 1), 'KFold', KFolds);
    bopts.AcquisitionFunctionName = "expected-improvement-plus";
    bopts.Optimizer = "bayesopt";
    tmpMdl = feval(bestLearner.func, tab, response, ...
            "OptimizeHyperparameters","all", ...
            "HyperparameterOptimizationOptions", bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);
    
elseif startsWith(ml, "lm") % ---------------------------------------------
    type = "linear"; rmode = "off";        
    if ml == "lm_interaction"
        type = "interactions";
    elseif ml == "lm_robust"
        rmode = "on";
    end

    tmpMdl = fitlm(tab, type, 'RobustOpts', rmode, ...
        "CategoricalVars", inopts.CategoricalPredictors);

elseif startsWith(ml, "gam") % --------------------------------------------
    hp = "none";
    if ml == "gam_auto"
        hp = "auto";
    end
    tmpMdl = fitrgam(tab, response, 'OptimizeHyperparameters', hp, ...
        'HyperparameterOptimizationOptions', bopts, ...
        "CategoricalPredictors", inopts.CategoricalPredictors);

elseif startsWith(ml, "gpr_") % Gaussian process regression (GPR) ---------
    if ml == "gpr_exponential"
        krnl = "exponential";
    elseif ml == "gpr_rationalquadratic"
        krnl = "rationalquadratic";
    elseif ml == "gpr_squaredexponential"
        krnl = "squaredexponential";
    elseif ml == "gpr_matern52"
        krnl = "matern52";
    end

    tmpMdl = fitrgp(tab, response, 'BasisFunction', 'constant', ...
        'KernelFunction', krnl, 'Standardize', true, ...
        'OptimizeHyperparameters', "none", ...
        'HyperparameterOptimizationOptions', bopts, ...
        "CategoricalPredictors", inopts.CategoricalPredictors); 

elseif startsWith(ml, "tree_") % binary decision tree ---------------------
    if ml == "tree_fine"
        minLeafSize = 4;
    elseif ml == "tree_medium"
        minLeafSize = 12;
    elseif ml == "tree_coarse"
        minLeafSize = 36;
    end

    tmpMdl = fitrtree(tab, response, 'MinLeafSize', minLeafSize, ...
        'Surrogate', 'off', "CategoricalPredictors", inopts.CategoricalPredictors);

elseif startsWith(ml, "ensemble_") % ensemble of learners -----------------
    template = templateTree('MinLeafSize', 8);
    if ml == "ensemble_boosted"
        tmpMdl = fitrensemble(tab, response, 'Method', 'LSBoost', ...
            'NumLearningCycles', 30, 'Learners', template, 'LearnRate', 0.1, ...
            'OptimizeHyperparameters', "none", ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);

    elseif ml == "ensemble_bagged"
        tmpMdl = fitrensemble(tab, response, 'Method', 'Bag', ...
            'NumLearningCycles', 30, 'Learners', template, ...
            'OptimizeHyperparameters', "none", ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);

    elseif ml == "ensemble_auto"
        bopts.CVPartition = cvpartition(size(tab.(response), 1), 'KFold', KFolds);
        tmpMdl = fitrensemble(tab, response, 'OptimizeHyperparameters', "all", ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);
    end

elseif startsWith(ml, "svm_") % support vector machine --------------------
    if ml == "svm_auto"
        bopts.CVPartition = cvpartition(size(tab.(response), 1), 'KFold', KFolds);
        tmpMdl = fitrsvm(tab, response, ...
            'OptimizeHyperparameters', "all", ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);

    else
        responseScale = iqr(tab.(response));
        if ~isfinite(responseScale) || responseScale == 0.0
            responseScale = 1.0;
        end
        boxConstraint = responseScale/1.349;
        epsilon = responseScale/13.49;
        
        kscale = "auto";
        if ml == "svm_linear"
            krnl = "linear"; porder = [];
        elseif ml == "svm_quadratic"
            krnl = "polynomial"; porder = 2;
        elseif ml == "svm_cubic"
            krnl = "polynomial"; porder = 3;
        elseif ml == "svm_finegaussian"
            krnl = "gaussian"; porder = []; kscale = 0.5;
        elseif ml == "svm_mediumgaussian"
            krnl = "gaussian"; porder = []; kscale = 2;
        elseif ml == "svm_coarsegaussian"
            krnl = "gaussian"; porder = []; kscale = 8;
        end
    
        tmpMdl = fitrsvm(tab, response, 'KernelFunction', krnl, ...
            'PolynomialOrder', porder, 'KernelScale', kscale, ...
            'BoxConstraint', boxConstraint, 'Epsilon', epsilon, ...
            'Standardize', true, ...
            'OptimizeHyperparameters', "none", ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);
    end

elseif startsWith(ml, "nnet_") % neural network regression ----------------
    if ml == "nnet_narrow"
        lsize = 10;
    elseif ml == "nnet_medium"
        lsize = 25;
    elseif ml == "nnet_wide"
        lsize = 100;
    elseif ml == "nnet_bi"
        lsize = [10, 10];
    elseif ml == "nnet_tri"
        lsize = [10 10 10];
    end
    
    if ml == "nnet_auto"
        params = hyperparameters('fitrnet', tab(:, features), tab.(response));
        params(1).Range = [1 5];
        params(10).Optimize = true;
        params(11).Optimize = true;
        for ii = 7:11
            params(ii).Range = [1 400];
        end
        bopts.CVPartition = cvpartition(size(tab.(response), 1), 'KFold', KFolds);
        tmpMdl = fitrnet(tab, response, 'OptimizeHyperparameters', params, ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);
    else
        tmpMdl = fitrnet(tab, response, 'LayerSizes', lsize, ...
            'Activations', 'relu', 'Lambda', 0, 'IterationLimit', 1000, ...
            'Standardize', true, ...
            'OptimizeHyperparameters', "none", ...
            'HyperparameterOptimizationOptions', bopts, ...
            "CategoricalPredictors", inopts.CategoricalPredictors);
    end
end

% train on full model
if startsWith(ml, "lm")
    fmdl = fitlm(ds, type, 'RobustOpts', rmode, "CategoricalVars", inopts.CategoricalPredictors);
else
    % see https://se.mathworks.com/matlabcentral/answers/481300-how-to-re-train-a-model-optimized-by-bayesian-optimization-on-new-data
    if isa(tmpMdl, "RegressionLinear")
        hps = tmpMdl.ModelParameters;
        rmfis = ["BatchIndex", "LossFunction", "InitialBeta",...
            "InitialBias", "Method", "Type", "VerbosityLevel"];
        hps.Beta = hps.InitialBeta;
        hps.Bias = hps.InitialBias;
        fis = string(fieldnames(hps));
        for k = 1:numel(fis)
            if isempty(hps.(fis(k)))
                rmfis = union(rmfis, fis(k));
            end
        end
        hps = rmfield(hps, rmfis);
        hps.CategoricalPredictors = inopts.CategoricalPredictors;
        hps2 = namedargs2cell(hps);
        fmdl = fitrlinear(ds, response, hps2{:});
    else
        tmp = classreg.learning.FitTemplate.makeFromModelParams(tmpMdl.ModelParameters);
        fmdl = fit(tmp, ds, char(response));
    end

end

% Create the result struct with predict function
pExtract = @(t) t(:, features);
pFnc = @(x) predict(tmpMdl, x);
mdl.predictFcn = @(x) pFnc(pExtract(x));
mdl.RequiredVariables = features;
mdl.model = tmpMdl;
mdl.full = fmdl;

% Perform cross-validation ------------------------------------------------
if ~startsWith(ml, "lm") && ...
    ~any(strcmp(class(mdl.model), ["RegressionLinear", "RegressionKernel"]))

    pModel = crossval(mdl.model, 'KFold', KFolds);
    try
        vRMSE = sqrt(kfoldLoss(pModel, 'LossFun', 'mse'));
    catch
        vRMSE = sqrt(kfoldLoss(pModel));
    end
    yhat = kfoldPredict(pModel);

    % for test set (to be reported for the best model based on validation CV)
    tRMSE = sqrt(mdl.model.loss(testSet, response, "LossFun", 'mse'));
    % mdl.model = mdl.model.compact;
    % mdl.full = mdl.full.compact;

    return
end

% for fitlm
cvp = cvpartition(size(tab.(response), 1), 'KFold', KFolds);

if strcmp(class(mdl.model), "RegressionLinear") % keep optimized parameters
    % by default in fitrauto, Learner, Regularization and Lambda are
    % optimized
    lopts.Learner = mdl.model.Learner;
    lopts.Lambda = mdl.model.Lambda;
    lopts.Regularization = mdl.model.ModelParameters.Regularization;
    lopts.Solver = mdl.model.ModelParameters.Solver;
elseif strcmp(class(mdl.model), "RegressionKernel")
    lopts.Learner = mdl.model.Learner;
    lopts.Lambda = mdl.model.Lambda;
    lopts.NumExpansionDimensions = mdl.model.NumExpansionDimensions;
    lopts.KernelScale = mdl.model.KernelScale;
    lopts.Epsilon = mdl.model.ModelParameters.Epsilon;
end

% Initialize the predictions to the proper sizes
yhat = tab.(response);
for fold = 1:KFolds

    % train a new model for each fold
    if startsWith(ml, "lm")
        tmpMdl = fitlm(tab(cvp.training(fold), :), type, ...
            'RobustOpts', rmode, "CategoricalVars", inopts.CategoricalPredictors);
%     elseif startsWith(ml, "svm_")
%         tmpMdl = fitrsvm(tab(cvp.training(fold), :), response, ...
%             'KernelFunction', krnl, 'PolynomialOrder', porder, ...
%             'KernelScale', kscale, 'BoxConstraint', boxConstraint, ...
%             'Epsilon', epsilon, 'Standardize', true, ...
%             'OptimizeHyperparameters', "none", ...
%             'HyperparameterOptimizationOptions', bopts);
    elseif strcmp(class(mdl.model), "RegressionLinear") % "auto"
         tmpMdl = fitrlinear(tab(cvp.training(fold), :), response, ...
             "Learner", lopts.Learner, ...
             "Regularization", lopts.Regularization, ...
             "Lambda", lopts.Lambda, ...
             "Solver", lopts.Solver, ...
             'OptimizeHyperparameters', "none", ...
             'HyperparameterOptimizationOptions', bopts, ...
             "CategoricalPredictors", inopts.CategoricalPredictors);
    elseif strcmp(class(mdl.model), "RegressionKernel") % "auto"
        if strcmp(lopts.Epsilon, 'auto')
            tmpMdl = fitrkernel(tab(cvp.training(fold), :), response, ...
             "Learner", lopts.Learner, ...
             "NumExpansionDimensions", lopts.NumExpansionDimensions, ...
             "Lambda", lopts.Lambda, ...
             "KernelScale", lopts.KernelScale, ...
             'OptimizeHyperparameters', "none", ...
             'HyperparameterOptimizationOptions', bopts, ...
             "CategoricalPredictors", inopts.CategoricalPredictors);
        else
            tmpMdl = fitrkernel(tab(cvp.training(fold), :), response, ...
             "Learner", lopts.Learner, ...
             "NumExpansionDimensions", lopts.NumExpansionDimensions, ...
             "Lambda", lopts.Lambda, ...
             "KernelScale", lopts.KernelScale, ...
             "Epsilon", lopts.Epsilon, ...
             'OptimizeHyperparameters', "none", ...
             'HyperparameterOptimizationOptions', bopts, ...
             "CategoricalPredictors", inopts.CategoricalPredictors);
        end
    end

    % Compute validation predictions
    foldPredictions = tmpMdl.predict(tab(cvp.test(fold), 1:end-1));

    % Store predictions in the original order
    yhat(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
vRMSE = sqrt(sum((yhat - tab.(response)).^2, 'omitnan')/numel(tab.(response)));

% Compute test RMSE
testPreds = mdl.predictFcn(testSet(:, features));
tRMSE = sqrt(sum((testPreds - testSet.(response)).^2, 'omitnan')/numel(testSet.(response)));

mdl.model = mdl.model.compact;
end 