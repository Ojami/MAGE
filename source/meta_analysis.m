function [fixed_effect, random_effect, hetero, weights, sanityP] = ...
    meta_analysis(ES, SE_P, sm, varargin)
% Perform meta-analysis using "metagen" function of "meta" package in R.
% For input details: https://www.rdocumentation.org/packages/meta/versions/4.9-6/topics/metagen
% For info: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/random.html
% INPUTS:
%   - ES: Numeric vector of estimates. In case of OR, RR, etc. the
%         estimates must be in form of log OR (i.e. they must be
%         regression coefficients).
%   - SE_P: Numeric vector of standard errors; if not available, can be
%         estimated from X2 statistics of regression p-value:
%         X2 = (ES/SE)^2 --> Z = ES/SE --> SE = ES/Z
%         
%   - sm: A character string indicating underlying summary measure, e.g.,
%         "RD", "RR", "OR", "ASD", "HR", "MD", "SMD", or "ROM"
%         - "RD":  Risk difference
%         - "RR":  Risk ratio
%         - "OR":  Odds ratio (use for binary logistic regression)
%         - "ASD": Arcsine difference
%         - "HR":  Hazard Ratio
%         - "MD":  Mean Difference (linear regression with the same scale)
%         - "SMD": Standardized mean difference (linear regression for different scales)
%         - "ROM:  Ratio of Means
% 
% OPTIONAL INPUTS:  
%   - pflag:     A logical indicating if P-values are provided instead of 
%                SE (false)
%   - print:     A logical indicating whether meta-analyzed output should
%                be written to an Excel file (false).
%   - forest:    A logical to indicate if a forest plot is needed (false).
%   - forestopt: A struct containing drawing options for forest <UNDER
%                CONSTRUCTION>
% OUTPUTS:
%   -fixed_effect: a struct with fields: p, CI and estimate for fixed
%                  effect model, from metagen function
%   -random_effect: a struct with fields: p, CI and estimate for random
%                  effect model, from metagen function
%   -hetero: Heterogenity information.
%   -weights: Meta-analysis weights for both fixed- and random-effects
%             models.
%   -sanityP: p-values for each single study (for further confirmation of
%             meta-analysis)
% 
% Oveis Jamialahmadi, GU, Dec 2019.
%                         Nov 2020: Forset plot option was added.
% 
% NOTE: The only true statistical test for homogeneity is the Chi-Square
% test of Q. Tau-squared, the between-study variance, is derived from this
% Q and subsequently used to conduct random effects analysis. Under an
% inverse-variance model, there is no need to switch to a fixed effect
% model in the absence of heterogeneity, as the two results will be exactly
% identical. That is, the weights for studies will be based on the inverse
% of a study's variance, which only under a random effects model will
% include tau-squared. Hence, if that term is zero, the random effects and
% fixed effect model are the same. I-Squared above 50% can typically be
% interpreted as more than half of the total heterogeneity stems from
% between-study variance that cannot be explained by sampling error
% (alone).

% clc

p = inputParser;
p.CaseSensitive = false;
p.StructExpand = true;
validateIN = @(x) ( isnumeric(x) && ~isempty(x) && any(size(x) == 1) );
addRequired(p, 'ES', validateIN);
addRequired(p, 'SE_P', validateIN);
addRequired(p, 'sm', @(x) ( ~isempty(x) && isscalar(x) && (ischar(x) ||...
    isstring(x)) ));

addParameter(p, 'pflag', false, @islogical);
addParameter(p, 'print', false, @islogical);
addParameter(p, 'forest', false, @islogical);
addParameter(p, 'forestopt', [], @isstruct);

p.parse(ES, SE_P, sm, varargin{:});
p = p.Results;


if isrow(p.SE_P)
    p.SE_P = p.SE_P';
end
if isrow(p.ES)
    p.ES = p.ES';
end

% check for NaN values
fnan = isnan(p.ES);
if any(fnan)
    p.ES(fnan) = []; p.SE_P(fnan) = [];
    fprintf('%d study(ies) was(were) removed because of nan values.\n',...
        sum(fnan))
end

if p.pflag % P-values are provided in SE argument
    fprintf('SE was inferred from p-value.\n')
    Z = sqrt(pval2chi(p.SE_P, 1));
    p.SE_P = abs(p.ES./Z);
end

writematrix(p.ES, 'ES_MAT.txt') 
writematrix(p.SE_P, 'SE_MAT.txt') 

META_R(p.sm, pwd)
MATLAB2Rconnector('meta_analysis.r');

if exist('R2MAT_meta_analysis.txt', 'file')
    sanityP = double(split(readmatrix('original.pvlas.txt',...
        'OutputType', 'string'))); % for verification of estimated p-values 
    
    meta_raw = readtable('R2MAT_meta_analysis.txt');
    meta_raw = rows2vars(meta_raw);
    meta_raw.Properties.VariableNames = {'tmp', 'Estimate', 'Z-value', ...
            'CI Lower', 'CI Upper', 'P', 'Prediction interval'};
    meta_raw(:, 1) = [];
    meta_raw = varfun(@(x)double(string(x)), meta_raw);
    meta_raw.Properties.VariableNames = ...
        replace(meta_raw.Properties.VariableNames, 'Fun_', '');
    meta_raw.Properties.RowNames = {'Fixed-effect', 'Random-effects'};
    
    meta_raw.("Prediction interval") = string(meta_raw.("Prediction interval"));
    meta_raw.("Prediction interval") = ...
        repmat(compose("[%.2f, %.2f]", meta_raw.("Prediction interval")'), 2, 1);
    
    if p.print
        writetable(meta_raw, "meta_analysis.out"+rand(1)+".xlsx", ...
            'WriteRowNames', true)
    end

    fixed_effect.es = meta_raw.Estimate(1);
    fixed_effect.z = meta_raw.("Z-value")(1);
    fixed_effect.ci = [meta_raw.("CI Lower")(1), meta_raw.("CI Upper")(1)];
    fixed_effect.p = meta_raw.P(1);
    random_effect.es = meta_raw.Estimate(2);
    random_effect.z = meta_raw.("Z-value")(2);
    random_effect.ci = [meta_raw.("CI Lower")(2), meta_raw.("CI Upper")(2)];
    random_effect.p = meta_raw.P(2);
    
    % get weigths
    weights = readtable('R2MAT.weights.txt');
    weights = varfun(@(x)x/sum(x), weights);
    weights.Properties.VariableNames = ...
        regexprep(weights.Properties.VariableNames, '^Fun_', '');
    
    % get heterogeneity info
    hetero = readtable('R2MAT.hetero.txt', 'PreserveVariableNames', true);
    
    % get studies' CIs
    CI = readmatrix('R2MAT.ci.txt', 'Delimiter', '\t');
    if strcmpi(p.sm, "or")
        CI = exp(CI);
    end
    p.ci = CI;
    
    warning off
    delete('.RData')
    delete('meta_analysis.r.Rout')
    delete('original.pvlas.txt');
    delete('R2MAT_meta_analysis.txt')
    delete('R2MAT.weights.txt')
    delete('R2MAT.hetero.txt')
    delete('R2MAT.ci.txt')
    warning on
    
    if p.forest
        forest(p, meta_raw, weights, hetero)
    end

    delete('ES_MAT.txt')
    delete('SE_MAT.txt')
    
else
    error('R could not perform meta analysis. Check meta_analysis.r!')
end

end

%%
function META_R(sm, currDir)
currDir = strrep(currDir, '\', '/');
% read the data
fid = fopen('meta_analysis.r', 'w');
%create in a cell
Rcode = ({});
% write the content : here, the 12th line, new string.
Rcode{1} = ['setwd("',currDir,'")'];
Rcode{2} = 'library(meta)';
Rcode{5} = 'ES = as.matrix(read.table("ES_MAT.txt",sep = ";"))';
Rcode{6} = 'SE = as.matrix(read.table("SE_MAT.txt",sep = ";"))';
% Rcode{8} = ['out_data <- metagen(ES, SE, sm = "',upper(char(sm)),'",prediction=T)'];
Rcode{8} = ['out_data = tryCatch({metagen(ES, SE, sm = "', upper(char(sm)), '",prediction=T)}, error=function(err){metagen(ES, SE, sm = "', upper(char(sm)), '",prediction=T, method.tau = "PM")})'];
% MyText{8} = ['out_data <- metagen(ES, SE, sm = "',upper(char(sm)),'",method.tau = "PM")'];

if strcmpi(sm, "or")
    Rcode{10} = ['fixed_data <- c(exp(out_data$TE.fixed),',...
        'out_data$zval.fixed,exp(out_data$lower.fixed),',...
        'exp(out_data$upper.fixed),out_data$pval.fixed,'...
        'exp(out_data$lower.predict))'];
    Rcode{12} = ['random_data <- c(exp(out_data$TE.random),',...
        'out_data$zval.random,exp(out_data$lower.random),',...
        'exp(out_data$upper.random),out_data$pval.random,'...
        'exp(out_data$upper.predict))'];
else
    Rcode{10} = ['fixed_data <- c(out_data$TE.fixed,out_data$zval.fixed,',...
        'out_data$lower.fixed,out_data$upper.fixed,out_data$',...
        'pval.fixed, out_data$lower.predict)'];
    Rcode{12} = ['random_data <- c(out_data$TE.random,out_data$zval.',...
        'random,out_data$lower.random,out_data$upper.random,out_data',...
        '$pval.random, out_data$upper.predict)'];
end

% save fixed- / random- wieghts 
Rcode{13} = 'w = cbind(out_data$w.fixed, out_data$w.random)';
Rcode{14} = 'colnames(w) <- c("fixed","random")';
Rcode{15} = 'write.table(w,"R2MAT.weights.txt",sep ="\t",row.names = F,quote = F)';

% save I2, tau2 and Q stat (heterogenity info)
Rcode{16} = ['H = cbind(out_data$Q, out_data$pval.Q, out_data$tau2, ',...
    'out_data$se.tau2, out_data$lower.tau2, out_data$upper.tau2, ',...
    'out_data$I2, out_data$lower.I2, out_data$upper.I2)'];
Rcode{17} = ['colnames(H) = c("Q", "pval.Q", "tau2", "se.tau2", ',...
    '"lower.tau2", "upper.tau2", "I2", "lower.I2", "upper.I2")'];
Rcode{18} = 'write.table(H, "R2MAT.hetero.txt",sep ="\t", row.names = F, quote= F)';

Rcode{19} = ['write.table(cbind(fixed_data, random_data),',...
    '"R2MAT_meta_analysis.txt",row.names = F, quote= F)'];
Rcode{20} = 'write(out_data$pval, file = "original.pvlas.txt")';
Rcode{21} = ['write.table(cbind((out_data$lower), (out_data$upper)),',...
    ' "R2MAT.ci.txt", row.names = F, col.names = F, sep = "\t")'];


fprintf(fid,'%s\n',Rcode{:});
fclose(fid);
end