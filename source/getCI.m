function tab = getCI(tab, ciflag, check_or)
% Is a subfunction of effPlotter and rarePlotter functions for
% calculating/extracting confidence intervals of estimates.
% 
% check_or is logical vector indicating binary outcomes. This is an
% extension of binaryflag when summary stats in each table belong to
% different traits (mixture of binary and quantiative traits).

if ciflag
    ci = tab.ci;
    ci(ismissing(ci)) = deal("[0, 0]");
    tab.ci = [];
    missedData = ismissing(ci);
    ci = cellfun(@(x)sscanf(x, '[%f, %f]').', ci, 'uni', false);
    ci = vertcat(ci{:});
    eff = tab.b;
    eff(missedData) = [];
    
    f_roundoff = ci(:, 2) <= eff | ci(:, 1) >= eff | any(abs(ci - 1) < 1e-6, 2) | any(abs(ci - 0) < 1e-6, 2);
    if any(f_roundoff)
        se = tab.se;
        se(missedData) = [];
        for i = 1:numel(check_or)
            if f_roundoff(i)
                if check_or(i)
                    ci(i, :) = [exp(log(eff(i)) - se(i).*1.96), ...
                        exp(log(eff(i)) + se(i).*1.96)];
                else
                    ci(i, :) = [eff(i) - se(i).*1.96, ...
                        eff(i) + se(i).*1.96];
                end
            end
        end
    end
else
    ci = nan(numel(check_or), 2);
    se = tab.se;
    se(se == 0) = eps;
    eff = tab.b;
    eff(check_or) = log(eff(check_or));
    
    %@16JULY2024: estimate df and t_crit for CI
    t_crit = 1.96.*ones(height(tab), 1);
    if any(colnames(tab) == "p")
        t_stat = eff./se;
        
        for k = 1:numel(t_stat)
            if isnan(t_stat(k)), continue; end
            if ~any(colnames(tab) == "df"), continue; end
            if ismissing(tab.df(k)), continue; end
            t_crit(k) = tinv(1 - 0.05/2, tab.df(k));

            % fun = @(df) 2 * (1 - tcdf(abs(t_stat(k)), df)) - tab.p(k);
            % 
            % % if returns error, use 1.96 (chisq)
            % df = fzero(fun, 10, optimset("Display", "off"));
            % if ~isnan(df)
            %     t_crit(k) = tinv(1 - 0.05/2, df);
            % end
        end
    end

    for i = 1:numel(check_or)
        if check_or(i)
            ci(i, :) = [exp(eff(i) - se(i).*t_crit(i)),...
                exp(eff(i) + se(i).*t_crit(i))];
        else
            ci(i, :) = [eff(i) - se(i).*t_crit(i),...
                eff(i) + se(i).*t_crit(i)];
        end
    end
end

tab.ci_l = ci(:, 1); tab.ci_h = ci(:, 2);

end