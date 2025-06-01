function qq(pvalues, opts)
% qq is equivalent to qq function of qqman in R, and plots observed vs
% expected p-values
% Info: https://rdrr.io/cran/qqman/src/R/qq.R
% Oveis Jamialahmadi, GU, Oct 2019
% 
% @10/08/2021 support for showing signal lables was added.
% @11JULY2024 now uses textRepel to find optimal positions for each label.

arguments
    pvalues double {mustBeVector}
    opts.labels {mustBeText, mustBeVector} = ""
    opts.showinsig (1,1) logical = true % to show few top insignificant labels as well (shown in black)
    opts.markersize (1,1) double = 20
    opts.save (1,1) logical = true
    opts.savefig (1,1) logical = false
    opts.savename {mustBeTextScalar} = "qq_plot"
    opts.resolution (1,1) double = 300
    opts.factor (1,1) {mustBeGreaterThanOrEqual(opts.factor, 1)} = 1.15;
    opts.n (1,1) double {mustBeGreaterThan(opts.n, 0)} = 5; % max insignifcant top hits to label
    opts.ci (1,1) logical = true
    opts.plot {mustBeMember(opts.plot, ["scatter", "plot"]), mustBeTextScalar} = "plot"
    opts.log10 (1,1) logical = false

    % debugging

    % DEPCRECATED
    opts.factorScale = 0.05; % linearly increase the arrow length by this size
    opts.maxAng = 15; % max iterations for text objects adjustment
    opts.margin = 1e-3; % allowable overlap margin between text obj's Extent property
    opts.resetAngle = true; % reset angle counter if any if collision conditions (including maxAng) are met (i.e. modify angle again). If false, only increases array length if iterations passed maxAng. 
end

if all(opts.labels ~= "") 
    if numel(opts.labels) ~= numel(pvalues)
        error('qq: labels and pvalues must have the same size!')
    else
        opts.annotate = true;
    end
else
    opts.annotate = false;
end

if ~opts.log10
    pvalues(pvalues < eps*realmin) = eps*realmin;
end

nanidx = ismissing(pvalues); % remove nan from pvalues
pvalues(nanidx) = [];

if opts.annotate 
    opts.labels(nanidx) = [];

    if opts.log10
        [pvalues, srtidx] = sort(pvalues, "descend");
    else
        [pvalues, srtidx] = sort(pvalues);
    end
    opts.labels = opts.labels(srtidx);
else
    pvalues = sort(pvalues);
end

if opts.log10
    obsPvals = pvalues;
    pvalues = 10.^-pvalues;
else
    obsPvals = -log10(pvalues);
end
expPvals = -log10(ppoints(numel(obsPvals))).';
maxpoint = max(expPvals);

close gcf force
fig = figure("WindowState", "maximized");
tiledlayout(fig, 1, 1, 'Padding', 'tight');
h = nexttile; 
xstart = min(union(obsPvals, expPvals));
l = line([xstart maxpoint],[xstart maxpoint], 'color', 'r', 'LineWidth', 1.1);
hold on
if opts.plot == "plot"
    plot(expPvals, obsPvals, 'k.', 'MarkerSize', opts.markersize)
else
    scatter(expPvals, obsPvals, opts.markersize, 'k', 'filled',...
        'Marker', 'o', ...
        'MarkerFaceAlpha', 0.6)
end
hold off
set(h,'box','on')
set(h, 'FontSize', 11, 'FontWeight', 'bold', 'LineWidth', 1.1, 'TickDir',...
    'out', 'TickLength', [0.005 0])

xlabel('Expected -log_{10} (\itp)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Observed -log_{10} (\itp)', 'FontSize', 12, 'FontWeight', 'bold')
[l.XData(2), l.YData(2)] = deal(max([h.XLim]));
% axis('tight')

if opts.ci
    n = 1:numel(expPvals);
    ci_l = betainv(0.025, n, flip(n));
    ci_h = betainv(0.975, n, flip(n));
    patch([expPvals; flip(expPvals)],...
        -log10([ci_l.'; flip(ci_h.')]), 'r', 'FaceAlpha', 0.25, ...
        'EdgeAlpha', 0, 'Parent', h);
end

% Get inflation factor
L = getInflation(pvalues);
title(h, ['\lambda= ',num2str(L,'%.2f')], 'FontSize', 13, 'FontWeight', 'bold')

if opts.annotate 
    if numel(expPvals) < 3e4 % use Bonferroni
        opts.thresh = -log10(0.05/numel(expPvals));
    else
        opts.thresh = -log10(5e-8); % gwas significance level
    end

    % check BH FDR
    pFDR = mafdr(pvalues, 'BH', true);
    opts.fdr = -log10(max(pvalues(pFDR <= 0.05)));
    if isempty(opts.fdr); opts.fdr = inf; end

    % h = addLabel(h, [expPvals, obsPvals], opts);
    h = addLabel2(h, [expPvals, obsPvals], opts);
end

% Save the qq plot
if opts.save
    try
        exportgraphics(h.Parent, opts.savename + ".png", 'Resolution', opts.resolution)
    catch
        exportgraphics(h.Parent, matlab.lang.makeValidName(opts.savename) + ".png", 'Resolution', opts.resolution)
    end
end

if opts.savefig
    try
        savefig(h.Parent, opts.savename + ".fig");
    catch
        savefig(h.Parent.Parent, matlab.lang.makeValidName(opts.savename) + ".fig");
    end
end

end % END

%% subfunctions ===========================================================
function n = ppoints(n)
% ppoints function of R https://rdrr.io/r/stats/ppoints.html

if numel(n) > 1 
    n = numel(n);
end

if n <= 10
    a = 3/8;
else
    a = 1/2;
    % a = 0.995; % debug: for top 1% percentile 
end

if n > 0 
    n = ((1:n) - a)./(n + 1 - 2*a);
end

% pt = prctile(n, 1);
% n = linspace(min(n), pt, numel(n));

end

%% ------------------------------------------------------------------------
function lambda = getInflation(pvalues)
% inflationFactor computes genomic inflation factor for a set of input
% p-values
% x = chi2inv(1 - pvalues, 1);
x = pval2chi(pvalues, 1);
lambda = (median(x) / chi2inv(0.5, 1)); % 0.995 for top 1% percentile
end

%% ------------------------------------------------------------------------
function h = addLabel2(h, pos, opts)

pos2 = pos; % for textRepel
opts.labels2 = opts.labels;

% FDR hits in red
idxfdr = pos(:, 2) < opts.fdr;

idx = pos(:, 2) < opts.thresh;
if all(idx & idxfdr) % nothing remains to label
    if opts.showinsig && opts.n > 0
        % show first n insignificant results
        pos = pos(1:opts.n, :);
        opts.labels = opts.labels(1:opts.n); 
        opts.color = repmat("k", opts.n, 1);
    else
        return
    end
else
    
    opts.color = repmat("r", sum(~idxfdr), 1); % FDR
    opts.color(~idx) = "b"; % Bonferroni
    if sum(~(idx & idxfdr)) < opts.n && opts.showinsig % add some insignifcant data to show as well       
        pos = pos(1:opts.n, :);
        opts.labels = opts.labels(1:opts.n); 
        opts.color = [opts.color; repmat("k", opts.n-numel(opts.color), 1)];
    else
        pos(idx, :) = [];
        opts.labels(idx) = [];

        % @13MAY2025: max 20
        n_max = min(20, numel(opts.labels));
        pos = pos(1:n_max, :);
        opts.labels = opts.labels(1:n_max);
    end
end

% for textRepel, set the pos numel to three times its default value to avoid any collisions.
% if numel(opts.labels) < 20
%     fa = 2;
% else
%     fa = 3;
% end

N = min(size(pos2, 1), 2*size(pos, 1));
pos2 = pos2(1:N, :);
opts.labels2 = opts.labels2(1:N);
pos2 = [[0, 0]; pos2];
opts.labels2 = ["NA"; opts.labels2];

% get non-overlapping text labels 
trep = textRepel(pos2(:, 1), pos2(:, 2), opts.labels2, ...
    "MarkerSize", 4, "TextSize", 8, "maxIter", 1e4, "BoxPadding", 0.15);
[~, idx] = ismember(opts.labels, trep.z);
trep = trep(idx, :);

for m = 1:numel(opts.labels)
    ta = annotation("textarrow");
    ta.Parent = h;
    ta.FontSize = 10;
    ta.FontName = "Garamond";
    ta.String = opts.labels(m);
    ta.X = [trep.x_repel(m),pos(m, 1)];
    ta.Y = [trep.y_repel(m),pos(m, 2)];
    % ta.HorizontalAlignment = "center";
    % ta.VerticalAlignment = "cap";
    ta.HeadStyle = "none";
    ta.TextEdgeColor = "none";
    ta.Color = opts.color(m);
end

end

%% ------------------------------------------------------------------------
function h = addLabel(h, pos, opts)

% FDR hits in red
idxfdr = pos(:, 2) < opts.fdr;

idx = pos(:, 2) < opts.thresh;
if all(idx & idxfdr) % nothing remains to label
    if opts.showinsig && opts.n > 0
        % show first n insignificant results
        pos = pos(1:opts.n, :);
        opts.labels = opts.labels(1:opts.n); 
        opts.color = repmat("k", opts.n, 1);
    else
        return
    end
else
    opts.color = repmat("r", sum(~idxfdr), 1); % FDR
    opts.color(~idx) = "b"; % Bonferroni
    if sum(~(idx & idxfdr)) < opts.n && opts.showinsig % add some insignifcant data to show as well       
        pos = pos(1:opts.n, :);
        opts.labels = opts.labels(1:opts.n); 
        opts.color = [opts.color; repmat("k", opts.n-numel(opts.color), 1)];
    else
        pos(idx, :) = [];
        opts.labels(idx) = [];
    end
end

opts.labels = string(opts.labels);
[arr, te] = deal({});
opts.rotang = ones(numel(opts.labels), 1).*-91;
horalgn = repmat("right", numel(opts.labels), 1);
loweridx = pos(:, 1) > pos(:, 2);
opts.rotang(loweridx) = -opts.rotang(loweridx);
horalgn(loweridx) = "left";
% ang = (0);
drawnow
h.Units = 'inches';
% fratio = (h.Position(3) - h.Position(1))/(h.Position(4) - h.Position(2));
h.Units = 'normalized';

opts.labels(ismissing(opts.labels) | opts.labels == "") = "NA";
for i = 1:numel(opts.labels)
    arr{i}= annotation('arrow');
    arr{i}.Parent = h;
    [X, Y] = rotatePoint(pos(i, 1).*opts.factor, pos(i, 2), pos(i, 1), pos(i, 2), opts.rotang(i));
    X(2) = pos(i, 1);
    Y(2) = pos(i, 2);
    
    arr{i}.X = X;
    arr{i}.Y = Y;
    arr{i}.HeadStyle = 'none';
    arr{i}.Color = opts.color(i);
    
    te{i} = text(h, X(1), Y(1), opts.labels(i));
    te{i}.HorizontalAlignment = horalgn(i);
    te{i}.Interpreter = 'none';
    te{i}.Color = opts.color(i);
end

marginCollisionCheck(arr, te, opts);
end

%% ------------------------------------------------------------------------
function marginCollisionCheck(arr, te, opts)

% first fix ax limits
drawnow;
h = te{1}.Parent;
maxY = zeros(numel(te), 1);
for i = 1:numel(te)
    maxY(i) = te{i}.Extent(2) + te{i}.Extent(4);
end
maxY = (max(maxY(:)));
if maxY > h.YLim(2)
    h.YLim(2) = maxY;
end

% traversed/checked textboxes should have fixed positions and other
% textboxes cannot overlap these
opts.fixedpos = []; % fixed positions 
angcoeffs = -sign(opts.rotang);
for i = 1:numel(te)-1

    % moving direction for next text obj in relation to the current one
    opts.angcoeff = angcoeffs(i);

    % rotate text obj i untill it doesn't pass obj i+1, if so scale arrow
    % length instead of rotation.
    opts.nextcoeff = angcoeffs(i + 1);
    opts.nextx = arr{i + 1}.X; 
    opts.nexty = arr{i + 1}.Y; 
    
    if i > 1
        % in combbination with the condition above, the arrow should
        % ideally be somewhere between next/previous arrows. By having the
        % slope of previous arrow, we can recalculate the angle when
        % rescaling the length of current arrow (recalculation is necessary
        % because 'angs' returned from collisionCheck only corresponds to
        % unscaled arrow).
        opts.prevb = diff(arr{i-1}.Y)/diff(arr{i-1}.X);
        % if both text objs are not on the same side, reset ang to zero
        % (since they current text won't overlap with the previous one
        % anymore)
        opts.prevcoeff = angcoeffs(i - 1); 
    end

    angs = zeros(numel(i+1:numel(te)), 1);
    for j = i+1:numel(te)
        isflippedi = diff(arr{i}.Y) > 0; % to check if text obj pos is flipped with regards the to diagonal line
        isflippedj = diff(arr{j}.Y) > 0;

        [arr{i}, te{i}, arr{j}, te{j}, angs(j)] = ...
            collisionCheck(arr{i}, te{i}, arr{j}, te{j}, opts);

        isflippedi = (isflippedi + (diff(arr{i}.Y) > 0)) == 1;
        isflippedj = (isflippedj + (diff(arr{j}.Y) > 0)) == 1;
        if isflippedi
            angcoeffs(i) = -angcoeffs(i);
        end
        if isflippedj
            angcoeffs(j) = -angcoeffs(j);
        end
    end
    
    [~, idx] = max(abs(angs));
    opts.ang = angs(idx);

    drawnow;
    te{i}.Units = 'normalized';
    opts.fixedpos = [opts.fixedpos; te{i}.Extent];
end

drawnow;
for i = 1:numel(te)
    te{i}.Units = 'data';
end

end

%% ------------------------------------------------------------------------
function [arr2, te2, arr1, te1, ang] = collisionCheck(arr2, te2, arr1, te1, opts)

if nargin < 5
    opts.fixedpos = [];
end
opts.factor = 1;

[arr1, te1] = marginCheck(arr1, te1, false);
[arr2, te2] = marginCheck(arr2, te2, false);

% check overlapping of two text/arrow annotations
collide = collisionCond(te1, te2, true, opts.margin);

if ~collide
    ang = 0;
    return
end

% current textbox should be to te left of previous text box (more negative angle)
step = opts.angcoeff;
if isfield(opts, 'ang') && opts.ang ~= 0 
    ang = opts.ang;  
else
    ang = step; % counterclockwise(below)/clockwise(above)
    step = -step;  
end
opts.nextcoeff = opts.nextcoeff*opts.angcoeff;

if isfield(opts, 'prevcoeff') % to see if ang should be reset
    if opts.prevcoeff*opts.angcoeff < 0
        ang = opts.angcoeff;
        step = -step;
    end
end
% resetAngle = true;

Xorig = arr2.X; Yorig = arr2.Y;

% for new angle based on new scaled arrow coordinates
angfun = @(a, X1new, Y1new, b, c) (((Xorig(1)*c - Xorig(2))*cosd(a) +...
    (Yorig(1)*c - Yorig(2))*sind(a) + Xorig(2)) - X1new).*b + Y1new ...
    + (Xorig(1)*c - Xorig(2))*sind(a) - (Yorig(1)*c - Yorig(2))*cosd(a) - Yorig(2);
syms a
angfun_solve = @(X1new, Y1new, b, c)(((Xorig(1)*c - Xorig(2))*cosd(a) + ...
    (Yorig(1)*c - Yorig(2))*sind(a) + Xorig(2)) - X1new).*b + ...
    Y1new == -(Xorig(1)*c - Xorig(2))*sind(a) + (Yorig(1)*c -...
    Yorig(2))*cosd(a) + Yorig(2);

angcnt = 1;

% for maginCheck
if opts.angcoeff > 0 % above diagonal
    cwflag = false;
else
    cwflag = true;
end

while collide % change position of latest text/arrow 
    
    % to check if text obj pos is flipped with regards the to diagonal line
    isflipped2 = diff(arr2.Y) > 0; 

    drawnow;
    [arr2.X(1), arr2.Y(1)] = rotatePoint(opts.factor*Xorig(1), opts.factor*Yorig(1), Xorig(2), Yorig(2), ang);
    te2.Units = 'data';
    te2.Position(1:2) = [arr2.X(1), arr2.Y(1)];
    
    % if flipped, then change also horizontal place of text tag
    isflipped2 = (isflipped2 + (diff(arr2.Y) > 0)) == 1;
    if isflipped2
        te2.HorizontalAlignment = setdiff(["left", "right"], te2.HorizontalAlignment);
        opts.angcoeff = -opts.angcoeff;
    end

    % [arr1, te1] = marginCheck(arr1, te1, true);
    [arr2, te2, ~, hadoverlap] = marginCheck(arr2, te2, cwflag);

    if hadoverlap % next time change cw argument of marginCheck to flip if necessary
        [arr2, te2] = flipSide(arr2, te2);
        [arr2, te2,] = marginCheck(arr2, te2, ~cwflag);
    end

    % if new arr2.X(1) is less than next text obj (outer loop in
    % marginCollisionCheck) arrow X, then stop rotating. This is because
    % the rotation creates messy cross lines. For points above diagonal red
    % line, the condition to check is X_i(1) > min(X_i+1), and for points
    % below, it is X_i(1) < min(X_i+1). If this condition is met, the
    % increase the arrow length (same as the next condition)
    crosslineCond = false;
    if isfield(opts, 'ang') && opts.ang ~= 0 && opts.nextcoeff > 0 % check only objects on the same side of red line
        crosslineCond = polyxpoly(arr2.X, arr2.Y, arr1.X, arr1.Y);
        if ~isempty(crosslineCond); crosslineCond = true; else; crosslineCond = false; end
    end

    angcnt = angcnt + 1;
    if angcnt > opts.maxAng % don't change angle anymore, instead change size
        maxangCond = true;
    else
        maxangCond = false;
    end
    
    % if marginCheck modifies te2, this means it passes axis limits, and
    % this is a trap because rotatePint again pushes it away of limits. To
    % tackle this issue, we increase the arrow length
    drawnow;
    if abs(arr2.X(1) - te2.Position(1)) > eps || ...
            abs(arr2.Y(1) - te2.Position(2)) > eps || crosslineCond || maxangCond
        opts.factor = opts.factor + opts.factorScale;
        
        if isfield(opts, 'prevb') && (crosslineCond || maxangCond)
            % reset rotation angle (only once): since the length of arrow is
            % changing, we need to first estimate the size of new arrow and map
            % it to have the same slope as the previous one, and then estimate
            % the resetting angle.
            d = sqrt(diff(Xorig)^2 + diff(Xorig)^2); % length of original arrow (unscaled)
            d = 4*(d*opts.factor - 0.75*d); % scaled length of arrow (approx.);
            % solve for mapped points (same slope as previous arrow, and new
            % scaled length).
            if opts.prevb == inf
                opts.prevb = 1e3;
            elseif opts.prevb == -inf
                opts.prevb = -1e3;
            end
            eq1 = @(x,y) (Yorig(2)-y)./(Xorig(2)-x) - opts.prevb;
            eq2 = @(x,y) sqrt((Yorig(2)-y)^2 + (Xorig(2)-x)^2) - d;
            eqn = @(xv) [eq1(xv(1),xv(2)); eq2(xv(1),xv(2))];
            op = optimoptions('fsolve', 'Display', 'off');
            newpoints = fsolve(eqn, [Xorig(1); Yorig(1)], op);
            % two solutions exist: above and below center. If current
            % solution doesn't correspond to current point, then generate a
            % new one
            if (step < 0 && newpoints(2) > Yorig(1))...
                    || step > 0 && newpoints(2) < Yorig(1) 
                [newpoints(1), newpoints(2)] = rotatePoint(newpoints(1), ...
                    newpoints(2), Xorig(2), Yorig(2),-180);
            end
            step = opts.angcoeff;
            ang = double(vpasolve(angfun_solve(newpoints(1), newpoints(2), ...
                opts.prevb, opts.factor), [-360, 0]));
            
            % check if ang rotation flips the side of text obj
            isflipped =  Yorig(1) > Yorig(2);
            [~, yr] = rotatePoint(opts.factor*Xorig(1), opts.factor*Yorig(1), Xorig(2), Yorig(2), ang);
            isflipped = (isflipped + (yr > Yorig(2)) > 0) == 1;
            if isflipped
                ang = double(vpasolve(angfun_solve(newpoints(1), newpoints(2), ...
                opts.prevb, opts.factor), [0, 360]));
            end
            ang = ang + step; % to cancel out change at the end of while loop
            
            % reset the counter (change the ang again)
            if maxangCond && opts.resetAngle
                angcnt = 1;
                step = step./2;
            end 
            
%         else
%             if resetAngle
%                 if isfield(opts, 'ang') && opts.ang ~= 0 
%                     step = opts.angcoeff;
%                     ang = opts.ang + step;
%                 else
%                     step = 1;
%                     ang = 0;
%                 end
%                 resetAngle = false;
%             end
        end
    end
    
    % forbiddedn positions: those textboxes which have already been checked/traversed 
    if ~isempty(opts.fixedpos)
        collide = collisionCond(te2, opts.fixedpos, false, opts.margin);
    else
        collide = false;
    end

    collide = collisionCond(te1, te2, true, opts.margin) | collide | crosslineCond;
    if ~collide
        break
    end
    ang = ang - step;
end

% move text label to lef/right of arrow's origin
if (arr2.X(1) > arr2.X(2))
   te2.HorizontalAlignment = "left";
else
    te2.HorizontalAlignment = "right";
end
if ~isempty(opts.fixedpos)
    collide = collisionCond(te2, opts.fixedpos, false, opts.margin);
else
    collide = false;
end
collide = collisionCond(te1, te2, true, opts.margin) | collide;
if collide % reset text label horizontal alignment
    te2.HorizontalAlignment = setdiff(["left", "right"], te2.HorizontalAlignment);
end

% compute the equivalent angle (i.e opts.factor = 1) for the next text obj
if opts.factor > 1 
    Xnew = arr2.X;
    Ynew = arr2.Y;

    % current scaled (opts.factor) line equation
    dXnew = diff(Xnew); if ~dXnew; dXnew = 0.1; end
    b = diff(Ynew)/dXnew;
    % eq = @(x) (x - Xnew(1)).*b + Ynew(1);

    % solve for angle with Xnew/Ynew as Xrot/Yrot
    ang = fzero(@(x)angfun(x, Xnew(1), Ynew(1), b, 1), -2);
end

end

%% ------------------------------------------------------------------------
function [arr, te, ang, hadoverlap] = marginCheck(arr, te, cw)
drawnow;
te.Units = 'normalized';
ext = te.Extent;
if ((ext(1) + ext(3)) >= 1.01) || ext(1) <= 0 % cuts X-axis
    xax = true;
else
    xax = false;
end
if (ext(2) + ext(4)) >= 1.01 || ext(2) <= 0 % cuts Y-axis
    yax = true;
else
    yax = false;
end
te.Units = 'data';

step = 2;
if cw
    ang = step; % clockwise
else
    ang = -step; % counterclockwise
end

hadoverlap = false;
while any([xax, yax])
    hadoverlap = true;
    [arr.X(1), arr.Y(1)] = rotatePoint(arr.X(1), arr.Y(1), arr.X(2), arr.Y(2), ang);
    te.Position(1:2) = [arr.X(1), arr.Y(1)];

    drawnow;
    if diff(arr.X) > 0 && diff(arr.Y) < 0 
        te.HorizontalAlignment = 'right';
    end

    drawnow;
    te.Units = 'normalized';
    ext = te.Extent;
    drawnow;
    if (ext(1) + ext(3)) < 1 && (ext(2) + ext(4)) < 1 && ext(1) > 0 && ext(2) > 0
        break
%     else % center adjust the text and check again
%         drawnow
%         hal = te.HorizontalAlignment;
%         te.HorizontalAlignment = 'center';
%         ext = te.Extent;
%         if (ext(1) + ext(3)) < 1 && (ext(2) + ext(4)) < 1 
%             break
%         else
%             drawnow
%             te.HorizontalAlignment = hal;
%         end
    end
    te.Units = 'data';
    if cw
        ang = ang + step;
    else % ccw
        ang = ang - step;
    end

end
end

%% ------------------------------------------------------------------------
function [arr, te] = flipSide(arr, te)
drawnow;
[arr.X(1), arr.Y(1)] = rotatePoint(arr.X(1), arr.Y(1), arr.X(2), arr.Y(2), 180);
te.Units = 'data';
te.Position(1:2) = [arr.X(1), arr.Y(1)];
te.Units = 'normalized';
te.HorizontalAlignment = setdiff(["left", "right"], te.HorizontalAlignment);
end

%%
function [Xrot, Yrot] = rotatePoint(X, Y, Xc, Yc, ang)
Xrot =  (X-Xc)*cosd(ang) + (Y-Yc)*sind(ang) + Xc;
Yrot = -(X-Xc)*sind(ang) + (Y-Yc)*cosd(ang) + Yc;
end

%%
function collide = collisionCond(te1, te2, check, thresh)

if nargin < 4
    thresh = 0.01;
end
if nargin < 3
    check = true;
end
drawnow;
te1.Units = 'normalized';

if check 
    te2.Units = 'normalized';
    ext.one = te1.Extent; ext.two = te2.Extent;
else
    ext.one = te1.Extent; ext.two = te2; % te2 is Extent property (can be matrix for different text boxes)
end

for i = 1:size(ext.two, 1)
    if ext.one(1) < ext.two(i, 1) % smaller X/Y values must be used
        collide = ext.one(1) + ext.one(3) >= ext.two(i, 1) + thresh;
    else
        collide = ext.two(i, 1) + ext.two(i, 3) >= ext.one(1) + thresh;
    end
    
    if ext.one(2) < ext.two(i, 2)
        collide = collide & (ext.one(2) + ext.one(4) >= ext.two(i, 2) + thresh);
    else
        collide = collide & (ext.two(i, 2) + ext.two(i, 4) >= ext.one(2) + thresh);
    end

    if collide, break, end % one (overlap) true is enough
end

end
