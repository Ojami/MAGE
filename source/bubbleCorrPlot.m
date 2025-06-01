function bubbleCorrPlot(r, opts)

arguments
    r {mustBeA(r, 'double')} % correlation matrix
    opts.labels {mustBeText} % corresponds to matrix r
end

disp('WARNING: under development')

f = figure;
ax = axes(f);

% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(ax, x(:), y(:), 400.*abs(r(:)), r(:), 'filled')

% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines

% show labels
if isfield(opts, 'labels')
    text(ax, 1:n, (n:-1:1) + 0.5, opts.labels, 'HorizontalAlignment', 'right')
    text(ax, (1:n) + 0.5, repmat(n + 1, n, 1), opts.labels, ...
        'HorizontalAlignment', 'right', 'Rotation', 270)
end

colorbar(ax);
ax.Visible = 'off';
ax.Position(4) = ax.Position(4)*0.9;
axis(ax, 'equal')
colormap('jet')

end