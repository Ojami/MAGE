function coords = getMarkerBounds(ax)
% gets marker width in data unit from axes 'ax' handle.

% position in points
drawnow;
units = ax(1).Parent.Units;
ax(1).Parent.Units = "points";
axpos = ax(1).Parent.Position;
ax(1).Parent.Units = units;
Xlim = ax(1).Parent.XLim;
Ylim = ax(1).Parent.YLim;

for k = 1:numel(ax)

    if ~isUnderlyingType(ax(k), 'matlab.graphics.chart.primitive.Line')
        continue
    end

    if ax(k).Visible == "off"
        continue
    end
    
    % marker center
    x = ax(k).XData;
    y = ax(k).YData;
    sz = ax(k).MarkerSize;
    w = sz*diff(Xlim)/axpos(3);
    h = sz*diff(Ylim)/axpos(4);

    if any(strcmp(ax(k).Marker, ["square", "s"]))
        % adjust the boundaries: circile --> square
        h = h*sqrt(2)/2;
        w = w*sqrt(2)/2;
    end

    coords.w(k, 1) = w;
    coords.h(k, 1) = h;
    coords.x(k, 1:2) = [x - w/2, x + w/2];
    coords.y(k, 1:2) = [y - h/2, y + h/2];
end

end % END