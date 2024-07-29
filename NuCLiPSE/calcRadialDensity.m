% WARNING: technically overestimates by attributing points beyond the
function [radialDensity, area] = calcRadialDensity(points, axesLengths, inc)
    n_r = floor(1/inc);
    radialDensity = zeros(n_r,1);
    rad = inc:inc:1;
    % first innermost ellipse
    radialDensity(1) = sum(isInEllipse(points, axesLengths*inc))/ellipseArea(axesLengths*inc);
    % next n rings
    for i=2:n_r-1
        nPoints_outer = sum(isInEllipse(points, axesLengths*rad(i)));
        nPoints_inner = sum(isInEllipse(points, axesLengths*rad(i-1)));
        area_outer = ellipseArea(axesLengths*rad(i));
        area_inner = ellipseArea(axesLengths*rad(i-1));
        radialDensity(i) = (nPoints_outer - nPoints_inner)/(area_outer-area_inner);
    end
    % final ring includes points outside as well
    nPoints = size(points,1)-nPoints_outer;
    area = ellipseArea(axesLengths)-area_outer;
    radialDensity(n_r) = nPoints/area;
    clearvars n_r rad nPoints_outer nPoints_inner area_outer area_inner nPoint area
end

function area = ellipseArea(axesLengths)
    area = axesLengths(1)*axesLengths(2)*pi/4;
end

function idx = isInEllipse(points, axesLengths)
    idx = ((points(:,1)/(axesLengths(1)/2)).^2) + ((points(:,2)/(axesLengths(2)/2)).^2) - 1 < 0;
end