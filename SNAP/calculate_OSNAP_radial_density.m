% WARNING: technically overestimates by attributing points beyond the 9th
% ring as included in 10th region
function [radial_density, area] = calculate_OSNAP_radial_density(points, axes_lengths, inc)
    n_r = floor(1/inc);
    radial_density = zeros(n_r,1);
    rad = inc:inc:1;
    % first innermost ellipse
    radial_density(1) = sum(is_in_ellipse(points, axes_lengths*inc))/ellipse_area(axes_lengths*inc);
    % next n rings
    for i=2:n_r-1
        n_points_outer = sum(is_in_ellipse(points, axes_lengths*rad(i)));
        n_points_inner = sum(is_in_ellipse(points, axes_lengths*rad(i-1)));
        area_outer = ellipse_area(axes_lengths*rad(i));
        area_inner = ellipse_area(axes_lengths*rad(i-1));
        radial_density(i) = (n_points_outer - n_points_inner)/(area_outer-area_inner);
    end
    % final ring includes points outside as well
    n_points = size(points,1)-n_points_outer;
    area = ellipse_area(axes_lengths)-area_outer;
    radial_density(n_r) = n_points/area;
    clearvars n_r rad n_points_outer n_points_inner area_outer area_inner n_point area
end

function area = ellipse_area(axes_lengths)
    area = axes_lengths(1)*axes_lengths(2)*pi/4;
end

function idx = is_in_ellipse(points, axes_lengths)
    idx = ((points(:,1)/(axes_lengths(1)/2)).^2) + ((points(:,2)/(axes_lengths(2)/2)).^2) - 1 < 0;
end

