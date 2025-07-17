% -------------------------------------------------------------------------
% calculate_OSNAP_radial_density.m
% -------------------------------------------------------------------------
% Calculates the density for number of points located in a series of
% rings defined by evenly segmenting the major/minor axes of a shape and
% fitting elliptical boundaries.
% *ANY POINTS BEYOND THE PENULTIMATE ELLIPSE ARE ASSIGNED TO THE OUTERMOST
%  REGION
% Example on how to use it:
%   [radial_density, area] = calculate_OSNAP_radial_density(...
%                                       locs, 
%                                       [1.6E4,1.8E2], 
%                                       0.1)
% -------------------------------------------------------------------------
% Input:
%   points: [Nx2] array of coordinates to determine whether they fall
%           inside the interior or periphery region.
%   axes_lengths: [1x2] array consisting of the  [major axis length,
%                 minor axis length] of the original shape of interest
%   inc: Percentage of the major/minor axis length to subdivide the
%        ellipses that define the bounds of the rings. The number of rings
%        is N_r = floor(1/inc)
% Output:
%   radial_density: [N_rx1] array where each value corresponds to the
%                   calculated density for rings 1 to N_r
%   area: [N_rx1] where each value is the area for rings 1 to N_r
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% ------------------------------------------------------------------------- 
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


