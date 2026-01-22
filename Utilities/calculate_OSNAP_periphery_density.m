% -------------------------------------------------------------------------
% calculate_OSNAP_periphery_density.m
% -------------------------------------------------------------------------
% Calculates the density for number of points located at the periphery and
% interior regions where the periphery is defined using a percent area
% threshold cutoff. 
% *THIS CODE ASSUMES THAT ALL POINTS FALL WITHIN THE POLGON
%
% Example on how to use it:
%   [interior_density, periphery_density] = calculate_OSNAP_periphery_density(points,polygon,periphery_thresh)
% -------------------------------------------------------------------------
% Input:
%   points: [Nx2] array of coordinates to determine whether they fall
%           inside the interior or periphery region.
%   polygon: Polyshape object of the shape of interest (e.g. nucleus
%            boundary)
%   periphery_thresh: Cutoff for what percentage of the polygon to consider
%                     as the periphery. Ex. 0.15 results in a periphery
%                     region with an area equivalent to 15% of the original
%                     polygon shape
% Output:
%   interior_density: Number of points per unit area for the interior
%                     region
%   periphery_density: Number of points per unit area for the periphery
%                     region
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function [interior_density, periphery_density]  = calculate_OSNAP_periphery_density(points,polygon,periphery_thresh)
    polygon = translate(polygon,-(max(polygon.Vertices) - range(polygon.Vertices)/2)); % ensure that the geometric center of the polygon (not centroid) is at the origin
    interior = scale(polygon,1-periphery_thresh);
    % interior_n = sum(insidepoly(points(:,1),points(:,2),interior.Vertices(:,1),interior.Vertices(:,2)));
    interior_n = sum(inpoly(points,interior.Vertices));
    interior_density = interior_n/area(interior);
    periphery =  addboundary(polygon,interior.Vertices(:,1),interior.Vertices(:,2));
    periphery_density = (size(points,1) - interior_n)/area(periphery);
end


