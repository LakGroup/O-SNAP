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
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
%%
function [interior_density, periphery_density]  = calculate_OSNAP_periphery_density(points,polygon,periphery_thresh)
    % Ensure that the geometric center of the polygon (not centroid) is at the origin
    polygon = translate(polygon,-(max(polygon.Vertices) - range(polygon.Vertices)/2));
    interior = scale(polygon,1-periphery_thresh);
    if size(points,1) > 0
        % Calculate interior density
        interior_n = sum(inpoly(points,interior.Vertices));
        interior_density = interior_n/area(interior);
        % Calculate periphery density
        periphery =  addboundary(polygon,interior.Vertices(:,1),interior.Vertices(:,2));
        periphery_density = (size(points,1) - interior_n)/area(periphery);
    else
        interior_density = 0;
        periphery_density = 0;
    end
end