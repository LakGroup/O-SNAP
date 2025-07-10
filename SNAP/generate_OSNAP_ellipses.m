% -------------------------------------------------------------------------
% generate_OSNAP_ellipses.m
% -------------------------------------------------------------------------
% Creates the elliptical polyshapes that are used to visualize and
% facilitate radial modeling and analyis of samples (nucleus boundaries) as
% concentric ellipses.
%
% Example on how to use it:
%   p = generate_OSNAP_ellipses(1.6E4,1.8E2,0.1)
% -------------------------------------------------------------------------
% Input:
%   major_axis: Length of major axis of polygon [nm]
%   minor_axis: Length of minor axis of polygon [nm]
%   inc: Percentage of the major/minor axis length to subdivide the
%        ellipses that define the bounds of the rings. The number of rings
%        is N_r = floor(1/inc)
% Output:
%   p: Series of polyshapes that approximate ellipses, acting as boundaries
%      for subsequent radial anaylsis
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
function p = generate_OSNAP_ellipses(major_axis, minor_axis, inc)
    inc_t = 30;
    t = 2*pi/inc_t:2*pi/inc_t:2*pi; % parametric variable for ellipse
    rad = inc:inc:1; % size of radius
    x = (major_axis/2)*cos(t);
    y = (minor_axis/2)*sin(t);
    p = repmat(polyshape,numel(rad),1);
    for r=1:numel(rad)
        % define outer boundary
        x_out = rad(r)*x;
        y_out = rad(r)*y;
        % generate ring
        if r == 1
            p(r) = polyshape(x_out,y_out);
        else
            % define inner boundary
            x_in = (rad(r)-inc)*x;
            y_in = (rad(r)-inc)*y;
            p(r) = polyshape({x_in, x_out},{y_in,y_out});
        end
    end
end


