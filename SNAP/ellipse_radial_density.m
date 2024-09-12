function p = ellipse_radial_density(major_axis, minor_axis, inc)
    t = inc:inc:2*pi; % parametric variable for ellipse
    rad = inc:inc:1; % size of radius
    x = (major_axis/2)*cos(t); y = (minor_axis/2)*sin(t);
    p = repmat(polyshape,length(rad),1);
    for r=1:length(rad)
        % define outer boundary
        x_out = rad(r)*x; y_out = rad(r)*y;
        % generate ring
        if r == 1
            p(r) = polyshape(x_out,y_out);
        else
            % define inner boundary
            x_in = (rad(r)-inc)*x; y_in = (rad(r)-inc)*y;
            p(r) = polyshape({x_in, x_out},{y_in,y_out});
        end
    end
end

