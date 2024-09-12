function [interior_density, periphery_density]  = calculate_periphery_density(points,polygon,periphery_thresh)
    polygon = translate(polygon,-(max(polygon.Vertices) - range(polygon.Vertices)/2)); % ensure that the geometric center of the polygon (not centroid) is at the origin
    interior = scale(polygon,1-periphery_thresh);
    interior_n = sum(insidepoly(points(:,1),points(:,2),interior.Vertices(:,1),interior.Vertices(:,2)));
    interior_density = interior_n/area(interior);
    periphery =  addboundary(polygon,interior.Vertices(:,1),interior.Vertices(:,2));
    periphery_density = (size(points,1) - interior_n)/area(periphery);
end

