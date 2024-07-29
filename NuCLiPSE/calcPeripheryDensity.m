function [interiorDensity, peripheryDensity]  = calcPeripheryDensity(points,polygon,peripheryThresh)
    polygon = translate(polygon,-(max(polygon.Vertices) - range(polygon.Vertices)/2)); % ensure that the geometric center of the polygon (not centroid) is at the origin
    interior = scale(polygon,1-peripheryThresh);
    interiorN = sum(insidepoly(points(:,1),points(:,2),interior.Vertices(:,1),interior.Vertices(:,2)));
    interiorDensity = interiorN/area(interior);
    periphery =  addboundary(polygon,interior.Vertices(:,1),interior.Vertices(:,2));
    peripheryDensity = (size(points,1) - interiorN)/area(periphery);
end