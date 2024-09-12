function [threshold, storm_data] = get_percentile_threshold(locs,pct)
    x = locs(:,1);
    y = locs(:,2);
    
    dt = delaunay_triangulation(x,y);
    [vertices,connections] = voronoi_diagram(dt);
    voronoi_cells = cellfun(@(x) vertices(x,:),connections,'uni',0);
    voronoi_areas = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'uni',0);
    voronoi_areas = vertcat(voronoi_areas{:});
    
    % cal density
    density = 1./voronoi_areas;
    storm_data = [locs(:,1:2),density];
    threshold = prctile(storm_data(:,3),pct);
%     fprintf('Now processing %s -- Data Size: %d \n', data{cell_idx}.name, length(locs(:,1)))
    % fprintf('density threshold is %d\n', threshold)
end

