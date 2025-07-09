function plot_OSNAP_ellipses_locs_and_clusters(file_path, max_radial_loc_density, max_radial_dbscan_cluster_density,options)
arguments
    file_path string
    max_radial_loc_density double
    max_radial_dbscan_cluster_density double
    options.overwrite = true;
end
    vars_to_load = {'locs_norm','x_length','y_length',...
        'components','boundary'...
        'radial_loc_density','radial_dbscan_cluster_density',...
        'interior_dbscan_cluster_center','periphery_dbscan_cluster_center',...
        'voronoi_areas_all'};
    data = load_variables_OSNAP(file_path,vars_to_load);
    [save_dir,name,~] = fileparts(file_path);
    save_dir = replace(save_dir,"SNAP_nucleus_data","SNAP_radial_density_images");
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
    %% plot
    plot_file = fullfile(save_dir,name+"_radial_loc_density.png");
    if ~exist(plot_file,'file') || options.overwrite
        plot_OSNAP_ellipses(data.radial_loc_density,...
            data.locs_norm,...
            data.voronoi_areas_all,...
            data.x_length, data.y_length, max_radial_loc_density, false, plot_file);
    end
    plot_file = fullfile(save_dir,name+"_radial_dbscan_cluster_density.png");
    if ~exist(plot_file,'file') || options.overwrite
        locs = [data.periphery_dbscan_cluster_center; data.interior_dbscan_cluster_center];
        plot_OSNAP_ellipses(data.radial_dbscan_cluster_density,...
            normalize_points(locs,data.components,data.boundary),...
            zeros(size(locs,1),1),...
            data.x_length, data.y_length, max_radial_dbscan_cluster_density, true, plot_file);
    end
end

function points_norm = normalize_points(points,components,boundary)
    rotated_boundary = boundary*components;
    points_norm = points*components - min(rotated_boundary) - range(rotated_boundary)/2;
end

