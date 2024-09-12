function plot_radial_density_batch(file_path, max_radial_density, max_radial_hetero_density)
    vars_to_load = {'locs_norm','x_length','y_length',...
        'components','boundary'...
        'radial_density','radial_hetero_density',...
        'hetero_center','voronoi_areas_all'};
    data = load_variables(file_path,vars_to_load);
    [save_dir,name,~] = fileparts(file_path);
    save_dir = replace(save_dir,"voronoi_data","nuclei_images");
    %% plot
    plot_file = fullfile(save_dir,name+"_radial_density.png");
    % if ~exist(plot_file,'file')
        plot_radial_density(data.radial_density, data.locs_norm, data.voronoi_areas_all, data.x_length, data.y_length, max_radial_density, 1, plot_file);
    % end
    
    plot_file = fullfile(save_dir,name+"_radial_hetero_density.png");
    % if ~exist(plot_file,'file')
        plot_radial_density(data.radial_hetero_density,...
            normalize_points(data.hetero_center,data.components,data.boundary),...
            zeros(size(data.hetero_center,1),1), data.x_length, data.y_length, max_radial_hetero_density, 10, plot_file);
    % end
end

function points_norm = normalize_points(points,components,boundary)
    rotated_boundary = boundary*components;
    points_norm = points*components - min(rotated_boundary) - range(rotated_boundary)/2;
end

