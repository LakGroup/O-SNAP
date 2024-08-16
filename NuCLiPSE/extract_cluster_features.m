function extract_cluster_features(filepath)
    if ~has_variables(filepath,{'cluster_n_locs','cluster_area','cluster_density','cluster_gyration_R'})
        vars_to_load = {'x','y','voronoi_areas_all',...
            'area_thresholds','min_number_of_localizations','clusters'};
        data = load_variables(filepath, vars_to_load);
        n_thresholds = length(data.area_thresholds);
        data.cluster_n_locs = cell(n_thresholds,1);
        data.cluster_area = cell(n_thresholds,1);
        data.cluster_density = cell(n_thresholds,1);
        data.cluster_gyration_R = cell(n_thresholds,1);
        for i=1:n_thresholds
            clusters = data.clusters(:,i);
            n_cluster = max(clusters);
            filtered_points = [data.x data.y data.voronoi_areas_all clusters];
            filtered_points = filtered_points(clusters > 0,:);
            data.cluster_n_locs{i} = arrayfun(@(x) sum(filtered_points(:,4) == x), 1:n_cluster,'uni',1)';
            data.cluster_area{i} = arrayfun(@(x) sum(filtered_points(filtered_points(:,4) == x,3),'omitmissing'), 1:n_cluster,'uni',1)';
            data.cluster_density{i} = arrayfun(@(x) data.cluster_n_locs{i}(x)/data.cluster_area{i}(x), 1:n_cluster,'uni',1)';
            data.cluster_gyration_R{i} = arrayfun(@(x) sqrt(...
                mean((filtered_points(filtered_points(:,4) == x,1)...
                -mean(filtered_points(filtered_points(:,4) == x,1))).^2 ...
                + (filtered_points(filtered_points(:,4) == x,2) ...
                -mean(filtered_points(filtered_points(:,4) == x,2))).^2)), 1:n_cluster,'uni',1);
        end
        clearvars -except filepath data vars_to_load
        data = rmfield(data,vars_to_load);
        save_voronoi_data(filepath, data);
    end
end

