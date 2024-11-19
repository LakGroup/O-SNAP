function extract_cluster_features(filepath)
    % if ~has_variables(filepath,{'cluster_n_locs','cluster_area','cluster_density','cluster_gyration_R'})
        vars_to_load = {'x','y','voronoi_areas_all',...
            'area_thresholds','min_number_of_localizations','clusters'};
        data = load_variables(filepath, vars_to_load);
        n_thresholds = length(data.area_thresholds);
        data.cluster_n_locs = cell(n_thresholds,1);
        data.cluster_area = cell(n_thresholds,1);
        data.cluster_density = cell(n_thresholds,1);
        data.cluster_gyration_R = cell(n_thresholds,1);
        for j=1:n_thresholds
            clusters = data.clusters(:,j);
            n_cluster = max(clusters);
            points = [data.x data.y data.voronoi_areas_all clusters];
            points = points(clusters > 0,:); % filter out cluster=0 that represents non-clustered data
            points_by_cluster = splitapply(@(x){(x)},points(:,1:3),findgroups(points(:,4))); % Split the clusters into the corresponding groups
            data.cluster_n_locs{j} = cellfun('size',points_by_cluster,1); % use cellfun+length, faster from backwards compatibility
            data.cluster_area{j} = zeros(n_cluster,1);
            data.cluster_density{j} = zeros(n_cluster,1);
            data.cluster_gyration_R{j} = zeros(n_cluster,1);
            for i=1:n_cluster
                data.cluster_area{j}(i) = sum(points_by_cluster{i}(:,3));
                data.cluster_density{j}(i) = data.cluster_n_locs{j}(i)./data.cluster_area{j}(i);
                data.cluster_gyration_R{j}(i) = sqrt(mean(...
                    (points_by_cluster{i}(:,1)-mean(points_by_cluster{i}(:,1))).^2 ...
                    + (points_by_cluster{i}(:,2)-mean(points_by_cluster{i}(:,2))).^2));
            end
        end
        clearvars -except filepath data vars_to_load
        data = rmfield(data,vars_to_load);
        save_voronoi_data(filepath, data);
    % end
end


