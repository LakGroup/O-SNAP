function extract_cluster_features(filepath)
        vars_to_load = {'x','y','voronoi_areas_all',...
            'area_thresholds','min_number_of_localizations','voronoi_clusters'};
        data = load_variables(filepath, vars_to_load);
        n_thresholds = length(data.area_thresholds);
        data.voronoi_cluster_n_locs = cell(n_thresholds,1);
        data.voronoi_cluster_radius = cell(n_thresholds,1);
        data.voronoi_cluster_density = cell(n_thresholds,1);
        data.voronoi_cluster_gyration_radius = cell(n_thresholds,1);
        for j=1:n_thresholds
            voronoi_clusters = data.voronoi_clusters(:,j);
            n_cluster = max(voronoi_clusters);
            points = [data.x data.y data.voronoi_areas_all voronoi_clusters];
            points = points(voronoi_clusters > 0,:); % filter out voronoi_cluster=0 that represents non-clustered data
            points_by_cluster = splitapply(@(x){(x)},points(:,1:3),findgroups(points(:,4))); % Split the clusters into the corresponding groups
            data.voronoi_cluster_n_locs{j} = cellfun('size',points_by_cluster,1); % use cellfun+length, faster from backwards compatibility
            data.voronoi_cluster_radius{j} = zeros(n_cluster,1);
            data.voronoi_cluster_density{j} = zeros(n_cluster,1);
            data.voronoi_cluster_gyration_radius{j} = zeros(n_cluster,1);
            for i=1:n_cluster
                data.voronoi_cluster_radius{j}(i) = sqrt(sum(points_by_cluster{i}(:,3))/pi);
                data.voronoi_cluster_density{j}(i) = data.cluster_n_locs{j}(i)./data.voronoi_cluster_radius{j}(i);
                data.voronoi_cluster_gyration_radius{j}(i) = sqrt(mean(...
                    (points_by_cluster{i}(:,1)-mean(points_by_cluster{i}(:,1))).^2 ...
                    + (points_by_cluster{i}(:,2)-mean(points_by_cluster{i}(:,2))).^2));
            end
        end
        clearvars -except filepath data vars_to_load
        data = rmfield(data,vars_to_load);
        save_voronoi_data(filepath, data);
end


