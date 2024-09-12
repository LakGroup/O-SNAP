function generate_voronoi_clusters(filepath,options)
arguments
    filepath string
    options.area_threshold_arr double;
    options.min_number_of_localizations_arr double;
end
    if ~has_variables(filepath,{'area_thresholds','min_number_of_localizations','clusters'})
        vars_to_load = {'voronoi_areas_all','voronoi_neighbors'};
        data = load_variables(filepath,vars_to_load);
        n_area_threshold = length(options.area_threshold_arr);
        n_min_number_of_localizations = length(options.min_number_of_localizations_arr);
        data.area_thresholds = options.area_threshold_arr;
        data.min_number_of_localizations = options.min_number_of_localizations_arr;
        data.clusters = zeros(length(data.voronoi_areas_all),n_area_threshold*n_min_number_of_localizations);
        threshold_idx = 1;
        for k=1:n_area_threshold
            for j=1:n_min_number_of_localizations
                %% calculate cluster for given threshold
                keep_points = data.voronoi_areas_all <= options.area_threshold_arr(k);
                used_points = zeros(1,numel(keep_points));
                idx = find(keep_points);
                idx_clustered = cell(1,floor(numel(idx)/options.min_number_of_localizations_arr(j)));
                for i = 1:numel(idx)
                    if ~used_points(idx(i))
                        neighbors_idx = data.voronoi_neighbors(:,idx(i));
                        neighbors_idx = neighbors_idx(~isnan(neighbors_idx));
                        % get the neighbors with area below the area threshold surrounding the seed-point
                        current_seed = neighbors_idx(keep_points(neighbors_idx));
                        current_seed = current_seed(keep_points(current_seed));
                        if ~isempty(current_seed)
                            size_one = 0;
                            size_two = length(current_seed);
                            % find all connected neighbors above threshold
                            while size_two ~= size_one
                                neighbors_current_seed = data.voronoi_neighbors(:,current_seed);
                                neighbors_current_seed = neighbors_current_seed(~isnan(neighbors_current_seed));
                                size_one = size_two;
                                current_seed = unique(vertcat(neighbors_current_seed,current_seed));
                                current_seed = current_seed(keep_points(current_seed));
                                size_two = length(current_seed);
                            end
                            if length(current_seed) >= options.min_number_of_localizations_arr(j)
                                idx_clustered{i} = current_seed;
                            end
                            used_points(current_seed) = 1;
                        end
                    end   
                end
                % Constructing the segmented clusters
                idx_clustered = idx_clustered(~cellfun('isempty',idx_clustered));
                for i=1:numel(idx_clustered)
                    data.clusters(idx_clustered{i},threshold_idx) = i;
                end
                %%
                clearvars keep_points used_points idx_clustered idx
                threshold_idx = threshold_idx+1;
            end
        end
        clearvars -except filepath data vars_to_load
        data = rmfield(data,vars_to_load);
        save_voronoi_data(filepath, data);
    end
end



