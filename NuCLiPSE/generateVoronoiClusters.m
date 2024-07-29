function generateVoronoiClusters(filepath,area_threshold_arr,min_number_of_localizations_arr)
    if ~hasVars(filepath,{'area_thresholds','min_number_of_localizations','clusters'})
        varsToLoad = {'voronoi_areas_all','voronoi_neighbors'};
        data = loadVars(filepath,varsToLoad);
        n_area_threshold = length(area_threshold_arr);
        n_min_number_of_localizations = length(min_number_of_localizations_arr);
        data.area_thresholds = area_threshold_arr;
        data.min_number_of_localizations = min_number_of_localizations_arr;
        data.clusters = zeros(length(data.voronoi_areas_all),n_area_threshold*n_min_number_of_localizations);
        threshold_idx = 1;
        for k=1:n_area_threshold
            for j=1:n_min_number_of_localizations
                %% calculate cluster for given threshold
                keep_points = data.voronoi_areas_all <= area_threshold_arr(k);
                counter = 0;
                number_of_points = size(data.voronoi_areas_all,1);
                used_points = zeros(number_of_points,1);
                idx_clustered = [];
                for i = 1:number_of_points
                    if keep_points(i) && ~used_points(i)
                        % get the neighbors with area below the area threshold surrounding the seed-point
                        seed_neighbors = data.voronoi_neighbors{i}(keep_points(data.voronoi_neighbors{i}));
                        if ~isempty(seed_neighbors)
                            size_one = 0;
                            size_two = length(seed_neighbors);
                            % find all connected neighbors above min_loc threshold
                            while size_two ~= size_one
                                size_one = length(seed_neighbors);
                                idx_all = unique(cell2mat(data.voronoi_neighbors(seed_neighbors)));
                                if ~any(intersect(idx_all,seed_neighbors))
                                    seed_neighbors = sort([idx_all;seed_neighbors]);
                                else
                                    seed_neighbors = idx_all;
                                end
                                seed_neighbors = seed_neighbors(keep_points(seed_neighbors));
                                size_two = length(seed_neighbors);
                            end
                        else
                            seed_neighbors = i;
                        end
                        used_points(seed_neighbors) = 1;
                        if length(seed_neighbors) >= min_number_of_localizations_arr(j)
                            counter = counter+1;
                            idx_clustered{counter} = seed_neighbors;
                        end
                        clearvars seed_neighbors idx_all
                    end    
                end
                if ~isempty(idx_clustered)
                    idx_not_clustered = vertcat(idx_clustered{:});
                    idx_not_clustered = setxor(1:length(data.voronoi_neighbors),idx_not_clustered);
                    idx = 1:length(data.voronoi_neighbors);
                    idx(idx_not_clustered) = -1;
                    for i = 1:length(idx_clustered)
                        idx(idx_clustered{i}) = i;
                    end
                    [~,~,data.clusters(:,threshold_idx)] = unique(idx);
                    data.clusters(:,threshold_idx) = data.clusters(:,threshold_idx)-1;
                else
                    data.clusters(:,threshold_idx) = zeros(number_of_points,1);
                end
                %%
                clearvars keep_points counter used_points idx_clustered idx
                threshold_idx = threshold_idx+1;
            end
        end
        clearvars -except filepath data varsToLoad
        data = rmfield(data,varsToLoad);
        saveVoronoiAnalysisData(filepath, data);
    end
end

