% -------------------------------------------------------------------------
% extract_OSNAP_voronoi_cluster_features.m
% -------------------------------------------------------------------------
% Calculates features related to multi-scale Voronoi clustering for a given
% sample (nucleus). Saves information to the MAT analysis file specific to 
% that sample.
%
% Example on how to use it:
%   extract_OSNAP_voronoi_cluster_features("D:\Analysis\ExperimentA\Rep1\sample_KO_1.mat",...
%                   [17.78   31.62   56.23  100.0],...
%                   [5])
% -------------------------------------------------------------------------
% Input:
%   filepath: File path to specific sample (nucleus) to analyze.
%   area_threshold_arr: Area thresholds for multi-scale Voronoi clustering.
%                       Act as upper bound for Voronoi cell areas. Locs
%                       with areas under this value are filtered.
%   min_number_of_localizations: The minimum number of locs in a Voronoi 
%                                cluster. Clusters with fewer locs are
%                                filtered.
% Options:
%   overwrite: Flag on whether to ask user before overwriting data
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function extract_OSNAP_voronoi_cluster_features(filepath,area_threshold_arr,min_number_of_localizations,options)
arguments
    filepath string
    area_threshold_arr double;
    min_number_of_localizations double;
    options.overwrite logical = false;
end

data_vars = {...
    'area_thresholds',...
    'min_number_of_localizations',...
    'voronoi_clusters',...
    'voronoi_cluster_n_locs',...
    'voronoi_cluster_radius',...
    'voronoi_cluster_density',...
    'voronoi_cluster_gyration_radius'};

if has_variables_OSNAP(filepath,data_vars) && ~options.overwrite
    return
end

n_area_threshold = numel(area_threshold_arr);

data = load(filepath);
data.area_thresholds = area_threshold_arr;
data.min_number_of_localizations = min_number_of_localizations;
data.voronoi_clusters = zeros(length(data.voronoi_areas_all),n_area_threshold);
data.voronoi_cluster_n_locs = cell(n_area_threshold,1);
data.voronoi_cluster_radius = cell(n_area_threshold,1);
data.voronoi_cluster_density = cell(n_area_threshold,1);
data.voronoi_cluster_gyration_radius = cell(n_area_threshold,1);

for k=1:n_area_threshold
    %% calculate voronoi cluster for given area threshold
    keep_points = data.voronoi_areas_all <= area_threshold_arr(k);
    used_points = zeros(1,numel(keep_points));
    idx = find(keep_points);
    idx_clustered = cell(1,floor(numel(idx)/min_number_of_localizations));
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
                if length(current_seed) >= min_number_of_localizations
                    idx_clustered{i} = current_seed;
                end
                used_points(current_seed) = 1;
            end
        end   
    end
    % Constructing the segmented clusters
    idx_clustered = idx_clustered(~cellfun('isempty',idx_clustered));
    for i=1:numel(idx_clustered)
        data.voronoi_clusters(idx_clustered{i},k) = i;
    end
    %%
    clearvars keep_points used_points idx_clustered idx

    %% characterize clusters
    n_cluster = max(data.voronoi_clusters(:,k));
    points = [data.x data.y data.voronoi_areas_all data.voronoi_clusters(:,k)];
    points = points(data.voronoi_clusters(:,k) > 0,:); % filter out voronoi_cluster=0 that represents non-clustered data
    points_by_cluster = splitapply(@(x){(x)},points(:,1:3),findgroups(points(:,4))); % Split the clusters into the corresponding groups
    data.voronoi_cluster_n_locs{k} = cellfun('size',points_by_cluster,1); % use cellfun+length, faster from backwards compatibility
    data.voronoi_cluster_radius{k} = zeros(n_cluster,1);
    data.voronoi_cluster_density{k} = zeros(n_cluster,1);
    data.voronoi_cluster_gyration_radius{k} = zeros(n_cluster,1);
    for i=1:n_cluster
        data.voronoi_cluster_radius{k}(i) = sqrt(sum(points_by_cluster{i}(:,3))/pi);
        data.voronoi_cluster_density{k}(i) = data.voronoi_cluster_n_locs{k}(i)./data.voronoi_cluster_radius{k}(i);
        data.voronoi_cluster_gyration_radius{k}(i) = sqrt(mean(...
            (points_by_cluster{i}(:,1)-mean(points_by_cluster{i}(:,1))).^2 ...
            + (points_by_cluster{i}(:,2)-mean(points_by_cluster{i}(:,2))).^2));
    end
end

save_OSNAP_sample(filepath, data);

end




