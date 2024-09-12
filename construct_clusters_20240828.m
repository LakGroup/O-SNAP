function [data_clustered,data_not_clustered,area_threshold] = construct_clusters(data,area_threshold,min_number_of_localizations,type,waitbar_counter)
f = waitbar(0,['Please Wait...',num2str(waitbar_counter(1)),'/',num2str(waitbar_counter(2))]);
voronoi_areas = data.vor.voronoi_areas;
neighbors = data.vor.neighbors;
if type == 2    
    area_threshold = prctile(voronoi_areas,area_threshold);
end
keep_points = voronoi_areas <= area_threshold;
counter = 0;
number_of_points = size(voronoi_areas,1);
used_points = zeros(number_of_points,1);

waitbar(0.1,f,['Please Wait...',num2str(waitbar_counter(1)),'/',num2str(waitbar_counter(2))]);
for i = 1:number_of_points
    if keep_points(i) && ~used_points(i)
        % get the neighbors with area below the area threshold surrounding the seed-point
        seed_neighbors = neighbors{i}(keep_points(neighbors{i}));
        if ~isempty(seed_neighbors)
            size_one = 0;
            size_two = length(seed_neighbors);
            % find all connected neighbors above threshold
            while size_two ~= size_one
                size_one = length(seed_neighbors);
                idx_all = unique(cell2mat(neighbors(seed_neighbors)));
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
        if length(seed_neighbors) >= min_number_of_localizations
            counter = counter+1;
            idx_clustered{counter} = seed_neighbors;
        end
    end    
end

waitbar(0.3,f,['Please Wait...',num2str(waitbar_counter(1)),'/',num2str(waitbar_counter(2))]);
idx_not_clustered = vertcat(idx_clustered{:});
idx_not_clustered = setxor(1:length(data.vor.neighbors),idx_not_clustered);
idx = 1:length(data.vor.neighbors);
idx(idx_not_clustered) = -1;
for i = 1:length(idx_clustered)
    idx(idx_clustered{i}) = i;
end

waitbar(0.5,f,['Please Wait...',num2str(waitbar_counter(1)),'/',num2str(waitbar_counter(2))]);
data_vor(:,1) = data.vor.points(:,1);
data_vor(:,2) = data.vor.points(:,2);
data_vor(:,3) = data.vor.voronoi_areas;

I = idx == -1;

waitbar(0.7,f,['Voronoi Segmentation Not-clustered Data...',num2str(waitbar_counter(1)),'/',num2str(waitbar_counter(2))]);
data_not_clustered.x_data = data_vor(I,1);
data_not_clustered.y_data = data_vor(I,2);
data_not_clustered.area = zeros(length(data_not_clustered.x_data),1);
data_not_clustered.name = [data.name,'_not_clustered_threshold_',num2str(area_threshold),'_min_points_',num2str(min_number_of_localizations)];
data_not_clustered.type = 'loc_list';

waitbar(0.9,f,['Voronoi Segmentation Clustered Data...',num2str(waitbar_counter(1)),'/',num2str(waitbar_counter(2))]);
idx(I) = [];
data_vor(I,:) = [];
% data_voronoi_cells = data.vor.voronoi_cells;
% data_voronoi_cells(I) = [];

%clusters = loc_list_find_clusters(data_vor,idx);
clusters = loc_list_find_clusters_inside(data_vor,idx);
clusters = vertcat(clusters{:});
%poly_shape = find_polyshape(clusters_vor_cells);

data_clustered.x_data = clusters(:,1);
data_clustered.y_data = clusters(:,2);
data_clustered.area = clusters(:,3);
data_clustered.name = [data.name,'_clusters_threshold_',num2str(area_threshold),'_min_points_',num2str(min_number_of_localizations)];
data_clustered.type = 'loc_list';
close(f)
end

function clusters = loc_list_find_clusters_inside(data,idx)
idx_unique = unique(idx);

for i = 1:length(idx_unique)
    I = idx == idx_unique(i);  
    %clusters_vor_cells{i} = data_vor_cell(I);
    clusters{i}(:,1:2) = data(I,1:2);
    clusters{i}(:,3) = sum(data(I,3));
    clear I
end
end

function poly_shape = find_polyshape(data)
f = waitbar(0,'claculating polyshapes');
for i = 1:length(data)
    for j = 1:length(data{i})
        polygons = polyshape(data{i}{j}(:,1),data{i}{j}(:,2));
        if j==1
            poly_union = polygons;
        elseif j>1
            poly_union = union(poly_union,polygons);
        end        
    end
    poly_shape{i} = poly_union;
    clear poly_union
    waitbar(i/length(data),f,'claculating polyshapes');
end
close(f)
end