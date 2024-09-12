function generate_voronoi_data(filepath, options)
arguments
    filepath string
    options.min_log_vor_density double = 0;
    options.max_log_vor_density double = 3;
    options.plot logical = true;
end
if options.plot
    vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density','faces','vertices'};
else
    vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density'};
end
if ~has_variables(filepath,vor_data_vars,"verbose",0)
    data = load_variables(filepath, {'x','y'});
    x = data.x;
    y = data.y;
    % calculate voronoi segmentation
    dt = delaunayTriangulation(x,y);
    [vertices,connections] = voronoiDiagram(dt);
    voronoi_areas_all = nan(numel(connections),1);
    for i = 1:numel(connections)
        tmp = vertices(connections{i},:);
        if ~any(isnan(tmp) | isinf(tmp))
            voronoi_areas_all(i) = polyarea(tmp(:,1),tmp(:,2));
        end
    end
    % clean data
    idx = ~isnan(voronoi_areas_all);
    voronoi_areas = voronoi_areas_all(idx);
    % connections = connections(idx);
    reduced_voronoi_areas = voronoi_areas / mean(voronoi_areas);
    reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
    % find neighbors
    connectivity_list = dt.ConnectivityList;
    attached_triangles = vertexAttachments(dt);
    neighbors = cell(numel(attached_triangles),1);
    for i = 1:numel(attached_triangles)
        neighbors{i} = connectivity_list(attached_triangles{i},:);
        neighbors{i} = unique(neighbors{i});
        neighbors{i}(neighbors{i}==i) = [];
    end
    % return relevant info
    data.voronoi_areas_all = voronoi_areas_all;
    data.voronoi_areas = voronoi_areas;
    data.voronoi_neighbors = uneven_cell2mat(neighbors);
    data.reduced_log_voronoi_density = reduced_log_voronoi_density;
    data.faces = uneven_cell2mat(connections)';
    data.vertices = vertices;
else
    data = load_variables(filepath, {'x','y','voronoi_areas_all','faces','vertices'});
end
    %% plot voronoi segmentation with reduced densities
    if options.plot
        reduced_voronoi_areas_all = data.voronoi_areas_all/mean(data.voronoi_areas_all,'omitnan');
        reduced_log_voronoi_density_plot = log10(1./reduced_voronoi_areas_all);
        idx = reduced_log_voronoi_density_plot >= options.min_log_vor_density;
        reduced_log_voronoi_density_plot = reduced_log_voronoi_density_plot(idx);  
        reduced_log_voronoi_density_plot(reduced_log_voronoi_density_plot>options.max_log_vor_density) = options.max_log_vor_density;
        % colorbar setup
        scale_bar_size_x = 2;
        scale_bar_size_y = 0.5;
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        map = (colormap(CT));
        f = figure('visible','off','name','Log Voronoi density Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
        patch('Vertices',data.vertices,'Faces',data.faces(idx,:),'FaceVertexCData',reduced_log_voronoi_density_plot,'FaceColor','flat','edgecolor','none')    
        axis equal
        xlim([min(data.x) max(data.x)])
        ylim([min(data.y) max(data.y)])    
        h = colorbar();    
        pbaspect([1 1 1])
        h.Position = [0.9 0.1 0.02 0.8];
        h.TickLabelInterpreter = 'latex';
        h.FontSize = 14;
        h.Color = 'w';
        log_vor_density_ticks = unique(sort([0 0.1 0.25, 0.41, 0.7, 1, 1.7, 2, options.min_log_vor_density, options.max_log_vor_density]));
        log_vor_density_ticks = log_vor_density_ticks(all([(log_vor_density_ticks >= options.min_log_vor_density); (log_vor_density_ticks <= options.max_log_vor_density)]));
        h.Ticks = log_vor_density_ticks;
        h.TickLabels = arrayfun(@(elem) sprintf('%.2f',elem), log_vor_density_ticks,  'uni',  0);
        ylabel(h, 'Reduced Log Voronoi polygon density $(\textrm{nm}^{-2})$','FontSize',14,'interpreter','latex');
        caxis([options.min_log_vor_density options.max_log_vor_density]);
        pixels_um_x = scale_bar_size_x*1000;
        pixels_um_y = scale_bar_size_y*1000;
        rectangle('Position',[min(data.x) min(data.y) pixels_um_x pixels_um_y],'facecolor','w')
        set(gca,'colormap',map,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
        set(f, 'InvertHardCopy', 'off'); 
        [dir_name,name,~] = fileparts(filepath);
        map_dir = replace(dir_name,"voronoi_data","nuclei_images");
        if ~exist(map_dir, 'dir')
            mkdir(map_dir)
        end
        saveas(f,fullfile(map_dir,name+"_reduced.png"));
        close();
    end
    clearvars -except data filepath
    save_voronoi_data(filepath, data);
end

% change structure of neighbors from cell array to matrix
function arr_new = uneven_cell2mat(arr)
    arr_new = nan(max(cellfun('length',arr)),length(arr));
    for i=1:length(arr)
        value = arr{i};
        arr_new(1:numel(value),i) = value;
    end
end
