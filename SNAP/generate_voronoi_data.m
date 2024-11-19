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
    reduced_voronoi_areas = voronoi_areas ./ mean(voronoi_areas);
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
        [dir_name,name,~] = fileparts(filepath);
        map_dir = replace(dir_name,"voronoi_data","nuclei_images");
        if ~exist(map_dir, 'dir')
            mkdir(map_dir)
        end
        %% setup
        % colorbar 
        scale_bar_size_x = 2;
        scale_bar_size_y = 0.5;
        nm_to_um_x = scale_bar_size_x*1000;
        nm_to_um_y = scale_bar_size_y*1000;
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        map = (colormap(CT));
        %% reduced voronoi density
        reduced_voronoi_areas_all = data.voronoi_areas_all/mean(data.voronoi_areas_all,'omitnan');
        reduced_log_voronoi_density_plot = log10(1./reduced_voronoi_areas_all);
        idx = reduced_log_voronoi_density_plot >= options.min_log_vor_density;
        reduced_log_voronoi_density_plot = reduced_log_voronoi_density_plot(idx);  
        reduced_log_voronoi_density_plot(reduced_log_voronoi_density_plot>options.max_log_vor_density) = options.max_log_vor_density;
        log_vor_density_ticks = unique(sort([0 0.1 0.25, 0.41, 0.7, 1, 1.7, 2, 3, options.min_log_vor_density, options.max_log_vor_density]));
        log_vor_density_ticks = log_vor_density_ticks(all([(log_vor_density_ticks >= options.min_log_vor_density); (log_vor_density_ticks <= options.max_log_vor_density)]));
        f_reduced = figure('visible','off','name','Reduced Log Voronoi Density Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
        set(gca,'colormap',map,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
        set(f_reduced, 'InvertHardCopy', 'off'); 
        patch('Vertices',data.vertices,'Faces',data.faces(idx,:),'FaceVertexCData',reduced_log_voronoi_density_plot,'FaceColor','flat','edgecolor','none')    
        axis equal
        xlim([min(data.x) max(data.x)])
        ylim([min(data.y) max(data.y)])    
        h = colorbar('Position',[0.9 0.1 0.02 0.8],...
            'TickLabelInterpreter','latex',...
            'FontSize',14,...
            'Color','w',...
            'Ticks',log_vor_density_ticks,...
            'TickLabels',arrayfun(@(elem) sprintf('%.2f',elem), log_vor_density_ticks,  'uni',  0));    
        pbaspect([1 1 1])
        ylabel(h, 'Reduced $\textrm{Log}_{10}$ Voronoi polygon density $(\textrm{nm}^{-2})$','FontSize',14,'interpreter','latex');
        caxis([options.min_log_vor_density options.max_log_vor_density]);
        rectangle('Position',[min(data.x) min(data.y) nm_to_um_x nm_to_um_y],'facecolor','w')
        saveas(f_reduced,fullfile(map_dir,name+"_reduced.png"));
        close();
        %% voronoi area
        log_voronoi_area_plot = log10(data.voronoi_areas_all);
        idx = log_voronoi_area_plot >= options.min_log_vor_density;
        log_voronoi_area_plot = log_voronoi_area_plot(idx);  
        log_voronoi_area_plot(log_voronoi_area_plot>options.max_log_vor_density) = options.max_log_vor_density;
        log_vor_area_ticks = unique(sort([0.5, 0.7, 1, 1.7, 2, 3]));
        log_vor_area_ticks = log_vor_area_ticks(all([(log_vor_area_ticks >=0.5); (log_vor_area_ticks <= 3)]));
        f = figure('visible','off','name','Log Voronoi Density Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
        set(gca,'colormap',flipud(map),'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
        set(f, 'InvertHardCopy', 'off'); 
        patch('Vertices',data.vertices,'Faces',data.faces(idx,:),'FaceVertexCData',log_voronoi_area_plot,'FaceColor','flat','edgecolor','none')    
        axis equal
        xlim([min(data.x) max(data.x)])
        ylim([min(data.y) max(data.y)])    
        h = colorbar('Position',[0.9 0.1 0.02 0.8],...
            'TickLabelInterpreter','latex',...
            'FontSize',14,...
            'Color','w',...
            'Ticks',log_vor_area_ticks,...
            'TickLabels',arrayfun(@(elem) sprintf('%.2f',elem), log_vor_area_ticks,  'uni',  0));    
        pbaspect([1 1 1])
        ylabel(h, '$\textrm{Log}_{10}$ Voronoi polygon area $(\textrm{nm}^{2})$','FontSize',14,'interpreter','latex');
        caxis([0.5 3]);
        rectangle('Position',[min(data.x) min(data.y) nm_to_um_x nm_to_um_y],'facecolor','w')
        saveas(f,fullfile(map_dir,name+".png"));
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
