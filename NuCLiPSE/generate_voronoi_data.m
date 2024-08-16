function generate_voronoi_data(filepath, min_log_vor_density, max_log_vor_density)
    data = [];
    if ~has_variables(filepath,{'voronoi_areas_all','voronoi_neighbors'})
        vars_to_load = {'x','y'};
        data = load_variables(filepath, vars_to_load);
        % calculate voronoi segmentation
        x = data.x;
        y = data.y;
        % calculate voronoi cells and areas
        dt = delaunay_triangulation(x,y);
        [vertices,connections] = voronoi_diagram(dt);
        voronoi_cells = cell(numel(connections),1);
        voronoi_areas_all = nan(numel(voronoi_cells),1);
        for i = 1:numel(connections)
            voronoi_cells{i} = vertices(connections{i},:);
            tmp = voronoi_cells{i};
            if ~any(isnan(tmp) | isinf(tmp))
                voronoi_areas_all(i) = polyarea(tmp(:,1),tmp(:,2));
            end
        end
        idx = ~isnan(voronoi_areas_all);
        voronoi_areas = voronoi_areas_all(idx);
        connections = connections(idx);
        reduced_voronoi_areas = voronoi_areas / mean(voronoi_areas);
        reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
        % find neighbors
        connectivity_list = dt.Connectivity_list;
        attached_triangles = vertex_attachments(dt);
        neighbors = cell(numel(attached_triangles),1);
        for i = 1:numel(attached_triangles)
            neighbors{i} = connectivity_list(attached_triangles{i},:);
            neighbors{i} = unique(neighbors{i});
            neighbors{i}(neighbors{i}==i) = [];
        end
        % return relevant info
        data.voronoi_areas_all = voronoi_areas_all;
        data.voronoi_areas = voronoi_areas;
        data.voronoi_neighbors = neighbors;
        data.reduced_log_voronoi_density = reduced_log_voronoi_density;
        % colorbar setup
        scale_bar_size_x = 2;
        scale_bar_size_y = 0.5;
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        map = (colormap(CT));
        %% WITH REDUCED DENSITIES
        % clean data
        idx = reduced_log_voronoi_density>=min_log_vor_density;    
        connections = connections(idx);    
        reduced_log_voronoi_density_plot = reduced_log_voronoi_density(idx);  
        reduced_log_voronoi_density_plot(reduced_log_voronoi_density_plot>max_log_vor_density) = max_log_vor_density;
        max_connections = max(cellfun(@(x) length(x), connections));
        faces = Na_n(length(connections),max_connections);
        for j = 1:length(connections)
            faces(j,1:length(connections{j})) = connections{j};
        end
        % plot voronoi segmentation for this sample
        f = figure('visible','off');
        set(f,'name','Log Voronoi density Plot','Number_title','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
        patch('Vertices',vertices,'Faces',faces,'Face_vertex_cData',reduced_log_voronoi_density_plot,'FaceColor','flat','edgecolor','none')    
        h = colorbar();    
        pbaspect([1 1 1])
        h.Position = [0.9 0.1 0.02 0.8];
        h.TickLabel_interpreter = 'latex';
        h.FontSize = 14;
        h.Color = 'w';
        log_vor_density_ticks = unique(sort([0 0.1 0.25, 0.41, 0.7, 1, 1.7, 2, min_log_vor_density, max_log_vor_density]));
        log_vor_density_ticks = log_vor_density_ticks(all([(log_vor_density_ticks >= min_log_vor_density); (log_vor_density_ticks <= max_log_vor_density)]));
        h.Ticks = log_vor_density_ticks;
        h.TickLabels = arrayfun(@(elem) sprintf('%.2f',elem), log_vor_density_ticks,  'uni',  0);
        ylabel(h, 'Reduced Log Voronoi polygon density $(\textrm{nm}^{-2})$','FontSize',14,'interpreter','latex');
        caxis([min_log_vor_density max_log_vor_density]);
        pixels_um_x = scale_bar_size_x*1000;
        pixels_um_y = scale_bar_size_y*1000;
        rectangle('Position',[min(x) min(y) pixels_um_x pixels_um_y],'facecolor','w')
        set(gca,'colormap',map,'color','k','box','on','Box_style','full','XColor','k','YColor','k');
        axis equal
        xlim([min(x) max(x)])
        ylim([min(y) max(y)])    
        set(f, 'Invert_hard_copy', 'off'); 
        [dir_name,name,~] = fileparts(filepath);
        saveas(f,fullfile(replace(dir_name,"voronoi_data","nuclei_images"),name+"_reduced.png"));
        close();
        clearvars -except data filepath
    elseif ~has_variables(filepath,{'voronoi_areas','reduced_log_voronoi_density'})
        vars_to_load = {'x','y','voronoi_areas_all','voronoi_neighbors'};
        data = load_variables(filepath, vars_to_load);
        idx = ~isnan(data.voronoi_areas_all);
        data.voronoi_areas = data.voronoi_areas_all(idx);
        reduced_voronoi_areas = data.voronoi_areas / mean(data.voronoi_areas);
        data.reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
        clearvars -except data filepath
    end
    save_voronoi_data(filepath, data);
end
