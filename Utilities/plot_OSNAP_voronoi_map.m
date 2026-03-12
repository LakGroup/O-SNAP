% Plot the voronoi maps
function plot_OSNAP_voronoi_map(filepath,data,options)
    [dir_name,name,~] = fileparts(filepath);
    map_dir = replace(dir_name,"_nucleus_data","_nucleus_images");
    if ~exist(map_dir, 'dir')
        mkdir(map_dir)
    end
    % Setup colorbar 
    scale_bar_size_x = 2;
    scale_bar_size_y = 0.5;
    nm_to_um_x = scale_bar_size_x*1000;
    nm_to_um_y = scale_bar_size_y*1000;
    CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
    % Plot reduced voronoi density
    reduced_voronoi_areas_all = data.voronoi_areas_all/mean(data.voronoi_areas_all,'omitnan');
    reduced_log_voronoi_density_plot = log10(1./reduced_voronoi_areas_all);
    idx = reduced_log_voronoi_density_plot >= options.min_log_vor_density;
    reduced_log_voronoi_density_plot = reduced_log_voronoi_density_plot(idx);  
    reduced_log_voronoi_density_plot(reduced_log_voronoi_density_plot>options.max_log_vor_density) = options.max_log_vor_density;
    log_vor_density_ticks = unique(sort(options.min_log_vor_density:0.5:options.max_log_vor_density));
    log_vor_density_ticks = log_vor_density_ticks(all([(log_vor_density_ticks >= options.min_log_vor_density); (log_vor_density_ticks <= options.max_log_vor_density)]));
    f_reduced = figure('visible','off','name','Reduced Log Voronoi Density Map','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
    set(f_reduced, 'InvertHardCopy', 'off'); 
    patch('Vertices',data.vertices,'Faces',data.faces(idx,:),'FaceVertexCData',reduced_log_voronoi_density_plot,'FaceColor','flat','edgecolor','none')    
    set(gca,'colormap',colormap(CT),'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
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
    ylabel(h, 'Reduced $\textrm{Log}_{10}$(Voronoi polygon density) [$\textrm{nm}^{-2}$]','FontSize',14,'interpreter','latex');
    caxis([options.min_log_vor_density options.max_log_vor_density]);
    rectangle('Position',[min(data.x) min(data.y) nm_to_um_x nm_to_um_y],'facecolor','w')
    saveas(f_reduced,fullfile(map_dir,name+"_reduced.png"));
    % Plot Voronoi area
    log_voronoi_area_plot = log10(data.voronoi_areas_all);
    idx = log_voronoi_area_plot >= options.min_log_vor_area;
    log_voronoi_area_plot = log_voronoi_area_plot(idx);  
    log_voronoi_area_plot(log_voronoi_area_plot>options.max_log_vor_area) = options.max_log_vor_area;
    log_vor_area_ticks = unique(sort(options.min_log_vor_area:0.5:options.max_log_vor_area));
    log_vor_area_ticks = log_vor_area_ticks(all([(log_vor_area_ticks >= options.min_log_vor_area); (log_vor_area_ticks <= options.max_log_vor_area)]));
    f = figure('visible','off','name','Log Voronoi Density Map','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
    set(gca,'colormap',flipud(colormap(CT)),'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
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
    ylabel(h, '$\textrm{Log}_{10}$(Voronoi polygon area) [$\textrm{nm}^{2}$]','FontSize',14,'interpreter','latex');
    caxis([options.min_log_vor_area options.max_log_vor_area]);
    rectangle('Position',[min(data.x) min(data.y) nm_to_um_x nm_to_um_y],'facecolor','w')
    saveas(f,fullfile(map_dir,name+".png"));
    close(f_reduced);
    close(f);
end