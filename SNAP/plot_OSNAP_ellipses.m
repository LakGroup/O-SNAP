function plot_OSNAP_ellipses(radial_density, points, voronoi_areas, x_length, y_length, max_density, is_cluster, save_name)
    %% options
    min_area = 0.5;
    max_area = 3;
    %% plot rings
    c = parula;
    n_r = size(radial_density,1);
    p = generate_OSNAP_ellipses(x_length, y_length, 1/n_r);
    d_r = [radial_density/max_density; 0];
    % d_r = [log10(radial_density/max_density); -2];
    % d_r(~isfinite(d_r)) = -2;
    d_rc = round((256-1)*(d_r - min(d_r))/range(d_r)+1);
    f = figure('visible','on','Position',[616,80,850,900],'color','k');
    %% plot rings (bounds)
    h_rings_edge = gca;
    for i=1:n_r
        p_temp = p(i);
        p_temp.Vertices(:,2) = p_temp.Vertices(:,2)+y_length;
        ring_edge = plot(p_temp,"FaceAlpha",0,"EdgeAlpha",0.4,"EdgeColor",c(d_rc(i),:));
        if i==n_r
            ring_edge.LineStyle="--";
        end
        hold on
        clearvars p_temp
    end
    %% plot voronoi densities
    h_scatter = axes;
    if is_cluster
        scatter(points(:,1),points(:,2)+y_length,15,'y','filled','MarkerEdgeColor','k');
    else
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        map = (colormap(CT));
        idx = all([~isnan(voronoi_areas)  log10(voronoi_areas) >= min_area],2);
        points = points(idx,:);
        voronoi_areas = voronoi_areas(idx);
        % downsample when plotting
        points = points(1:10:end,:);
        voronoi_areas = voronoi_areas(1:10:end);
        % color data scales with voronoi area
        c_data = log10(voronoi_areas);
        c_data(c_data > max_area) = max_area;
        % size data adjusted so that point size scales with voronoi area,
        %    to an extent
        s_data = c_data;
        s_data = 5*(s_data-min(s_data))/(max(s_data)-min(s_data))+0.5;
        % plot
        scatter(points(:,1),points(:,2)+y_length,s_data,c_data,'filled','MarkerEdgeColor','none');
        clim(h_scatter,[min_area max_area]);
        colormap(h_scatter,flipud(map));
    end
    %% plot rings (filled)
    h_rings_filled = axes;
    for i=1:n_r
        plot(p(i),"FaceColor",c(d_rc(i),:),"FaceAlpha",1,"EdgeAlpha",0);
        hold on
    end
    colormap(h_rings_filled,"parula");
    %% plot barplot
    h_bar = axes;
    b_width = x_length/(2*n_r);
    b = bar(b_width/2:x_length/(2*n_r):(x_length-b_width)/2,y_length*radial_density/(2*max_density),1);
    b.FaceColor = 'flat';
    b.LineWidth = 1;
    b.CData = c(d_rc(1:n_r),:);
    %% save axis information
    axpos = [0.06,0.17,0.55,0.8];
    %% colorbar for barplot
    cb = colorbar('Color','w','Location','eastoutside','FontSize',24);
    ylabel(cb, {'Density [nm^{-2}]'},'HorizontalAlignment','center')
    xlabel("Percent of major/minor axis")
    xticklabels(["","","","","","0%","20%","40%","60%","80%","100%"])
    set(gca,'FontSize',24);
    c_max_label = round(max_density,ceil(-log10(max_density)));
    c_max_tick = (c_max_label*range(cb.Limits)/max_density+cb.Limits(1));
    cb.Limits = [cb.Limits(1) c_max_tick];
    n_ticks = 4;
    cb.Ticks = cb.Limits(1):range(cb.Limits)/n_ticks:cb.Limits(2);
    cb.TickLabels = compose("%.2e",0:c_max_label/n_ticks:c_max_label);
    colormap(h_bar,"parula");
    %% adjust axes geometry
    xticks(-x_length/2:x_length/10:x_length/2)
    set(h_rings_edge,'position',axpos,'PlotBoxAspectRatio', [1 2*y_length/x_length 1],'color','k','xlim',[-x_length/2 x_length/2],'ylim',[-y_length/2 3*y_length/2],'xtick',[],'ytick',[]);
    set(h_scatter,'position',axpos,'PlotBoxAspectRatio',[1 2*y_length/x_length 1],'color','none','xlim',[-x_length/2 x_length/2],'ylim',[-y_length/2 3*y_length/2],'xtick',[],'ytick',[]);
    set(h_rings_filled,'position',axpos,'PlotBoxAspectRatio', [1 2*y_length/x_length 1],'color','none','xlim',[-x_length/2 x_length/2],'ylim',[-y_length/2 3*y_length/2],'xtick',[],'ytick',[]);
    set(h_bar,'position',axpos,'PlotBoxAspectRatio', [1 2*y_length/x_length 1],'color','none','xlim',[-x_length/2 x_length/2],'ylim',[-y_length/2 3*y_length/2],'ytick',[]);
    set(h_bar,'XGrid','on','GridColor','w','XColor','w');
    %% save
    if ~isempty(save_name)
        set(f,'InvertHardcopy','off');
        print(save_name,'-dpng','-r300')
        close(f);
    end
end

