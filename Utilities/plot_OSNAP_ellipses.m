% -------------------------------------------------------------------------
% plot_OSNAP_ellipses.m
% -------------------------------------------------------------------------
% Renders the ellipses and color codes them based on the density. Also
% provides a plot of the original localization/cluster coordinates above
% the ellipse plot.
% -------------------------------------------------------------------------
% Input:
%   radial_density: [1xn_r] array of the calculated localization density
%                   corresponding to each ring of the radial analysis
%   points: [Nx2] array of 2D coordinate positions to plot
%   x_length: Major axis length of respective nucleus [nm]
%   y_length: Minor axis length of respective nucleus [nm]
% Options:
%   max_density: Maximum density cutoff (across all samples being plotted)
%                to scale colors of ellipses comparably across samples
%   c_data_raw: Information to color-code values. Designed for voronoi areas
%   save_path: The file path to save to
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
function plot_OSNAP_ellipses(radial_density, points, x_length, y_length,options)
arguments
    radial_density double
    points double
    x_length double
    y_length double
    options.max_density double = 3;
    options.c_data_raw double = [];
    options.save_path string = [];
end
    %% options
    min_area = 0.5;
    max_area = 3;
    %% plot rings
    c = parula;
    n_r = size(radial_density,1);
    p = generate_OSNAP_ellipses(x_length, y_length, 1/n_r);
    d_r = [radial_density/options.max_density; 0];
    % d_r = [log10(radial_density/options.max_density); -2];
    % d_r(~isfinite(d_r)) = -2;
    d_rc = round((256-1)*(d_r - min(d_r))/range(d_r)+1);
    if isempty(options.save_path)
        f = figure('visible','on','Position',[616,80,850,900],'color','k');
    else
        f = figure('visible','off','Position',[616,80,850,900],'color','k');
    end
    %% plot rings (bounds)
    h_rings_edge = gca;
    for i=1:n_r
        p_temp = p(i);
        p_temp.Vertices(:,2) = p_temp.Vertices(:,2)+y_length;
        ring_edge = plot(p_temp,"FaceAlpha",0,"EdgeAlpha",0.4,"LineWidth",2,"EdgeColor",'w');%c(d_rc(i),:));
        if i==n_r
            ring_edge.LineStyle="--";
        end
        hold on
        clearvars p_temp
    end
    %% plot voronoi densities
    h_scatter = axes;
    if isempty(options.c_data_raw)
        scatter(points(:,1),points(:,2)+y_length,15,'y','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
    else
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        map = (colormap(CT));
        idx = all([~isnan(options.c_data_raw)  log10(options.c_data_raw) >= min_area],2);
        points = points(idx,:);
        options.c_data_raw = options.c_data_raw(idx);
        % downsample when plotting
        points = points(1:10:end,:);
        options.c_data_raw = options.c_data_raw(1:10:end);
        % color data scales with voronoi area
        c_data = log10(options.c_data_raw);
        c_data(c_data > max_area) = max_area;
        % size data adjusted so that point size scales with voronoi area,
        %    to an extent
        s_data = c_data;
        s_data = 5*(s_data-min(s_data))/(max(s_data)-min(s_data))+0.5;
        % plot
        scatter(points(:,1),points(:,2)+y_length,s_data,c_data,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
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
    b = bar(b_width/2:x_length/(2*n_r):(x_length-b_width)/2,y_length*radial_density/(2*options.max_density),1);
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
    c_max_label = round(options.max_density,ceil(-log10(options.max_density)));
    c_max_tick = (c_max_label*range(cb.Limits)/options.max_density+cb.Limits(1));
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
    if ~isempty(options.save_path)
        set(f,'InvertHardcopy','off');
        print(options.save_path,'-dpng','-r300')
        close(f);
    end
end


