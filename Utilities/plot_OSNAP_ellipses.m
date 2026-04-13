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
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
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
    %% Options
    min_area = 0.5;
    max_area = 3;
    %% Plot rings
    n_r = size(radial_density,1);
    p = generate_OSNAP_ellipses(x_length, y_length, 1/n_r);
    d_r = min([radial_density/options.max_density; 1],1);
    d_rc = floor(255*d_r+1);
    if isempty(options.save_path)
        f = figure('visible','on','Position',[616,80,850,900],'color','k');
    else
        f = figure('visible','off','Position',[616,80,850,900],'color','k');
    end
    tiledlayout(2,1,"tilespacing","compact");
    %% Plot rings (bounds)
    ax = nexttile();
    for i=1:n_r
        p_temp = p(i);
        p_temp.Vertices(:,2) = p_temp.Vertices(:,2);
        ring_edge = plot(p_temp,"FaceAlpha",0,"EdgeAlpha",0.4,"LineWidth",2,"EdgeColor",'w');%c(d_rc(i),:));
        if i==n_r
            ring_edge.LineStyle="--";
        end
        hold on
        clearvars p_temp
    end
    %% Plot voronoi densities
    if isempty(options.c_data_raw)
        s = scatter(points(:,1),points(:,2),25,'w','filled','MarkerEdgeColor','k');
        legend(s,"domain centroid","location","eastoutside","fontsize",16,"color","k","TextColor","w")
    else
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        colormap(ax, flipud(CT));
        idx = all([~isnan(options.c_data_raw)  log10(options.c_data_raw) >= min_area],2);
        for i=1:2
            loc_diff = max(points(:,i))+min(points(:,i));
            points(:,i) = points(:,i) - loc_diff/2;
        end
        points = points(idx,:);
        options.c_data_raw = options.c_data_raw(idx);
        % Downsample when plotting
        points = points(1:10:end,:);
        options.c_data_raw = options.c_data_raw(1:10:end);
        % Color data scales with voronoi area
        c_data = log10(options.c_data_raw);
        c_data(c_data > max_area) = max_area;
        % Size data adjusted so that point size scales with voronoi area,
        % to an extent
        s_data = c_data;
        s_data = 5*(s_data-min(s_data))/(max(s_data)-min(s_data))+0.5;
        % Plot
        scatter(points(:,1),points(:,2),s_data,c_data,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
        clim([min_area max_area]);
        cb_1 = colorbar('Color','w','Location','eastoutside','FontSize',20,...
            'Ticks',min_area:0.5:max_area,...
            'TickLabels',arrayfun(@(elem) sprintf('10^{%g}',elem), min_area:0.5:max_area,'uni',0),...
            'tickdirection','out');
        ylabel(cb_1, {'Voronoi polygon area [nm^{2}]'},'HorizontalAlignment','center')
    end
    set(ax,'color','k','xlim',[-x_length/2 x_length/2],'ylim',[-y_length/2 y_length/2],'xticklabels',[],'ytick',[],...
        'xtick',linspace(-x_length/2,x_length/2,11),'XGrid','on','GridColor','w','XColor','none');
    daspect([1 1 1])
    disp("")
   
    %% Plot rings (filled)
    ax=nexttile();
    c = interp1(linspace(0,1,7),[0.0 0.0 0.0; 0.1412 0.1804 0.5059; 0.1804 0.4157 0.6824; 0.2275 0.6196 0.749; 0.4627 0.7765 0.7412; 0.749 0.898 0.7373; 1.0 1.0 0.8],linspace(0,1,256),'linear','extrap');
    for i=1:n_r
        plot(p(i),"FaceColor",c(d_rc(i),:),"FaceAlpha",1,"EdgeAlpha",0);
        hold on
    end
    colormap(ax,c);
    %% Plot barplot
    b_width = x_length/(2*n_r);
    b = bar(b_width/2:x_length/(2*n_r):(x_length-b_width)/2,(y_length/2)*(radial_density/options.max_density),1);
    b.FaceColor = 'flat';
    b.LineWidth = 1;
    b.CData = c(d_rc(1:n_r),:);
    %% Colorbar for barplot
    cb_2 = colorbar('Color','w','Location','eastoutside','FontSize',20,'tickdirection','out');
    xlabel("Ring Index",'Color','w')
    set(ax,'xtick',linspace(x_length/20,x_length/2,10),'xticklabels',compose("%g",1:10));
    c_max_label = round(options.max_density,ceil(-log10(options.max_density)));
    c_max_tick = (c_max_label*range(cb_2.Limits)/options.max_density+cb_2.Limits(1));
    cb_2.Limits = [cb_2.Limits(1) c_max_tick];
    n_ticks = 4;
    cb_2.Ticks = cb_2.Limits(1):range(cb_2.Limits)/n_ticks:cb_2.Limits(2);
    cb_2.TickLabels = compose("%.1e",0:c_max_label/n_ticks:c_max_label);
    daspect([1 1 1]),
    set(ax,'FontSize',14,'color','k','xlim',[-x_length/2 x_length/2],'ylim',[-y_length/2 y_length/2],'ytick',[],...
        'XGrid','on','GridColor','w','XColor','w','XTickLabelRotation',30); 
    ax.XLabel.FontSize = 22;
    if isempty(options.c_data_raw)
        ylabel(cb_2, {'Radial Domain Density [nm^{-2}]'},'HorizontalAlignment','center')
        cb_1.Position(1) = cb_2.Position(1);
        cb_1.Label.Position(1) = cb_2.Label.Position(1);
    else
        ylabel(cb_2, {'Radial Localization Density [nm^{-2}]'},'HorizontalAlignment','center')
    end
    %% Save
    if ~isempty(options.save_path)
        set(f,'InvertHardcopy','off');
        print(options.save_path,'-dpng','-r300')
        close(f);
    end
end


