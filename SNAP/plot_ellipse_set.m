function plot_ellipse_set(radial_density, points, voronoi_areas, x_length, y_length, max_density, is_cluster, save_name)
    %% plot rings
    c = jet;
    n_r = size(radial_density,1);
    p = ellipse_radial_density(x_length, y_length, 1/n_r);
    d_r = [radial_density/max_density; 0];
    % d_r = [log10(radial_density/max_density); -2];
    % d_r(~isfinite(d_r)) = -2;
    d_rc = round((256-1)*(d_r - min(d_r))/range(d_r)+1);
    f = figure('visible','off','Position',[616,420,700,500]);
    hax = axes;
    yyaxis(hax,'left');
    for i=1:n_r
        plot(p(i),"FaceColor",c(d_rc(i),:),"FaceAlpha",1,"EdgeAlpha",0);
        hold on
    end
    %% plot voronoi densities
    if is_cluster
        scatter(points(:,1),points(:,2)+y_length,5,'k','filled','MarkerEdgeColor','none');
    else
        colormap jet;
        idx = ~isnan(voronoi_areas);
        points = points(idx,:);
        reduced_voronoi_density = log10(1./voronoi_areas(idx));
        % downsample when plotting
        scatter(points(1:10:end,1),points(1:10:end,2)+y_length,1,reduced_voronoi_density(1:10:end),'filled','MarkerEdgeColor','none');
    end
    hold on
    %%
    xlim([-x_length/2 x_length/2])
    ylim([-y_length/2 3*y_length/2])
    xticks(-x_length/2:x_length/10:x_length/2)
    xticklabels(string(-1:1/5:1))
    hax.YAxis(1).Visible = 'off';
    hax.XGrid = 'on';
    %% plot barplot
    yyaxis(hax,'right');
    b_width = x_length/(2*n_r);
    b=bar(b_width/2:x_length/(2*n_r):(x_length-b_width)/2,radial_density,1);
    b.FaceColor = 'flat';
    b.LineWidth = 1;
    b.CData = c(d_rc(1:n_r),:);
    ylim(max_density*[-1 3])
    hax.YAxis(2).Visible='off';
    colormap jet
    cb = colorbar();
    cb_height = gca().Position(4)/4;
    cb.Position=[0.83,gca().Position(2)+cb_height,0.03,cb_height];
    ylabel(cb, {'Density','(nm^{-2})'},...
        'Rotation',0,'HorizontalAlignment','center')
    c_max_label = round(max_density,ceil(-log10(max_density)));
    c_max_tick = (c_max_label*range(cb.Limits)/max_density+cb.Limits(1));
    cb.Limits = [cb.Limits(1) c_max_tick];
    n_ticks = 4;
    cb.Ticks = cb.Limits(1):range(cb.Limits)/n_ticks:cb.Limits(2);
    cb.TickLabels = compose("%.2e",0:c_max_label/n_ticks:c_max_label);
    cb.Label.Position = [0.5 range(cb.Limits)*0.6 0];
    hax.PlotBoxAspectRatio = [1 2*y_length/x_length 1];
    saveas(f,save_name);
    close(f);
end

