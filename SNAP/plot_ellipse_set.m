function plot_ellipse_set(radial_density, points, voronoi_areas, x_length, y_length, max_density, is_heterochromatin, save_name)
    %% plot rings
    c = jet;
    n_r = size(radial_density,1);
    p = ellipse_radial_density(x_length, y_length, 1/n_r);
    d_r = [log10(radial_density/max_density); -2];
    d_r(~isfinite(d_r)) = -2;
    d_rc = round(rescale(d_r, 1, 256));
    f = figure('visible','off');
    for i=1:n_r
        plot(p(i),"FaceColor",c(d_rc(i),:),"FaceAlpha",1,"EdgeAlpha",0);
        hold on
    end
    hold on
    %% plot voronoi densities
    if is_heterochromatin
        scatter(points(:,1),points(:,2)+y_length,5,'k','filled','MarkerEdgeColor','none');
    else
        colormap jet;
        idx = ~isnan(voronoi_areas);
        points = points(idx,:);
        reduced_voronoi_density = log10(1./voronoi_areas(idx));
        % downsample when plotting
        scatter(points(1:10:end,1),points(1:10:end,2)+y_length,1,reduced_voronoi_density(1:10:end),'filled','MarkerEdgeColor','none');
    end
    %%
    xlim([-x_length/2 x_length/2])
    ylim([-y_length/2 3*y_length/2])
    xticks(-x_length/2:x_length/10:x_length/2)
    xticklabels(string(-1:1/5:1))
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XGrid = 'on';
    axis equal
    saveas(f,save_name);
    close(f);
end

