function plotRadialDensity(radialDensity, points, voronoi_areas, x_length, y_length, maxDensity, pointSize, saveName)
    %% plot rings
    c = jet;
    n_r = size(radialDensity,1);
    p = ellipseRadialDensity(x_length, y_length, 1/n_r);
    d_r = [log10(radialDensity/maxDensity); -2];
    d_r(~isfinite(d_r)) = -2;
    d_rc = round(rescale(d_r, 1, 256));
    f = figure('visible','off');
    for i=1:n_r
        plot(p(i),"FaceColor",c(d_rc(i),:),"FaceAlpha",1,"EdgeAlpha",0);
        hold on
    end
    hold on
    %% plot voronoi densities
    colormap jet;
    % remove NaNs
    idx = ~isnan(voronoi_areas);
    points = points(idx,:);
    reduced_voronoi_density = log10(1./voronoi_areas(idx));
    % downsample when plotting
    scatter(points(1:10:end,1),points(1:10:end,2)+y_length,pointSize,reduced_voronoi_density(1:10:end),'filled','MarkerEdgeColor','none');
    %%
    xlim([-x_length/2 x_length/2])
    ylim([-y_length/2 3*y_length/2])
    xticks(-x_length/2:x_length/10:x_length/2)
    xticklabels(string(-1:1/5:1))
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XGrid = 'on';
    axis equal
    saveas(f,saveName);
    close(f);
end