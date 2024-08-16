function plot_percentile_violin_plot(save_dir, percentile_data, percentile)
    v_data = [];
    v_cat = [];
    for i=1:length(percentile_data)
        i_data = percentile_data{1, i}.(['prctile_' num2str(percentile)]);
        v_data = [v_data; i_data'];
        v_cat = [v_cat; repmat(string(percentile_data{1, i}.group), length(i_data), 1)];
    end
    figure('visible','off');
    violin_plot(v_data, cellstr(v_cat),'EdgeColor',[0 0 0], 'BoxColor',[0 0 0]);
    ylabel('Log Voronoi density');
    ylim([0 5]);
    title([num2str(percentile) 'th percentile'])
    saveas(gcf,fullfile(save_dir,['voronoi_density_violin_' num2str(percentile) '.png']));
end
