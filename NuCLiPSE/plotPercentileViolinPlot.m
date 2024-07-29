function plotPercentileViolinPlot(saveDir, percentileData, percentile)
    vData = [];
    vCat = [];
    for i=1:length(percentileData)
        iData = percentileData{1, i}.(['prctile_' num2str(percentile)]);
        vData = [vData; iData'];
        vCat = [vCat; repmat(string(percentileData{1, i}.group), length(iData), 1)];
    end
    figure('visible','off');
    violinplot(vData, cellstr(vCat),'EdgeColor',[0 0 0], 'BoxColor',[0 0 0]);
    ylabel('Log Voronoi Density');
    ylim([0 5]);
    title([num2str(percentile) 'th percentile'])
    saveas(gcf,fullfile(saveDir,['VoronoiDensityViolin_' num2str(percentile) '.png']));
end