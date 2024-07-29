function plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData, suffix, nbins)
%% set up for cdf plots
n_groups = length(groups);
n_samps = length(voronoiData);
[~, groupIdx] = sortSampleToGroups(groups,voronoiData);
name2color = zeros(n_samps,3);
name2style = cell(n_samps);
linewidth = 2;
group_binned_data = cell(n_groups,1);
% determine min/max for x-axis
min_val_cdf = min(cellfun(@(x) prctile(x.reduced_log_voronoi_density, 1), voronoiData, 'uni',1));
max_val_cdf = max(cellfun(@(x) prctile(x.reduced_log_voronoi_density, 99.5), voronoiData, 'uni',1));
for g=1:n_groups
    nSampInGroup = length(groupIdx{g});
    if nSampInGroup > 0
        group_binned_data{g} = zeros(nSampInGroup,nbins);
        for i=1:nSampInGroup
            % organize styles
            idx = groupIdx{g}(i);
            name2color(idx,:) = repmat(group_colors(g,:),length(idx),1);
            name2style{idx} = group_styles{g};
            % fill a 2D array containing the binned log voronoi density of each sample in each row
            group_binned_data{g}(i,:) = histcounts(voronoiData{idx}.reduced_log_voronoi_density,'Normalization','cumcount','BinEdges',linspace(min_val_cdf,max_val_cdf,nbins+1));
        end
    end
end
% plot average cdfs
f = figure('visible','off');
f.Position = [0 0 1067 600];
x = linspace(min_val_cdf,max_val_cdf,nbins);
L = gobjects(1,n_groups);
grid on;
for g=1:n_groups
    nSampInGroup = length(groupIdx{g});
    if nSampInGroup < 1
        continue;
    end
    norm_data = bsxfun(@rdivide, group_binned_data{g}, max(group_binned_data{g}, [], 2)); % normalized based on row-max
    y = mean(norm_data,1);
    s = std(norm_data,1);
    patch([x fliplr(x)], [(y-s) fliplr(y+s)], group_colors(g,:),...
        'FaceAlpha',0.05,'LineWidth',0.2,'LineStyle',':',...
        'EdgeColor',group_colors(g,:),'EdgeAlpha','0.2');
    hold on;
    norm_data = bsxfun(@rdivide, group_binned_data{g}, max(group_binned_data{g}, [], 2)); % normalized based on row-max
    y = mean(norm_data,1);
    plot(x, y, 'LineStyle',group_styles(g),'LineWidth',linewidth,'Color',group_colors(g,:));
    hold on;
    L(g) = plot(nan, nan, 'color', group_colors(g,:), 'LineStyle',group_styles(g),'LineWidth',linewidth);
    hold on;
end
validGroups = logical(cellfun(@(x) length(x), groupIdx,'uni',1));
legend(L(validGroups), groups(validGroups),'Location','best'); hold on;
xlabel('Voronoi Density (nm^{-2})','FontSize',24);
xtickVals = (0:1:6);
set(gca, 'XLim', [min(xtickVals) max(xtickVals)], 'XTick',xtickVals, 'XTickLabel', arrayfun(@(elem) sprintf('10^{%d}',elem), xtickVals,'uni',0));
ylabel('Probability','FontSize',24,'interpreter','latex');
xlim([min(xtickVals) max(xtickVals)]);
ylim([0,1])
title({'Empirical CDF','Average'},'Interpreter','none','FontSize',9);
hold off
if strcmp(suffix,'')
    saveas(gcf,fullfile(workDir, 'VoronoiDensityCDF.png'));
else
    saveas(gcf,fullfile(workDir, ['VoronoiDensityCDF_' suffix '.png']));
end