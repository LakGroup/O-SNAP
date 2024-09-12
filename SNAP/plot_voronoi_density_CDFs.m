function plot_voronoi_density_CDFs(work_dir, groups, group_styles, group_colors, voronoi_data, suffix, nbins)
%% set up for cdf plots
n_groups = length(groups);
n_samples = length(voronoi_data);
[~, group_idx] = sort_sample_to_groups(groups,voronoi_data);
name2color = zeros(n_samples,3);
name2style = cell(n_samples);
linewidth = 2;
group_binned_data = cell(n_groups,1);
% determine min/max for x-axis
min_val_cdf = min(cellfun(@(x) prctile(x.reduced_log_voronoi_density, 1), voronoi_data, 'uni',1));
max_val_cdf = max(cellfun(@(x) prctile(x.reduced_log_voronoi_density, 99.5), voronoi_data, 'uni',1));
for g=1:n_groups
    n_samp_in_group = length(group_idx{g});
    if n_samp_in_group > 0
        group_binned_data{g} = zeros(n_samp_in_group,nbins);
        for i=1:n_samp_in_group
            % organize styles
            idx = group_idx{g}(i);
            name2color(idx,:) = repmat(group_colors(g,:),length(idx),1);
            name2style{idx} = group_styles{g};
            % fill a 2D array containing the binned log voronoi density of each sample in each row
            group_binned_data{g}(i,:) = histcounts(voronoi_data{idx}.reduced_log_voronoi_density,'Normalization','cumcount','BinEdges',linspace(min_val_cdf,max_val_cdf,nbins+1));
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
    n_samp_in_group = length(group_idx{g});
    if n_samp_in_group < 1
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
valid_groups = logical(cellfun(@(x) length(x), group_idx,'uni',1));
legend(L(valid_groups), groups(valid_groups),'Location','best'); hold on;
xlabel('Voronoi density (nm^{-2})','FontSize',24);
xtick_vals = (0:1:6);
set(gca, 'XLim', [min(xtick_vals) max(xtick_vals)], 'XTick',xtick_vals, 'XTickLabel', arrayfun(@(elem) sprintf('10^{%d}',elem), xtick_vals,'uni',0));
ylabel('Probability','FontSize',24,'interpreter','latex');
xlim([min(xtick_vals) max(xtick_vals)]);
ylim([0,1])
title({'Empirical CDF','Average'},'Interpreter','none','FontSize',9);
hold off
if strcmp(suffix,'')
    saveas(gcf,fullfile(work_dir, 'voronoi_density_CDF.png'));
else
    saveas(gcf,fullfile(work_dir, ['voronoi_density_CDF_' suffix '.png']));
end

