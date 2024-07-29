function percentileData = plotVoronoiDensityPercentiles(workDir, groups, voronoiData, percentileData)
%% set up
n_groups = length(groups);
% create list of percentiles
percentileFields = fieldnames(percentileData{1});
percentileFields = percentileFields(3:end);
percentiles = cellfun(@(x) str2double(extractAfter(x,'prctile_')),percentileFields);
n_percentiles = length(percentiles);
% each cell of y is a different groupi's percentile data 
% row: each sample // col: percentiles
y_percentile = cell(1,n_groups);
yl_min = -50;
yl_max = 150;
% get data
for g=1:n_groups
    g_y = strcmp(groups{g}, cellfun(@(x) x.group, percentileData, 'uni', 0));
    nSampInGroup = length(percentileData{g_y}.sampleNames);
    if nSampInGroup > 0
        y_percentile{g_y} = cell2mat(cellfun(@(x) (percentileData{g_y}.(x)), percentileFields, 'uni',0))';
    end
end
[~, groupIdx] = sortSampleToGroups(groups,voronoiData);

%% plot (k=1: w/ outliers; k=2: without outliers)
for k=1:2
    if k==1
        savePath = fullfile(workDir,'PercentileValues.png');
        titleName = 'Log Voronoi Density Percentiles';
    else
        savePath = fullfile(workDir,'PercentileValues_outliersRemoved.png');
        titleName = {'Log Voronoi Density Percentiles','(Outliers Removed)'};
    end
    % prepare percentile plot
    f_percentile = figure('visible','off');
    f_percentile.Position = [0 0 1680 920];
    c = lines(n_groups);
    % offset for x positions so that data between groups does not overlap
    x_dif = 0.8*(1:n_groups)/(n_groups); x_dif= x_dif-mean(x_dif);
    % plot the values
    for g=1:n_groups
        if size(groupIdx{g},2) < 1
            continue;
        end
        % plot every sample
        g_y = strcmp(groups{g}, cellfun(@(x) x.group, percentileData, 'uni', 0));
        x = reshape(repmat(1:n_percentiles,size(groupIdx{g},2),1).',1,[]) + x_dif(g);
        if k==1
            for i=1:n_percentiles
                y_percentile{g_y}(:,i) = filloutliers(percentileData{g_y}.(percentileFields{i}),NaN);
            end
        else
            for i=1:n_percentiles
                    y_percentile{g_y}(:,i) = percentileData{g_y}.(percentileFields{i});
            end
        end
        plot(x, reshape(y_percentile{g_y}.',1,[]), '.', 'MarkerSize', 7, 'Color',c(g,:),'DisplayName',char(groups(g)));
        hold on;
        % plot the mean of the group
        if contains(groups(g),"Pos")
            p = plot((1:n_percentiles) + x_dif(g), mean(y_percentile{g_y},1,'omitnan'), '+', 'MarkerSize', 10, 'LineWidth', 3,'Color',c(g,:));
        else
            p = plot((1:n_percentiles) + x_dif(g), mean(y_percentile{g_y},1,'omitnan'), 'x', 'MarkerSize', 10, 'LineWidth', 3,'Color',c(g,:));
        end
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    % finish plot
    legend('Location','northwest','FontSize',18);
    ax = gca;
    ax.XTick = 1:n_percentiles;
    ax.XTickLabels = arrayfun(@num2str, percentiles, 'UniformOutput', 0);
    ax.XAxis.FontSize = 16;
    ax.YTickLabels = arrayfun(@(elem) sprintf('10^{%d}',elem), ax.YTick,'uni',0);
    ax.YAxis.FontSize = 16;
    title(titleName,'FontSize',36)
    xlabel('Percentile','FontSize',24);
    ylabel("Voronoi Density (nm^{-2})",'FontSize',24,'interpreter','tex');
    grid on;
    saveas(gcf,savePath);
end
%% plot barplots
percentile_thresh = [40, 50, 80];
if all(cellfun(@(x) ~isempty(x), y_percentile,'uni',1))
    %% find diff_density between LowColVI and HighColVI
    group_a = 'LowColVI';
    group_b = 'HighColVI';
    delta_density_P = zeros(length(percentile_thresh));
    for i=1:length(percentile_thresh)
        a_data = percentileData{strcmp(cellfun(@(elem) elem.group, percentileData, 'uni',0),group_a)};
        b_data = percentileData{strcmp(cellfun(@(elem) elem.group, percentileData, 'uni',0),group_b)};
        a_val = arrayfun(@(elem) 10^(elem), a_data.(['prctile_' num2str(percentile_thresh(i))]));
        b_val = arrayfun(@(elem) 10^(elem), b_data.(['prctile_' num2str(percentile_thresh(i))]));
        a_mean = mean(a_val);
        b_mean = mean(b_val);
        delta_density_P(i) = sign(b_mean - a_mean)*(b_mean)*100/a_mean;
    end
    for i=1:size(delta_density_P,1)
        figure('visible','off');
        bar(delta_density_P(i));
        ylim([yl_min yl_max]);
        ylabel('\Delta Voronoi Polygon density (%)','Interpreter','tex');
        title({'\Delta Voronoi Density (High ColVI/Low ColVI)',[num2str(percentile_thresh(i)) '^{th} Percentile']},'FontSize',14,'Interpreter','tex')
        set(gca,'Units','normalized','OuterPosition',[0 0 1 1]);
        saveas(gcf,fullfile(workDir,['VoronoiDensityDiff_ColVI_' num2str(percentile_thresh(i)) '.png']));
    end
    % %% find diff_density between P(N) and P0
    % p0 = 'P0';
    % pn = {'P2','P6'};
    % treatments = {'CMNeg','CMPos'};
    % col = reshape(compose('%s-%s', string(pn)',string(treatments)),1,length(treatments)*length(string(pn)));
    % delta_density_P = zeros(length(percentile_thresh),size(col,2));
    % for i=1:length(percentile_thresh)
    %     for j=1:size(col,2)
    %         group_a = [p0 '-' treatments{ceil(j/size(pn,2))}];
    %         group_b = col{j};
    %         a_data = percentileData{strcmp(cellfun(@(elem) elem.group, percentileData, 'uni',0),group_a)};
    %         b_data = percentileData{strcmp(cellfun(@(elem) elem.group, percentileData, 'uni',0),group_b)};
    %         a_val = arrayfun(@(elem) 10^(elem), a_data.(['prctile_' num2str(percentile_thresh(i))]));
    %         b_val = arrayfun(@(elem) 10^(elem), b_data.(['prctile_' num2str(percentile_thresh(i))]));
    %         a_mean = mean(a_val);
    %         b_mean = mean(b_val);
    %         % disp([group_a ': ' num2str(a_mean) '   ' group_b ': ' num2str(b_mean) ' = ' num2str((b_mean - a_mean)*100/a_mean)])
    %         delta_density_P(i,j) = sign(b_mean - a_mean)*(b_mean)*100/a_mean;
    %     end
    % end
    % for i=1:size(delta_density_P,1)
    %     figure('visible','off');
    %     c = lines(4);
    %     c = c(2:4,:);
    %     x = categorical(treatments);
    %     b = bar(x, reshape(delta_density_P(i,:),numel(pn),[])');
    %     for j=1:length(b)
    %         b(j).FaceColor = c(j,:);
    %     end
    %     % yl_max = max(abs(ylim));
    %     % yl_min = -yl_max
    %     ylim([yl_min yl_max]);
    %     ylabel('\Delta Voronoi Polygon density (%)','Interpreter','tex');
    %     title({'\Delta Voronoi Density (P_N/P_0)',[num2str(percentile_thresh(i)) '^{th} Percentile']},'FontSize',14,'Interpreter','tex')
    %     legend(pn,'Location','northwest'); 
    %     set(gca,'Units','normalized','OuterPosition',[0 0 1 1]);
    %     saveas(gcf,fullfile(workDir,['VoronoiDensityDiff_passage_' num2str(percentile_thresh(i)) '.png']));
    % end
    % %% find diff_density between CMNeg and CMPos
    % pn = {'P0','P2','P6'};
    % treatments = {'CMNeg','CMPos'};
    % delta_density_CM = zeros(length(percentile_thresh),size(pn,2));
    % for i=1:length(percentile_thresh)
    %     for j=1:size(pn,2)
    %         group_a = [pn{j} '-' treatments{1}];
    %         group_b = [pn{j} '-' treatments{2}];
    %         a_data = percentileData{strcmp(cellfun(@(elem) elem.group, percentileData, 'uni',0),group_a)};
    %         b_data = percentileData{strcmp(cellfun(@(elem) elem.group, percentileData, 'uni',0),group_b)};
    %         a_val = arrayfun(@(elem) 10^(elem), a_data.(['prctile_' num2str(percentile_thresh(i))]));
    %         b_val = arrayfun(@(elem) 10^(elem), b_data.(['prctile_' num2str(percentile_thresh(i))]));
    %         a_mean = mean(a_val);
    %         b_mean = mean(b_val);
    %         % disp([group_a ': ' num2str(a_mean) '   ' group_b ': ' num2str(b_mean) ' = ' num2str((b_mean - a_mean)*100/a_mean)])
    %         delta_density_CM(i,j) = sign(b_mean - a_mean)*(b_mean)*100/a_mean;
    %     end
    % end
    % for i=1:size(delta_density_CM,1)
    %     figure('visible','off');
    %     c = lines(4);
    %     x = categorical(pn);
    %     b = bar(x, delta_density_CM(i,:));
    %     b.FaceColor = 'flat';
    %     for j=1:length(x)
    %         b.CData(j,:) = c(j,:);
    %     end
    %     % yl_max = max(abs(ylim));
    %     % yl_min = -yl_max
    %     ylim([yl_min yl_max]);
    %     ylabel('\Delta Voronoi Polygon density (%)','Interpreter','tex');
    %     title({'\Delta Voronoi Density (CM+/CM-)', [num2str(percentile_thresh(i)) '^{th} Percentile']},'FontSize',14,'Interpreter','tex');
    %     set(gca,'Units','normalized','OuterPosition',[0 0 1 1]);
    %     saveas(gcf,fullfile(workDir,['VoronoiDensityDiff_treatment_' num2str(percentile_thresh(i)) '.png']));
    % end
end
%% plot violin plot
for i=1:length(percentile_thresh)
    plotPercentileViolinPlot(workDir, percentileData, percentile_thresh(i))
end
end