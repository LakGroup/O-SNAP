% -------------------------------------------------------------------------
% plot_OSNAP_volcano.m
% -------------------------------------------------------------------------
% Creates a volcano plot to display the magnitude and significance of
% O-SNAP feature changes between two conditions
%
% Example on how to use it:
%   plot_OSNAP_volcano(feature_compare_data,"Control","Experimental")
% -------------------------------------------------------------------------
% Input:
%   feature_compare_data: A table containing all features and the values of
%                         the means for each group, fold change, p-value, 
%                         and adj. p-value
%   group_1: String identifier for group 1
%   group_2: String identifier for group 2 
% Options:
%   alpha: Significance testing threshold
%   fold_change_threshold: Signficant fold change threshold
%   save_path: Name of file location to save to, excluding extension
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
function plot_OSNAP_volcano(feature_compare_data,group_1,group_2,options)
arguments
    feature_compare_data table;
    group_1 string;
    group_2 string;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.save_path string = "";
end
marker_size = 20;
%% Find features that satisfies significance and fold change thresholds
% Change color depending in group increase/decrease
% Change size depending on if the feature is identified as significant
log2_fold_change_threshold = log2(options.fold_change_threshold);
idx_sig = all([abs(feature_compare_data.log2_fold_change)>log2_fold_change_threshold feature_compare_data.adj_p_value < options.alpha],2);
idx_sig_up = all([feature_compare_data.log2_fold_change>0 idx_sig],2);
idx_sig_down = all([feature_compare_data.log2_fold_change<0 idx_sig],2);
scatter_size = marker_size*ones(size(feature_compare_data,1),1);
scatter_size(idx_sig) = 2*marker_size*ones(sum(idx_sig),1);
scatter_color = repmat([0.5 0.5 .5],size(feature_compare_data,1),1);
scatter_color(idx_sig_up,:) = repmat([1 0 0],sum(idx_sig_up),1);
scatter_color(idx_sig_down,:) = repmat([0 0 1],sum(idx_sig_down),1);
%% Plot
f = figure('visible','off'); hold on;
% Plot thresholds
p1=plot([log2_fold_change_threshold log2_fold_change_threshold],[NaN NaN],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2);
p2=plot(-[log2_fold_change_threshold log2_fold_change_threshold],[NaN NaN],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2);
p3=plot([NaN NaN], [-log10(options.alpha) -log10(options.alpha)],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2);
% Plot features
scatter(feature_compare_data.log2_fold_change,-log10(feature_compare_data.adj_p_value),scatter_size,scatter_color,'filled','MarkerEdgeColor','k');
% Plot labels
xlabel("Log_{2}(Fold Change)","FontSize",14)
ylabel("-Log_{10}(adj. p-value)","FontSize",14)
title(replace(group_1 + " vs " + group_2,"_","-"),'FontSize',14,'FontWeight','normal');
x_bounds = xlim();
y_bounds = ylim();
ax=gca();
set(ax,'FontSize',14)
marker_group_1 = scatter(NaN,NaN,marker_size*2,"filled","b","MarkerEdgeColor","k");
marker_group_2 = scatter(NaN,NaN,marker_size*2,"filled","r","MarkerEdgeColor","k");
legend([marker_group_1,marker_group_2],...
    sprintf("Increased in %s",replace(group_1,"_","-")),...
    sprintf("Increased in %s",replace(group_2,"_","-")),...
    "FontSize",10,"Location","Best");
% Adjust thresholds
set(ax,'YLimMode','manual')
set(ax,'XLimMode','manual')
p1.YData = y_bounds;
p2.YData = y_bounds;
p3.XData = x_bounds;
%% Save results
if options.save_path ~= ""
    save_path = sprintf("%s_%s_%s",options.save_path,group_1,group_2);
    savefig(save_path+".fig");
    exportgraphics(f,save_path+".png","Resolution",300,"BackgroundColor","w");
else
    set(f,'visible','on');
end
end

