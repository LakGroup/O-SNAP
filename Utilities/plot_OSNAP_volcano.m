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
%% Find features that satisfies significance and fold change thresholds
% Change color depending in group increase/decrease
% Change size depending on if the feature is identified as significant
log2_fold_change_threshold = log2(options.fold_change_threshold);
idx_sig = all([abs(feature_compare_data.log2_fold_change)>log2_fold_change_threshold feature_compare_data.adj_p_value < options.alpha],2);
idx_sig_up = all([feature_compare_data.log2_fold_change>0 idx_sig],2);
idx_sig_down = all([feature_compare_data.log2_fold_change<0 idx_sig],2);
scatter_size = 20*ones(size(feature_compare_data,1),1);
scatter_size(idx_sig) = 40*ones(sum(idx_sig),1);
scatter_color = repmat([0.5 0.5 .5],size(feature_compare_data,1),1);
scatter_color(idx_sig_up,:) = repmat([1 0 0],sum(idx_sig_up),1);
scatter_color(idx_sig_down,:) = repmat([0 0 1],sum(idx_sig_down),1);
%% Plot
f = figure('visible','off'); hold on;
% Plot features
scatter(feature_compare_data.log2_fold_change,-log10(feature_compare_data.adj_p_value),scatter_size,scatter_color,'filled','MarkerEdgeColor','k');
% Plot thresholds
plot([log2_fold_change_threshold log2_fold_change_threshold],ylim(),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
plot(-[log2_fold_change_threshold log2_fold_change_threshold],ylim(),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
plot(xlim(), [-log10(options.alpha) -log10(options.alpha)],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
% Plot labels
ylabel("-log_{10}(adj. p-value)","FontSize",18)
xlabel("log_{2}(Fold Change)","FontSize",18)
title(replace("\color{blue}" + group_1 + " \color{black}vs " + "\color{red}" + group_2,"_","-"),...
    'FontSize',20);
set(gca,'FontSize',14)
%% Save results
if options.save_path ~= ""
    save_path = sprintf("%s_%s_%s",options.save_path,group_1,group_2);
    savefig(save_path+".fig");
    saveas(f,save_path+".png");
end
end

