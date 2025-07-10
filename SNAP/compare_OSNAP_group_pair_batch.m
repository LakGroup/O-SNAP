% -------------------------------------------------------------------------
% compare_OSNAP_group_pair_batch.m
% -------------------------------------------------------------------------
% Performs a pair-wise, fold-change analysis for all combinations of
% phenotype pairs present in the O-SNAP feature table.
%
% Example on how to use it:
%   feature_comparisons = compare_OSNAP_group_pair_batch(feature_data)
% -------------------------------------------------------------------------
% Input:
%   feature_data: O-SNAP feature table
% Output:
%   feature_comparisons: A cell array where each cell contains a struct 
%                        array on the pair-wise feature comparison results
%            - groups: The two phenotypes being compared
%            - feature_table: Table with the corresponding pair-wise
%              fold-change results
% Options:
%   alpha: Significance testing threshold
%   p_adj_method: Identifier for the adjustment method for
%                 multiple hypothesis testing
%   fold_change_threshold: Signficant fold change threshold
%   plot: Flag to plot results
%   save_path: Name of file location to save to, excluding extension
%   remove_outliers: Flag to consider or remove outliers when calculating
%                    fold changes
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function feature_comparisons = compare_OSNAP_group_pair_batch(feature_data,options)
arguments
    feature_data table
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.plot logical = true;
    options.save_path string = "";
    options.remove_outliers logical = false;
end
%% setup
groups = unique(feature_data.group);
group_pairs = nchoosek(groups,2);
n_group_pairs = size(group_pairs,1);
feature_comparisons = cell(1,n_group_pairs);
%% compare each combo of group pairs
for i=1:n_group_pairs
    feature_comparisons{i}.groups = group_pairs(i,:);
    feature_comparisons{i}.feature_table = compare_OSNAP_group_pair(feature_data,group_pairs(i,1),group_pairs(i,2),...
        "alpha",options.alpha,...
        "fold_change_threshold",options.fold_change_threshold,...
        "plot",options.plot,...
        "save_path",options.save_path, ...
        "remove_outliers",options.remove_outliers);
end
end
