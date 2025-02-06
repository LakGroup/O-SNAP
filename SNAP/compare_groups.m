function feature_comparisons = compare_groups(feature_data,options)
arguments
    feature_data table
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.plot logical = true;
    options.plot_labels logical = true;
    options.save_path string = "";
end
%% setup
groups = unique(feature_data.group);
group_pairs = nchoosek(groups,2);
n_group_pairs = size(group_pairs,1);
feature_comparisons = cell(1,n_group_pairs);

%% compare each combo of group pairs
for i=1:n_group_pairs
    feature_comparisons{i}.groups = group_pairs(i,:);
    feature_comparisons{i}.feature_table = compare_group_pair(feature_data,group_pairs(i,1),group_pairs(i,2),...
        "alpha",options.alpha,...
        "fold_change_threshold",options.fold_change_threshold,...
        "plot",options.plot,...
        "plot_labels",options.plot_labels,...
        "save_path",options.save_path);
end
end