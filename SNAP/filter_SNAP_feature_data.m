function feature_data_filtered = filter_SNAP_feature_data(feature_data,groups,replicates)
arguments
    feature_data table
    groups cell
    replicates cell
end
% filter for groups and replicates
fprintf("  Filtering feature data table for desired groups and replicates...\n")
feature_data_filtered = preprocess_SNAP_table(feature_data,...
    "groups",groups,"replicates",replicates,...
    "normalize",false,"remove_NaN",true,...
    "keep_rep_sample_info",true);
% end execution if not enough groups after filtering
if numel(unique(feature_data_filtered.group)) < 2
    ME = MException('SNAP:insufficient_group_number', ...
                    sprintf('Only one group found: %s',groups{1}));
    throw(ME)
end
end