function generate_OSNAP_feature_set(feature_IDs,feature_universe)
arguments
    feature_IDs cell
    feature_universe struct
end
    filepath = feature_universe.name+".xlsx";
    feature_set_names = fieldnames(feature_universe.feature_set);
    n_feature_sets = length(feature_set_names);
    feature_idxs = cell(n_feature_sets,1);
    for s=1:n_feature_sets
        idx = cellfun(@(x) find(strcmpi(feature_IDs,x)),feature_universe.feature_set.(feature_set_names{s}),'uni',1);
        dir = feature_universe.direction.(feature_set_names{s});
        feature_idxs{s} = (idx.*dir)';
    end
    max_n_features = max(cellfun(@(x) length(x), feature_idxs, 'uni',1));
    data_feature_set = table('Size',[n_feature_sets,2+max_n_features],...
        'VariableTypes',[{'string','string'}  repmat({'double'},1,max_n_features)]);
    for s=1:n_feature_sets
        data_feature_set{s,1} = feature_set_names(s);
        data_feature_set{s,2} = {NaN};
        data_feature_set{s,3:max_n_features+2} = [feature_idxs{s} NaN(1,max_n_features - length(feature_idxs{s}))];
    end
    writetable(data_feature_set,filepath,"WriteVariableNames",false,"WriteMode","overwritesheet");
    update_FS(char(filepath))
end
