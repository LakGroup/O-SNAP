function generate_SNAP_GS(feature_IDs,feature_universe,options)
arguments
    feature_IDs cell
    feature_universe struct
    options.filepath char = ['SNAP_' char(feature_universe.name) '_GS.xlsx']
end
    filepath = char(options.filepath);
    feature_set_names = fieldnames(feature_universe.universe);
    n_feature_sets = length(feature_set_names);
    feature_idxs = cell(n_feature_sets,1);
    for s=1:n_feature_sets
        idx = cellfun(@(x) strcmpi(feature_IDs,x), feature_universe.universe.(feature_set_names{s}), 'uni',0);
        feature_idxs{s} = find(any(cell2mat(idx), 2))';
    end
    max_n_features = max(cellfun(@(x) length(x), feature_idxs, 'uni',1));
    data_GS = table('Size',[n_feature_sets,2+max_n_features],...
        'VariableTypes',[{'string','string'}  repmat({'double'},1,max_n_features)]);
    for s=1:n_feature_sets
        data_GS{s,1} = feature_set_names(s);
        data_GS{s,2} = NaN;
        data_GS{s,3:max_n_features+2} = [feature_idxs{s} NaN(1,max_n_features - length(feature_idxs{s}))];
    end
    writetable(data_GS,filepath,"WriteVariableNames",false);
    update_GS(filepath)
end
