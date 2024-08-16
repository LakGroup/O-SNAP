function cluster_features_hierarchical(T, var_selected)
    %% get ground truth data
    X = table2array(normalize(T(:,4:end)));
    feature_names = T.Properties.VariableNames(:,4:end);
    %% correlation
    R = abs(corrcoef(X,'Rows','pairwise'));
    % remove nans
    nan_idx = all(isnan(R));
    R = R(~nan_idx,~nan_idx);
    feature_names = feature_names(~nan_idx);
    %% plot setup
    figure('Position', [10 10 1210 910]);
    tree = linkage(R,'ward');
    D = pdist(R);
    leaf_order = optimalleaforder(tree,D);
    dendrogram(tree,0,'Reorder',leaf_order,'ColorThreshold','default');
    ax = gca;
    featureLabels = replace(feature_names(str2num(ax.XTickLabel)),"_"," ");
    set(ax,"XTickLabel",featureLabels);
    var_sel_idx = arrayfun(@(x) strcmpi(x,ax.XTickLabel), var_selected, 'uni',0);
    var_sel_idx = [var_sel_idx{:}];
    var_sel_idx = any(var_sel_idx,2);
    ax.XTickLabel(var_sel_idx) = cellfun(@(x) ['\bf{' upper(x) '}'], ax.XTickLabel(var_sel_idx), 'uni',0);
    title({'Hierarchical clustering of features','|R|'})
end
