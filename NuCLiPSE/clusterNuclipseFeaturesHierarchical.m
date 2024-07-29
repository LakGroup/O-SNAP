function clusterNuclipseFeaturesHierarchical(T, varSelected)
    %% get ground truth data
    X = table2array(normalize(T(:,4:end)));
    featureNames = T.Properties.VariableNames(:,4:end);
    %% correlation
    R = abs(corrcoef(X,'Rows','pairwise'));
    %% plot setup
    figure()
    tree = linkage(R,'ward');
    D = pdist(R);
    leafOrder = optimalleaforder(tree,D);
    dendrogram(tree,0,'Reorder',leafOrder,'ColorThreshold','default');
    ax = gca;
    featureLabels = featureNames(str2num(ax.XTickLabel));
    set(ax,"XTickLabel",featureLabels);
    varSelIdx = arrayfun(@(x) strcmpi(x,ax.XTickLabel), varSelected, 'uni',0);
    varSelIdx = [varSelIdx{:}];
    varSelIdx = any(varSelIdx,2);
    ax.XTickLabel(varSelIdx) = cellfun(@(x) ['\bf{' upper(x) '}'], ax.XTickLabel(varSelIdx), 'uni',0);
    title({'Hierarchical clustering of features','|R|'})
end