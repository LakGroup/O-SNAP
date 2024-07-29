function plotNuclipse(T, groups, saveDir)
% varSelMethods = ["fscc hi2","fscmrmr","relieff"];
varSelMethods = ["fscmrmr"];
idx = any(reshape(cell2mat(arrayfun(@(x) strcmpi(T.Group,x), groups,'uni',0)),[],length(groups)),2);
T = T(idx,:);
[groups, ~, groupIdx] = unique(T.Group);
nVarsSel = length(groups)+1;
for m=1:length(varSelMethods)
    % set up
    varSelMethod = varSelMethods(m);
    figure(); 
    tiledlayout('flow')
    sgtitle(varSelMethod);
    %% VARIABLE SELECTED PCA
    varSelected = selectVars(T,nVarsSel,varSelMethod);
    nexttile
    plotPCA(T(:,[{'Group'} varSelected]));
    % SHUFFLED LABELS
    groups_shuffled = T.Group(randperm(length(T.Group)));
    T_shuffled = T(:,2:end);
    T_shuffled.Group = groups_shuffled;
    varSelected_shuffled = selectVars(T_shuffled, nVarsSel,varSelMethod);
    title("Variable Score - shuffled labels")
    nexttile
    plotPCA(T_shuffled(:,[{'Group'} varSelected_shuffled]))
    title("PCA - shuffled labels")
    savefig(fullfile(saveDir, "PCA_"+varSelMethod+"_"+join(groups,'_')+".fig"));
    pause(1)
    % hist plot selected features
    figure();
    nVarFac = factor(nVarsSel);
    if(nVarFac(end) > 4)
        nVarFac = factor(nVarsSel+1);
    end
    for k = 1:nVarsSel
        subplot(nVarFac(end), prod(nVarFac(1:end-1)), k)
        for g=1:length(groups)
            minBin = min(T{:, varSelected{k}});
            maxBin = max(T{:, varSelected{k}});
            histogram(T{groupIdx==g, varSelected{k}},...
                'BinEdges', minBin:(maxBin-minBin)/15:maxBin,...
                'DisplayStyle','stairs','Normalization','probability');
            hold on
        end
        xlabel(varSelected{k})
        ylim([0 0.75])
    end
    sgtitle(varSelMethod)
    savefig(fullfile(saveDir, "Hist_"+varSelMethod+"_"+join(groups,'_')+".fig"));
    save(fullfile(saveDir, "varsSelected_" +varSelMethod+ ".mat"),"varSelected")
end
try
    %% unsupervised hierarchical clustering (observations)
    clusterNuclipseObservationsHierarchical(T, varSelected);
    savefig(fullfile(saveDir, "ClusteringObs_"+join(groups,'_')+".fig"));
end
try
    %% unsupervised hierarchical clustering (features)
    clusterNuclipseFeaturesHierarchical(T, varSelected);
    savefig(fullfile(saveDir, "ClusteringFeat_"+join(groups,'_')+".fig"));
end
end
