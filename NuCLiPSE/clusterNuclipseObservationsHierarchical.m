function clusterNuclipseObservationsHierarchical(T, varSelected)
    %% plot setup
    figure();
    TL = tiledlayout(3,6,"TileSpacing","tight");
    %% PCA
    % ground truth
    nexttile(1)
    [~,scorePCA,~,~,~] = plotPCA(T(:,[{'Group'} varSelected]));
    % tree
    ax = nexttile(3);
    [h,tt]=plotTree(ax,scorePCA);
    % cluster
    [cIdxPCA, cColorPCA] = extractClustersFromTree(ax,h,tt);
    set(gca,"XTickLabel",[]);
    % plot clustering
    nexttile(2)
    plotClusters(scorePCA,cIdxPCA,cColorPCA);
    copyLabelsLimits(nexttile(2),nexttile(1));
    %% TSNE
    % ground truth
    nexttile(4)
    [scoreTSNE, ~] = plotTSNE(T(:,[{'Group'} varSelected]));
    % tree
    ax = nexttile(6);
    [h,tt]=plotTree(ax, scoreTSNE);
    % cluster
    [cIdxTSNE, cColorTSNE] = extractClustersFromTree(ax,h,tt);
    set(gca,"XTickLabel",[]);
    % plot clustering
    nexttile(5);
    plotClusters(scoreTSNE,cIdxTSNE,cColorTSNE);
    copyLabelsLimits(nexttile(5),nexttile(4));
    %% cluster composition
    clusterComposition(cIdxPCA,cColorPCA,T,TL,[1,3],7)
    clusterComposition(cIdxTSNE,cColorTSNE,T,TL,[1,3],10)
    %% show group breakdown
    clusterGroupDivision(cIdxPCA,cColorPCA,T,TL,[1,3],13)
    clusterGroupDivision(cIdxTSNE,cColorTSNE,T,TL,[1,3],16)
end