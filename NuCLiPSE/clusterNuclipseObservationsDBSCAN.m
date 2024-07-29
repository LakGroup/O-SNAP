%% get ground truth data
[~,~,gIdx]=unique(T.Group);
X = table2array(normalize(T(:,varSelected)));
%% plot setup
figure()
TL = tiledlayout(2,2,"TileSpacing","tight");
%% PCA
nexttile
[~,scorePCA,~,~,~] = plotPCA(T(:,varSelected),T.Group);
% cluster
cIdxPCA = dbscan(scorePCA,5,5);
nClustersPCA = length(unique(cIdxPCA(cIdxPCA>0)));
% plot clustering
plotClusters(scorePCA,cIdxPCA);
%% TSNE
nexttile
[scoreTSNE, ~] = plotTSNE(T(:,4:end),T.Group);
% cluster
cIdxTSNE = dbscan(scoreTSNE,20,5);
nClustersTSNE = length(unique(cIdxTSNE(cIdxTSNE>0)));
% plot clustering
plotClusters(scoreTSNE,cIdxTSNE);
%% show cluster composition
clusterComposition(cIdxPCA,T,TL,[1,1],3)
clusterComposition(cIdxTSNE,T,TL,[1,1],4)