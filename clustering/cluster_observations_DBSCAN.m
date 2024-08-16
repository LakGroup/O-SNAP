%% get ground truth data
[~,~,g_idx]=unique(T.group);
X = table2array(normalize(T(:,var_selected)));
%% plot setup
figure()
TL = tiledlayout(2,2,"TileSpacing","tight");
%% PCA
nexttile
[~,score_PCA,~,~,~] = plot_PCA(T(:,var_selected),T.group);
% cluster
c_idx_PCA = dbscan(score_PCA,5,5);
n_clusters_PCA = length(unique(c_idx_PCA(c_idx_PCA>0)));
% plot clustering
plot_clusters(score_PCA,c_idx_PCA);
%% TSNE
nexttile
[score_tSNE, ~] = plot_tSNE(T(:,4:end),T.group);
% cluster
c_idx_tSNE = dbscan(score_tSNE,20,5);
n_clusters_tSNE = length(unique(c_idx_tSNE(c_idx_tSNE>0)));
% plot clustering
plot_clusters(score_tSNE,c_idx_tSNE);
%% show cluster composition
cluster_composition(c_idx_PCA,T,TL,[1,1],3)
cluster_composition(c_idx_tSNE,T,TL,[1,1],4)
