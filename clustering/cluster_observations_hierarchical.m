function cluster_observations_hierarchical(T, var_selected)
    %% plot setup
    figure('Position', [10 10 1210 910])
    TL = tiledlayout(3,6,"TileSpacing","tight");
    %% PCA
    % ground truth
    nexttile(1)
    [~,score_PCA,~,~,~] = plot_PCA(T(:,[{'group'} var_selected]));
    % tree
    ax = nexttile(3);
    [h,tt]=plot_tree(ax,score_PCA);
    % cluster
    [c_idx_PCA, cColor_PCA] = extract_clusters_from_tree(ax,h,tt);
    set(gca,"XTickLabel",[]);
    % plot clustering
    nexttile(2)
    plot_clusters(score_PCA,c_idx_PCA,cColor_PCA);
    copy_labels_limits(nexttile(2),nexttile(1));
    %% TSNE
    % ground truth
    nexttile(4)
    [score_tSNE, ~] = plot_tSNE(T(:,[{'group'} var_selected]));
    % tree
    ax = nexttile(6);
    [h,tt]=plot_tree(ax, score_tSNE);
    % cluster
    [c_idx_tSNE, cColor_tSNE] = extract_clusters_from_tree(ax,h,tt);
    set(gca,"XTickLabel",[]);
    % plot clustering
    nexttile(5);
    plot_clusters(score_tSNE,c_idx_tSNE,cColor_tSNE);
    copy_labels_limits(nexttile(5),nexttile(4));
    %% cluster composition
    cluster_composition(c_idx_PCA,cColor_PCA,T,TL,[1,3],7)
    cluster_composition(c_idx_tSNE,cColor_tSNE,T,TL,[1,3],10)
    %% show group breakdown
    cluster_group_division(c_idx_PCA,cColor_PCA,T,TL,[1,3],13)
    cluster_group_division(c_idx_tSNE,cColor_tSNE,T,TL,[1,3],16)
end
