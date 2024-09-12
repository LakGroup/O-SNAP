function plot_clusters(points,c_idx,cColor)
n_clusters = length(unique(c_idx(c_idx>0)));
%% plot clusters
for c=1:n_clusters
    idx = c_idx == c;
    scatter3(points(idx,1),points(idx,2),points(idx,3),10,'filled',...
        "MarkerFaceColor",cColor(c,:),...
        "DisplayName","Cluster " + c);
    hold on
end
legend("Clusters" + (1:n_clusters),'Location','southoutside')
title("Clusters")
set(gca, "XTickLabel",[]);
set(gca, "YTickLabel",[]);
set(gca, "ZTickLabel",[]);
end

