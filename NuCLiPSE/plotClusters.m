function plotClusters(points,cIdx,cColor)
nClusters = length(unique(cIdx(cIdx>0)));
%% plot clusters
for c=1:nClusters
    idx = cIdx == c;
    scatter3(points(idx,1),points(idx,2),points(idx,3),10,'filled',...
        "MarkerFaceColor",cColor(c,:),...
        "DisplayName","Cluster " + c);
    hold on
end
legend("Clusters" + (1:nClusters),'Location','southoutside')
title("Clusters")
set(gca, "XTickLabel",[]);
set(gca, "YTickLabel",[]);
set(gca, "ZTickLabel",[]);
end