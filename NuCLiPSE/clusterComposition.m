function clusterComposition(cIdx,cColor,T,TL,tileSpan,tileN)
nClusters = length(unique(cIdx(cIdx>0)));
if nClusters > 0
    TL_1 = tiledlayout(TL,1,nClusters);
    TL_1.Layout.Tile = tileN;
    TL_1.Layout.TileSpan = tileSpan;
    for c=1:nClusters
        inCluster = cIdx==c;
        S = table(T.Group,inCluster,'VariableNames',["Group","InCluster"]);
        S = pivot(S,Columns="InCluster",Rows="Group");
        nexttile(TL_1)
        pie(S.true);
        clusterTitle = title("Cluster " + c);
        clusterTitle.Color = cColor(c,:);
        clearvars S
    end
    TL_2 = tiledlayout(TL,1,nClusters);
    TL_2.Layout.Tile = tileN+2*nClusters;
    TL_2.Layout.TileSpan = [1 nClusters];
    lgd = legend(unique(T.Group),'Orientation',"horizontal");
    lgd.Layout.Tile = 'south';
end
end