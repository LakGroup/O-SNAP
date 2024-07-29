function clusterGroupDivision(cIdx,cColor,T,TL,tileSpan,tileN)
[groups,~,gIdx] = unique(T.Group);
nGroups = length(groups);
nClusters = length(unique(cIdx(cIdx>0)));
gColor = lines(nGroups);
if nClusters > 0
    TL_1 = tiledlayout(TL,1,nGroups);
    TL_1.Layout.Tile = tileN;
    TL_1.Layout.TileSpan = tileSpan;
    for g=1:nGroups
        inGroup = gIdx==g;
        S = table(cIdx,inGroup,'VariableNames',["Cluster","InGroup"]);
        S = pivot(S,Columns="InGroup",Rows="Cluster");
        nexttile(TL_1)
        p = pie(S.true);
        for c=1:nClusters
            p(c*2-1).FaceColor = cColor(c,:);
        end
        groupTitle = title(groups(g));
        groupTitle.Color = gColor(g,:);
        clearvars S
    end
    TL_2 = tiledlayout(TL,1,nGroups);
    TL_2.Layout.Tile = tileN+2*nGroups;
    TL_2.Layout.TileSpan = [1 nGroups];
    lgd = legend("Cluster " + (1:nClusters),'Orientation',"horizontal");
    lgd.Layout.Tile = 'south';
end
end