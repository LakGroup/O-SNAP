function cluster_group_division(c_idx,cColor,T,TL,tile_span,tile_n)
[groups,~,g_idx] = unique(T.group);
n_groups = length(groups);
n_clusters = length(unique(c_idx(c_idx>0)));
gColor = lines(n_groups);
if n_clusters > 0
    TL_1 = tiledlayout(TL,1,n_groups);
    TL_1.Layout.Tile = tile_n;
    TL_1.Layout.TileSpan = tile_span;
    for g=1:n_groups
        in_group = g_idx==g;
        S = table(c_idx,in_group,'VariableNames',["Cluster","In_group"]);
        S = pivot(S,Columns="In_group",Rows="Cluster");
        nexttile(TL_1)
        p = pie(S.true);
        for c=1:n_clusters
            p(c*2-1).FaceColor = cColor(c,:);
        end
        group_title = title(replace(groups(g),"_"," "));
        group_title.Color = gColor(g,:);
        clearvars S
    end
    TL_2 = tiledlayout(TL,1,n_groups);
    TL_2.Layout.Tile = tile_n+2*n_groups;
    TL_2.Layout.TileSpan = [1 n_groups];
    lgd = legend("Cluster " + (1:n_clusters),'Orientation',"horizontal");
    lgd.Layout.Tile = 'south';
end
end

