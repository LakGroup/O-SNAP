function cluster_composition(c_idx,cColor,T,TL,tile_span,tile_n)
n_clusters = length(unique(c_idx(c_idx>0)));
if n_clusters > 0
    TL_1 = tiledlayout(TL,1,n_clusters);
    TL_1.Layout.Tile = tile_n;
    TL_1.Layout.TileSpan = tile_span;
    for c=1:n_clusters
        in_cluster = c_idx==c;
        S = table(T.group,in_cluster,'VariableNames',["group","In_cluster"]);
        S = pivot(S,Columns="In_cluster",Rows="group");
        nexttile(TL_1)
        pie(S.true);
        cluster_title = title("Cluster " + c);
        cluster_title.Color = cColor(c,:);
        clearvars S
    end
    TL_2 = tiledlayout(TL,1,n_clusters);
    TL_2.Layout.Tile = tile_n+2*n_clusters;
    TL_2.Layout.TileSpan = [1 n_clusters];
    lgd = legend(replace(unique(T.group),"_"," "),'Orientation',"horizontal");
    lgd.Layout.Tile = 'south';
end
end
