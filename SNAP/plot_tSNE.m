function [score,loss] = plot_tSNE(T)
    rng default % for reproducibility
    %% normalize T
    groups = T.group;
    T_num = T(:,vartype('numeric'));
    X = table2array(normalize(T_num));
    [groups, ~, group_idx] = unique(groups);
    %% perform TSNE
    [score,loss] = tsne(X,'NumDimensions',3);
    %% plot
    sFaceColors = lines(length(groups));
    for g=1:length(groups)
        idx = group_idx == g;
        s = scatter3(score(idx,1),score(idx,2),score(idx,3),10,'filled');
        s.MarkerFaceColor = sFaceColors(g,:);
        hold on
    end
    legend(groups,'Location','southoutside','Interpreter','none')
    title('t-SNE')
    set(gca, "XTickLabel",[]);
    set(gca, "YTickLabel",[]);
    set(gca, "ZTickLabel",[]);
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
end

