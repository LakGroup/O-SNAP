function [coeff,score,latent,tsquared,explained] = plot_PCA(T)
    %% normalize T
    groups = T.group;
    T_num = T(:,vartype('numeric'));
    X = table2array(normalize(T_num));
    [groups, ~, group_idx] = unique(groups);
    %% perform PCA
    [coeff,score,latent,tsquared,explained] = pca(X,...
        'Rows','pairwise','NumComponents',3);
    %% plot
    sFaceColors = lines(length(groups));
    % Force each column of the coefficients to have a positive largest element.
    % This tends to put the large var vectors in the top and right halves of
    % the plot.
    [p,d] = size(coeff);
    [~,maxind] = max(abs(coeff),[],1);
    colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
    coeff = coeff .* colsign;
    % Scale the scores so they fit on the plot, and change the sign of
    % their coordinates according to the sign convention for the coefs.
    max_coef_len = sqrt(max(sum(coeff.^2,2)));
    score = (max_coef_len.*(score ./ max(abs(score(:))))).*colsign;
    for g=1:length(groups)
        idx = group_idx == g;
        s = scatter3(score(idx,1),score(idx,2),score(idx,3),10,'filled');
        s.MarkerFaceColor = sFaceColors(g,:);
        hold on
    end
    % biplot(coeff,'VarLabels',T_num.Properties.VariableNames);
    legend(groups,'Location','southoutside','Interpreter','none')
    title('PCA')
    set(gca, "XTickLabel",[]);
    set(gca, "YTickLabel",[]);
    set(gca, "ZTickLabel",[]);
    xlabel("PC1 (" +  sprintf("%4.2f",explained(1)) + "%)")
    ylabel("PC2 (" +  sprintf("%4.2f",explained(2)) + "%)")
    zlabel("PC3 (" +  sprintf("%4.2f",explained(3)) + "%)")
end

