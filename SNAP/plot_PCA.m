function [coeff,score,latent,tsquared,explained] = plot_PCA(T,options)
arguments
    T table
    options.color_by_group logical = true
    options.show_bio_reps logical = true
    options.marker_size = 20
end
    %% normalize T
    groups = T.group;
    T_num = T(:,vartype('numeric'));
    X = table2array(normalize(T_num));
    [groups, ~, group_idx] = unique(groups);
    %% perform PCA
    [coeff,score,latent,tsquared,explained] = pca(X,...
        'Rows','pairwise','NumComponents',3);
    %% plot
    marker_list = {'o','^','square','v','diamond','pentagram','hexagram'};
    if options.show_bio_reps && contains('biological_replicate',T.Properties.VariableNames)
        [bio_reps,~,bio_reps_idx] = unique(T.biological_replicate);
    end
    if ~options.color_by_group && exist(bio_reps,"var")
        s_face_colors = lines(length(bio_reps));
    else
        s_face_colors = lines(length(groups));
    end
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
    % only group information is availble
    if ~exist("bio_reps","var")
        for g=1:length(groups)
            idx = group_idx == g;
            s = scatter3(score(idx,1),score(idx,2),score(idx,3),options.marker_size,'filled');
            s.MarkerFaceColor = s_face_colors(g,:);
            hold on
        end
        legend(groups,'Location','southoutside','Interpreter','none')
    % color = groups; marker type = replicate
    elseif options.color_by_group
        for b=1:length(bio_reps)
            for g=1:length(groups)
                idx = all([(group_idx == g) (bio_reps_idx == b)],2);
                s = scatter3(score(idx,1),score(idx,2),score(idx,3),'filled');
                s.MarkerFaceColor = s_face_colors(g,:);
                s.Marker = marker_list(b);
                if b == 1
                    s.SizeData = options.marker_size;
                else
                    s.SizeData = 2*options.marker_size;
                end
                hold on
            end
        end
        legend(groups,'Location','southoutside','Interpreter','none')
    % color = replicate; marker type = color
    else
        for b=1:length(bio_reps)
            for g=1:length(groups)
                idx = all([(group_idx == g) (bio_reps_idx == b)],2);
                s = scatter3(score(idx,1),score(idx,2),score(idx,3),'filled');
                s.MarkerFaceColor = s_face_colors(b,:);
                s.Marker = marker_list(g);
                if g == 1
                    s.SizeData = options.marker_size;
                else
                    s.SizeData = 2*options.marker_size;
                end
                hold on
            end
        end
        legend(bio_reps,'Location','southoutside','Interpreter','none')
    end
    % biplot(coeff,'VarLabels',T_num.Properties.VariableNames);
    title('PCA')
    set(gca, "XTickLabel",[]);
    set(gca, "YTickLabel",[]);
    set(gca, "ZTickLabel",[]);
    xlabel("PC1 (" +  sprintf("%4.2f",explained(1)) + "%)")
    ylabel("PC2 (" +  sprintf("%4.2f",explained(2)) + "%)")
    zlabel("PC3 (" +  sprintf("%4.2f",explained(3)) + "%)")
end