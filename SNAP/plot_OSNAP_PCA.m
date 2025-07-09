function plot_OSNAP_PCA(coeff,score,explained,options)
arguments
    coeff double
    score double
    explained double
    options.save_path string = ""
    options.bio_reps_info string = []
    options.groups_info string = []
    options.color_by_group logical = true
    options.show_bio_reps logical = true
    options.marker_size = 40
end
    n_comps = size(coeff,2);
    if n_comps >=3
        coeff = coeff(:,1:3);
        score = score(:,1:3);
    elseif n_comps == 2
        score = score(:,1:2);
    elseif n_comps < 2
        return
    end
    %% plot
    marker_list = {'o','^','square','v','diamond','pentagram','hexagram'};
    if ~isempty(options.groups_info)
        [groups,~,groups_idx] = unique(options.groups_info);
    end
    if options.show_bio_reps && ~isempty(options.bio_reps_info)
        [bio_reps,~,bio_reps_idx] = unique(options.bio_reps_info);
    end
    if ~options.color_by_group && exist(bio_reps,"var")
        s_face_colors = lines(length(bio_reps));
    elseif ~isempty(options.groups_info)
        s_face_colors = lines(length(groups));
    else
        s_face_colors = [0 0 0];
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
    figure('visible','off','Position', [10 10 1210 910]);
    % only group information is availble
    if ~exist("bio_reps","var") && ~isempty(groups)
        for g=1:length(groups)
            idx = groups_idx == g;
            if n_comps>=3
                s = scatter3(score(idx,1),score(idx,2),score(idx,3),options.marker_size,'filled');
            elseif n_comps==2
                s = scatter(score(idx,1),score(idx,2),options.marker_size,s_face_colors(g,:),'filled');
            elseif n_comps == 1
                s = swarmchart(ones(sum(idx),1),score(idx,1),'filled','YJitter','density');
            end
            s.MarkerFaceColor = s_face_colors(g,:);
            hold on
        end
        legend(groups,'Location','best','Interpreter','none','FontSize',14)
    % color = groups; marker type = replicate
    elseif options.color_by_group && ~isempty(groups)
        for b=1:length(bio_reps)
            for g=1:length(groups)
                idx = all([(groups_idx == g) (bio_reps_idx == b)],2);
                if n_comps >=3
                    s = scatter3(score(idx,1),score(idx,2),score(idx,3),'filled');
                elseif n_comps ==2
                    s = scatter(score(idx,1),score(idx,2),'filled');
                elseif n_comps == 1
                    s = swarmchart(ones(sum(idx),1),score(idx,1),'filled','YJitter','density');
                end
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
    elseif ~isempty(groups)
        for b=1:length(bio_reps)
            for g=1:length(groups)
                idx = all([(groups_idx == g) (bio_reps_idx == b)],2);
                if n_comps >=3
                    s = scatter3(score(idx,1),score(idx,2),score(idx,3),'filled');
                elseif n_comps ==2
                    s = scatter(score(idx,1),score(idx,2),'filled');
                elseif n_comps == 1
                    s = swarmchart(ones(sum(idx),1),score(idx,1),'filled','YJitter','density');
                end
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
        legend(bio_reps,'Location','southoutside','Interpreter','none','FontSize',18)
    else 
        if n_comps >=3
            scatter3(score(:,1),score(:,2),score(:,3),'filled');
        else
            scatter(score(:,1),score(:,2),'filled');
        end
    end
    % biplot(coeff,'VarLabels',T_num.Properties.VariableNames);
    title('PCA')
    set(gca, "XTickLabel",[]);
    set(gca, "YTickLabel",[]);
    if n_comps >=3
        set(gca, "ZTickLabel",[]);
    end
    xlabel("PC1 (" +  sprintf("%4.2f",explained(1)) + "%)","FontSize",18)
    ylabel("PC2 (" +  sprintf("%4.2f",explained(2)) + "%)","FontSize",18)
    if n_comps >=3
        zlabel("PC3 (" +  sprintf("%4.2f",explained(3)) + "%)","FontSize",18)
    end
    axis tight
    if options.save_path ~= ""
        saveas(gcf,options.save_path+".png");
        savefig(options.save_path+".fig");
    end
end