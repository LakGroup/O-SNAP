%% plot ROC matrix
function plot_SNAP_ROC(ax,model_type,classifiers)
    classifiers_same_model = get_classifers_by_model_type(model_type,classifiers);
    % each column is a batch, row is a replicate classifier
    ROC_data_all = cellfun(@(x) x.PerformanceCurve,classifiers_same_model,'uni',0);
    groups = cellfun(@(x) x.group, ROC_data_all{1,1},'uni',1);
    n_groups = numel(groups);
    n_batch = size(ROC_data_all,2);
    n_x = 100;
    x_all = 0:1/n_x:1;
    y_all = cell(1,n_groups);
    AUC_all = zeros(n_batch,n_groups);
    for g=1:n_groups
        y_all{g} = zeros(n_x+1,n_batch);
    end
    c_data = lines(n_groups);
    for b=1:n_batch
        ROC_data_batch_all = [ROC_data_all{:,b}];
        for g=1:n_groups
            tmp = ROC_data_batch_all(g,:);
            n_tmp = numel(tmp);
            xy_all = cellfun(@(x) [x.x x.y],tmp,'uni',0);
            [~,tmp_x_longest_idx] = max(cellfun('length',xy_all));
            x = xy_all{tmp_x_longest_idx}(:,1);
            y = zeros(numel(x),n_tmp);
            for i=1:n_tmp
                [x_unique,~,x_unique_idx] = unique(xy_all{i}(:,1));
                y_unique=zeros(1,numel(x_unique));
                for j=1:numel(x_unique)
                    y_unique(j) = max(xy_all{i}(x_unique_idx==j,2));
                end
                y(:,i) = interpn(x_unique,y_unique,x);
            end
            xy = unique([x mean(y,2)],'rows');
            AUC_all(b,g) = median(cellfun(@(x) x.AUC,tmp,'uni',1));
            y_all{g}(:,b) = interp1(xy(:,1),xy(:,2),x_all);
        end        
    end
    for g=1:n_groups
        y = mean(y_all{g},2)';
        y_up = prctile(y_all{g},97.5,2)';
        y_down = prctile(y_all{g},2.5,2)';
        fill(ax,[x_all x_all(end:-1:1)],[y_up y_down(end:-1:1)],c_data(g,:),...
            'FaceAlpha',0.1,...
            'EdgeColor','none');
        hold on
        plot(ax,x_all,y_up,':','Color',c_data(g,:),'LineWidth',1);
        hold on
        plot(ax,x_all,y_down,':','Color',c_data(g,:),'LineWidth',1)
        hold on
        plot(ax,[0 x_all],[0 y],'Color',c_data(g,:),'LineWidth',2)
        hold on
    end
    legend_labels = compose("%s (AUC=%0.3f+/-%0.3f)",groups,mean(AUC_all,1)',std(AUC_all,0,1)');
    legend(ax,reshape([strings(n_groups,3) legend_labels]',1,[]),'Location','best','Interpreter','none');
    xlim(ax,[0 1])
    xlabel('False Positive Rate')
    ylim(ax,[0 1])
    ylabel('True Positive Rate')
    set(gca,"FontSize",16)
end