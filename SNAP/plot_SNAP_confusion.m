%% plot confusion matrix
function plot_SNAP_confusion(ax,model_type,classifiers)
    classifiers_same_model = get_classifers_by_model_type(model_type,classifiers);
    confusion_data_all = cellfun(@(x) confusionmat(string(x.TestResponse),string(x.TestPredictions)),classifiers_same_model,'uni',0);
    n_model = size(confusion_data_all,1);
    n_batch = size(confusion_data_all,2);
    n_groups = size(confusion_data_all{1});
    c_data = cell(n_model,n_batch);
    c_data_batch = zeros(n_groups(1),n_groups(2),n_batch);
    for j=1:n_batch
        for i=1:n_model
            c_data{i,j} = confusion_data_all{i,j}./sum(confusion_data_all{i,j},2); 
        end
        c_data_batch(:,:,j) = mean(reshape(cell2mat(c_data(:,j)'),n_groups(1),n_groups(2),[]),3);
    end
    confusion_data = reshape(c_data_batch,n_groups(1),n_groups(2),[]);
    confusion_rowsum = sum(confusion_data,2);
    confusion_data = confusion_data./confusion_rowsum;
    confusion_mean = mean(confusion_data,3);
    confusion_std = std(confusion_data,0,3);
    groups = sort(unique(classifiers_same_model{1}.TestResponse));
    n_groups = numel(groups);
    set(ax,'ydir','reverse')
    xlim(ax,[0 n_groups])
    set(ax,'TickLabelInterpreter','none')
    ylim(ax,[0 n_groups])
    xticks(ax,(1:n_groups)-0.5)
    yticks(ax,(1:n_groups)-0.5)
    ax.FontSize = 18;
    groups = replace(groups,"Ctrl","Control");
    xticklabels(ax, replace(groups,"_","-"));
    yticklabels(ax, replace(groups,"_","-"));
    ytickangle(ax,90);
    xlabel(ax,'Predicted Class','FontSize',18)
    ylabel(ax,'True Class','FontSize',18,'Rotation',90)
    color_diagonal = [0.00,0.45,0.74];
    color_offdiagonal = [0.85,0.33,0.10];
    title(ax,replace(model_type,"_"," "),'FontSize',12)
    for i_true=1:n_groups
        for i_pred=1:n_groups
            c_value = confusion_mean(i_true,i_pred);
            if i_true == i_pred
                c = color_diagonal+(1-color_diagonal)*(1-c_value);
            else
                c = color_offdiagonal+(1-color_offdiagonal)*(1-c_value);
            end
            % fprintf('i_pred:%.0f  i_true:%.0f  mean:%.4f std:%.4f\n',i_pred,i_true,c_value,confusion_std(i_true,i_pred))
            if c_value > 0.5
                t = 'w';
            else
                t = 'k';
            end
            rectangle(ax,...
                'Position',[i_pred-1 i_true-1 1 1],...
                'FaceColor',c)
            % text(ax,i_pred-0.5,i_true-0.5-0.08,...
            %     sprintf('%.0f \x00B1 %.0f',...
            %     confusion_mean(i_pred,i_true),...
            %     confusion_std(i_pred,i_true)),...
            %     'HorizontalAlignment','center',...
            %     'FontSize',12,...
            %     'Color',t);
            text(ax,i_pred-0.5,i_true-0.5,...
                sprintf('%.2f%%\n\x00B1\n%.2f%%',...
                c_value*100,...
                confusion_std(i_true,i_pred)*100),...
                'HorizontalAlignment','center',...
                'FontSize',18,...
                'Color',t);
        end
    end
    disableDefaultInteractivity(ax)
end