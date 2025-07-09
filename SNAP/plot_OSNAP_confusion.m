%% plot confusion matrix
function plot_OSNAP_confusion(ax,groups,confusion_counts,options)
arguments
    ax
    groups string
    confusion_counts double
    options.model_type_name = ""
end
    n_groups = numel(groups);
    n_batch = size(confusion_counts,3);
    row_counts = sum(confusion_counts,2);
    confusion_normalized = zeros(n_groups,n_groups,n_batch);
    confusion_mean = zeros(n_groups,n_groups);
    confusion_std = zeros(n_groups,n_groups);

    for j=1:n_batch
        confusion_normalized(:,:,j) = bsxfun(@rdivide, confusion_counts(:,:,j), sum(confusion_counts(:,:,j),2));
    end
    for j=1:n_batch
        for i_true=1:n_groups
            row_count = row_counts(i_true,:);
            for i_pred=1:n_groups
                confusion_mean(i_true,i_pred) = mean(confusion_normalized(i_true,i_pred,:),"Weights",row_count);
                if n_batch > 1
                    confusion_std(i_true,i_pred) = std(confusion_normalized(i_true,i_pred,:),row_count);
                end
            end
        end
    end

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
    if options.model_type_name ~= ""
        title(ax,replace(options.model_type_name,"_"," "),'FontSize',12)
    end
    for i_true=1:n_groups
        for i_pred=1:n_groups
            c_value = confusion_normalized(i_true,i_pred);
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
            if n_batch > 1
                text(ax,i_pred-0.5,i_true-0.5,...
                    sprintf('%.2f%%\n\x00B1\n%.2f%%',...
                    confusion_mean(i_true,i_pred)*100,...
                    confusion_std(i_true,i_pred)*100),...
                    'HorizontalAlignment','center',...
                    'FontSize',18,...
                    'Color',t);
            else
                text(ax,i_pred-0.5,i_true-0.5,...
                    sprintf('%.2f%%',...
                    confusion_mean(i_true,i_pred)*100),...
                    'HorizontalAlignment','center',...
                    'FontSize',18,...
                    'Color',t);
            end
        end
    end
    disableDefaultInteractivity(ax)
end