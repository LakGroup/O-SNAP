function feature_comparison_pair = compare_group_pair(feature_data,group_1,group_2,options)
    arguments
        feature_data table;
        group_1 string;
        group_2 string;
        options.alpha double = 0.05;
        options.p_adj_method char = 'BH'
        options.fold_change_threshold double = 2;
        options.plot_labels logical = true;
        options.save_path string = "";
        options.remove_outliers logical = false;
    end
    [T_num, group_idx] = prepare_voronoi_table_data(feature_data,"normalize",0,"numeric_only",1);
    mean_1 = zeros(1,size(T_num,2));
    mean_2 = zeros(1,size(T_num,2));
    p = zeros(1,size(T_num,2));
    for i=1:size(T_num,2)
        data_i = T_num{:,i};
        if options.remove_outliers
            data_group_1 = rmoutliers(data_i(strcmpi(group_idx,group_1)));
            data_group_2 = rmoutliers(data_i(strcmpi(group_idx,group_2)));
        else
            data_group_1 = data_i(strcmpi(group_idx,group_1));
            data_group_2 = data_i(strcmpi(group_idx,group_2));
        end
        mean_1(i) = mean(data_group_1,'omitmissing');
        mean_2(i) = mean(data_group_2,'omitmissing');
        [~,p(i)] = ttest2(data_group_1,data_group_2);
    end
    p_adj = pval_adjust(p, options.p_adj_method); % IMPORTANT: adjusts for multiple hypothesis testing
    %% Filter samples with no differential expression
    % note: features where directionality is represented in the sign of the
    %       value (ex: skewness) will not be preserved for the volcano plot
    %       analysis
    ratios = abs(mean_2./mean_1); 
    fold_change = ratios;
    fold_change(mean_2<mean_1) = -1./ratios(mean_2<mean_1);
    feature_comparison_pair = table(string(feature_data(:,4:end).Properties.VariableNames'),...
        mean_1',...
        mean_2',...
        fold_change',...
        log2(ratios'),...
        p',...
        p_adj',...
        'VariableNames', ...
        ["feature",...
        "mean_"+group_1,...
        "mean_"+group_2,...
        "fold_change",...
        "log2_fold_change",...
        "p_value",...
        "adj_p_value"...
        ]);
    if options.save_path ~= ""
        volcano_plot_SNAP(feature_comparison_pair,group_1,group_2,...
            "alpha",options.alpha,...
            "fold_change_threshold",options.fold_change_threshold,...
            "plot_labels",options.plot_labels,...
            "save_path",options.save_path);
    end
end

%% generate volcano plot
function volcano_plot_SNAP(feature_comparison_pair,group_1,group_2,options)
    arguments
        feature_comparison_pair table;
        group_1 string;
        group_2 string;
        options.alpha double = 0.05;
        options.fold_change_threshold double = 2;
        options.plot_labels logical = true;
        options.save_path string = "";
    end
    log2_fold_change_threshold = log2(options.fold_change_threshold);
    idx_sig = all([abs(feature_comparison_pair.log2_fold_change)>log2_fold_change_threshold feature_comparison_pair.adj_p_value < options.alpha],2);
    idx_sig_up = all([feature_comparison_pair.log2_fold_change>0 idx_sig],2);
    idx_sig_down = all([feature_comparison_pair.log2_fold_change<0 idx_sig],2);
    v_fig_size = 20*ones(size(feature_comparison_pair,1),1);
    v_fig_size(idx_sig) = 35*ones(sum(idx_sig),1);
    v_fig_color = repmat([0.5 0.5 .5],size(feature_comparison_pair,1),1);
    v_fig_color(idx_sig_up,:) = repmat([1 0 0],sum(idx_sig_up),1);
    v_fig_color(idx_sig_down,:) = repmat([0 0 1],sum(idx_sig_down),1);
    v_fig = figure();
    scatter(feature_comparison_pair.log2_fold_change,-log10(feature_comparison_pair.adj_p_value),v_fig_size,v_fig_color,'filled','MarkerEdgeColor','k');
    hold on
    plot([log2_fold_change_threshold log2_fold_change_threshold],ylim(),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    hold on
    plot(-[log2_fold_change_threshold log2_fold_change_threshold],ylim(),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    hold on
    plot(xlim(), [-log10(options.alpha) -log10(options.alpha)],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    ylabel("-log_{10}(adj. p-value)")
    xlabel("log_{2}(Fold Change)")
    title(replace("\color{blue}" + group_1 + " \color{black}vs " + "\color{red}" + group_2,"_","-"));
    if options.plot_labels
        hold on
        x = feature_comparison_pair{idx_sig,"log2_fold_change"}+0.025;
        y = -log10(feature_comparison_pair{idx_sig,"adj_p_value"});
        feature_name = replace(feature_comparison_pair{idx_sig,"feature"},"_"," ");
        for f=1:size(feature_name,1)
            if x(f) < 0
                t_color = 'b';
            elseif x(f) > 0
                t_color = 'r';
            else
                t_color = 'k';
            end
            ht = text(x(f),y(f),'    ','HorizontalAlignment','center',...
                'tag',num2str(f),'Color',t_color,'FontSize',10,'Interpreter','none');
            pointerBehavior.enterFcn = ...
                @(v_fig, cpp)set(findobj(v_fig,'tag',num2str(f)), ...
                'string', {feature_name(f),""});
            pointerBehavior.traverseFcn = [];
            pointerBehavior.exitFcn = ...
                @(v_fig, cpp)set(findobj(v_fig, 'tag',num2str(f)), ...
                'string', '    ');
            iptSetPointerBehavior(ht, pointerBehavior);
            iptPointerManager(v_fig, 'enable');
            hold on
        end
    end
    if options.save_path ~= ""
        save_path = sprintf("%s_%s_%s",options.save_path,group_1,group_2);
        savefig(save_path+".fig");
        saveas(gcf,save_path+".png");
        writetable(sortrows(feature_comparison_pair,"adj_p_value"),save_path);
    end
end
