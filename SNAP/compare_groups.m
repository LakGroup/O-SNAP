function [S,v_fig] = compare_groups(T,group_1,group_2,options)
    arguments
        T table;
        group_1 string;
        group_2 string;
        options.alpha double = 0.05;
        options.fold_change_threshold double = 2;
        options.plot logical = true;
        options.plot_labels logical = true;
        options.save_path string = "";
    end
    [T_num, group_idx] = prepare_voronoi_table_data(T,"normalize",0,"numeric_only",1);
    mean_1 = zeros(1,size(T_num,2));
    mean_2 = zeros(1,size(T_num,2));
    p = zeros(1,size(T_num,2));
    for i=1:size(T_num,2)
        data_i = T_num{:,i};
        data_group_1 = rmoutliers(data_i(strcmpi(group_idx,group_1)));
        data_group_2 = rmoutliers(data_i(strcmpi(group_idx,group_2)));
        mean_1(i) = mean(data_group_1,'omitmissing');
        mean_2(i) = mean(data_group_2,'omitmissing');
        [~,p(i)] = ttest2(data_group_1,data_group_2);
    end
    p_adj = pval_adjust(p, 'BH'); % IMPORTANT: adjust for multiple hypothesis testing
    %% Filter samples with no differential expression
    % note: features where directionality is represented in the sign of the
    %       value (ex: skewness) will not be preserved for the volcano plot
    %       analysis
    ratios = abs(mean_2./mean_1); 
    fold_change = ratios;
    fold_change(mean_2<mean_1) = -1./ratios(mean_2<mean_1);
    S = table(string(T(:,4:end).Properties.VariableNames'),...
        mean_1',...
        mean_2',...
        fold_change',...
        log2(ratios'),...
        p',...
        p_adj',...
        'VariableNames',...
        ["feature",...
        "mean_"+group_1,...
        "mean_"+group_2,...
        "fold_change",...
        "log2_fold_change",...
        "p_value",...
        "adj_p_value"...
        ]);
    if options.plot
        v_fig = volcano_plot(S,group_1,group_2,"alpha",options.alpha,"fold_change_threshold",options.fold_change_threshold,"plot_labels",options.plot_labels);
    else
        v_fig = [];
    end
    if ~strcmp(options.save_path, "")
        writetable(sortrows(S,"adj_p_value"),options.save_path);
    end
end

%% generate volcano plot
function v_fig = volcano_plot(S,group_1,group_2,options)
    arguments
        S table;
        group_1 string;
        group_2 string;
        options.alpha double = 0.05;
        options.fold_change_threshold double = 2;
        options.plot_labels logical = true;
    end
    log2_fold_change_threshold = log2(options.fold_change_threshold);
    idx_sig = all([abs(S.log2_fold_change)>log2_fold_change_threshold S.adj_p_value < options.alpha],2);
    idx_sig_up = all([S.log2_fold_change>0 idx_sig],2);
    idx_sig_down = all([S.log2_fold_change<0 idx_sig],2);
    v_fig_size = 20*ones(size(S,1),1);
    v_fig_size(idx_sig) = 35*ones(sum(idx_sig),1);
    v_fig_color = repmat([0.5 0.5 .5],size(S,1),1);
    v_fig_color(idx_sig_up,:) = repmat([1 0 0],sum(idx_sig_up),1);
    v_fig_color(idx_sig_down,:) = repmat([0 0 1],sum(idx_sig_down),1);
    v_fig = figure();
    scatter(S.log2_fold_change,-log10(S.adj_p_value),v_fig_size,v_fig_color,'filled','MarkerEdgeColor','k');
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
        x = S{idx_sig,"log2_fold_change"}+0.025;
        y = -log10(S{idx_sig,"adj_p_value"});
        feature_name = replace(S{idx_sig,"feature"},"_"," ");
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
end
