function [feature_compare_data,feature_compare_data_filtered] = compare_OSNAP_group_pair(feature_data,group_1,group_2,options)
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
    [T_num, group_idx] = preprocess_OSNAP_feature_data(feature_data,"normalize",false,"numeric_only",1);
    % note that all information stored in signs is erased
    % T_num = abs(T_num);
    mean_1 = zeros(size(T_num,2),1);
    mean_2 = zeros(size(T_num,2),1);
    p = zeros(size(T_num,2),1);
    fold_change = zeros(size(T_num,2),1);
    idx_group_1 = strcmpi(group_idx,group_1);
    idx_group_2 = strcmpi(group_idx,group_2);
    for i=1:size(T_num,2)
        data_i = T_num{:,i};
        if options.remove_outliers
            data_group_1 = rmoutliers(data_i(idx_group_1));
            data_group_2 = rmoutliers(data_i(idx_group_2));
        else
            data_group_1 = data_i(idx_group_1);
            data_group_2 = data_i(idx_group_2);
        end
        mean_1(i) = mean(data_group_1,'omitmissing');
        mean_2(i) = mean(data_group_2,'omitmissing');
        mean_signs = sign([mean_1(i) mean_2(i)]);
        if all(mean_signs==1) % both means have same sign (FC>1 when mean2>mean1)
            fold_change(i) = mean_2(i)/mean_1(i);
            [~,p(i)] = ttest2(data_group_1,data_group_2,'vartype','unequal');
        elseif all(mean_signs==-1) % both means negative (FC>1 when mean2<mean1)
            fold_change(i) = mean_1(i)/mean_2(i);
            [~,p(i)] = ttest2(data_group_1,data_group_2,'vartype','unequal');
        else % If one mean is positive, other is negative. should only happen w/ skewness features. 
            % In this case, just compare mean2 to mean 1
            % note: features where directionality is represented in the sign of the
            % value (ex: skewness) will not be preserved for the volcano plot analysis
            fold_change(i) = abs(mean_2(i))/abs(mean_1(i));
            [~,p(i)] = ttest2(abs(data_group_1),abs(data_group_2),'vartype','unequal');
        end
    end
    p_adj = pval_adjust(p, options.p_adj_method); % IMPORTANT: adjusts for multiple hypothesis testing
    %% Filter samples with no differential expression
    feature_compare_data = table(string(feature_data(:,4:end).Properties.VariableNames'),...
        mean_1,...
        mean_2,...
        fold_change,...
        log2(fold_change),...
        p,...
        p_adj,...
        'VariableNames', ...
        ["feature",...
        "mean_"+group_1,...
        "mean_"+group_2,...
        "fold_change",...
        "log2_fold_change",...
        "p_value",...
        "adj_p_value"...
        ]);
        feature_compare_data_filtered =  sortrows(feature_compare_data,"log2_fold_change");
        feature_compare_data_filtered = feature_compare_data_filtered(all([...
            abs(feature_compare_data_filtered.log2_fold_change)>log2(options.fold_change_threshold),...
            feature_compare_data_filtered.adj_p_value<options.alpha],2),:);
    if options.save_path ~= ""
        plot_OSNAP_volcano(feature_compare_data,group_1,group_2,...
            "alpha",options.alpha,...
            "fold_change_threshold",options.fold_change_threshold,...
            ..."plot_labels",options.plot_labels,...
            "save_path",options.save_path);
        writetable(sortrows(feature_compare_data,"log2_fold_change"),join([options.save_path,group_1,group_2],"_")+".csv");
        writetable(feature_compare_data_filtered,join([options.save_path,group_1,group_2,"filtered.csv"],"_"));
    end
end