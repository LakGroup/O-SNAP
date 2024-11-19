function plot_SNAP(T, groups, save_dir, options)
arguments
    T table
    groups cell
    save_dir string
    options.plot_PCA logical = true
    options.plot_shuffled logical = true
    options.plot_histograms logical = true
    options.cluster_observations logical = true
    options.cluster_features logical = true
end
% var_sel_methods = ["fscchi2","fscmrmr","relieff"];
vars = T(:,vartype('numeric')).Properties.VariableNames;
var_sel_methods = ["fscmrmr"];
n_vars_sel = length(groups)+1;
for m=1:length(var_sel_methods)
    % set up
    var_sel_method = var_sel_methods(m);
    %% VARIABLE SELECTED PCA
    [vars_ranked, importance_score] = select_features_SNAP(T,"method",var_sel_method);
    vars_sel = vars_ranked(1:n_vars_sel);
    if options.plot_PCA
        figure('Position', [10 10 1210 910]); 
        plot_var_selected(vars_ranked,importance_score,"n_select",n_vars_sel)
        savefig(fullfile(save_dir, "features_importance_"+var_sel_method+"_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "features_importance_"+var_sel_method+"_"+join(groups,'_')+".png"))
        figure('Position', [10 10 1210 910]);
        plot_PCA(T(:,[{'group','biological_replicate'} vars_sel]));
        savefig(fullfile(save_dir, "PCA_"+var_sel_method+"_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "PCA_"+var_sel_method+"_"+join(groups,'_')+".png"))
    end
    %% SHUFFLED LABELS
    if options.plot_PCA && options.plot_shuffled
        figure('Position', [10 10 1210 910]);
        T_shuffled = T(:,2:end);
        T_shuffled.group = T.group(randperm(length(T.group)));
        T_shuffled = T_shuffled(:,[end 1:end-1]);
        [vars_ranked_shuffled, importance_score] = select_features_SNAP(T_shuffled,"method",var_sel_method);
        plot_PCA(T_shuffled(:,[{'group'} vars_ranked_shuffled(1:n_vars_sel)]));
        title("PCA - shuffled labels - " + var_sel_method)
        savefig(fullfile(save_dir, "PCA_"+var_sel_method+"_"+join(groups,'_')+"_shuffled.fig"));
        saveas(gcf,fullfile(save_dir, "PCA_"+var_sel_method+"_"+join(groups,'_')+"_shuffled.png"))
    end
    % hist plot selected features
    if options.plot_histograms
        figure('Position', [10 10 1210 910]);
        n_var_fac = factor(n_vars_sel);
        if(n_var_fac(end) > 4)
            n_var_fac = factor(n_vars_sel+1);
        end
        for k = 1:n_vars_sel
            subplot(n_var_fac(end), prod(n_var_fac(1:end-1)), k)
            for g=1:length(groups)
                min_bin = min(T{:, vars_sel{k}});
                max_bin = max(T{:, vars_sel{k}});
                histogram(T{T.group==groups(g), vars_sel{k}},...
                    'BinEdges', min_bin:(max_bin-min_bin)/15:max_bin,...
                    'DisplayStyle','bar','Normalization','probability');
                hold on
            end
            xlabel(replace(vars_sel{k},"_"," "))
            ylim([0 0.75])
            if k==1
                legentry=repelem({''},1,length(groups));
                color_order = get(gca,'ColorOrder');
                for g=1:length(groups)
                   scatter(nan,nan,20,'square','filled','MarkerFaceColor',color_order(g,:));
                   legentry{g} = groups{g};
                end
                leg = legend(legentry);
                leg.FontSize = 16;
                leg.IconColumnWidth = 20;
            end
        end
        sgtitle(var_sel_method)
        savefig(fullfile(save_dir, "hist_"+var_sel_method+"_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "hist_"+var_sel_method+"_"+join(groups,'_')+".png"));
        save(fullfile(save_dir, "vars_selected_" +var_sel_method+ ".mat"),"vars_sel")
    end
end
if options.cluster_observations
    try
        %% unsupervised hierarchical clustering (observations)
        cluster_observations_hierarchical(T, vars_sel);
        savefig(fullfile(save_dir, "clustering_obs_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "clustering_obs_"+join(groups,'_')+".png"));
    catch ME
        disp(getReport(ME))
    end
end
if options.cluster_features
    try
        %% unsupervised hierarchical clustering (features)
        cluster_features_hierarchical(T, vars_sel);
        savefig(fullfile(save_dir, "clustering_feat_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "clustering_feat_"+join(groups,'_')+".png"));
    catch ME
        disp(getReport(ME))
    end
end
end

