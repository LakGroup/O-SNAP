function plot_NuCLiPSE(T, groups, save_dir, options)
arguments
    T table
    groups cell
    save_dir string
    options.plot_PCA logical = true
    options.plot_histograms logical = false
    options.cluster_observations logical = true
    options.cluster_features logical = false
end
% var_sel_methods = ["fscc hi2","fscmrmr","relieff"];
var_sel_methods = ["fscmrmr"];
idx = any(reshape(cell2mat(arrayfun(@(x) strcmpi(T.group,x), groups,'uni',0)),[],length(groups)),2);
T = T(idx,:);
[groups, ~, group_idx] = unique(T.group);
n_vars_sel = length(groups)+1;
for m=1:length(var_sel_methods)
    % set up
    var_sel_method = var_sel_methods(m);
    figure('Position', [10 10 1210 910]); 
    tiledlayout('flow')
    sgtitle(var_sel_method);
    %% VARIABLE SELECTED PCA
    var_selected = select_features(T,n_vars_sel,var_sel_method);
    nexttile
    plot_PCA(T(:,[{'group'} var_selected]));
    % SHUFFLED LABELS
    groups_shuffled = T.group(randperm(length(T.group)));
    T_shuffled = T(:,2:end);
    T_shuffled.group = groups_shuffled;
    var_selected_shuffled = select_features(T_shuffled, n_vars_sel,var_sel_method);
    title("Variable Score - shuffled labels")
    nexttile
    plot_PCA(T_shuffled(:,[{'group'} var_selected_shuffled]));
    title("PCA - shuffled labels")
    if ~options.plot_PCA
        close()
    else
        savefig(fullfile(save_dir, "PCA_"+var_sel_method+"_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "PCA_"+var_sel_method+"_"+join(groups,'_')+".png"))
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
                min_bin = min(T{:, var_selected{k}});
                max_bin = max(T{:, var_selected{k}});
                histogram(T{group_idx==g, var_selected{k}},...
                    'BinEdges', min_bin:(max_bin-min_bin)/15:max_bin,...
                    'DisplayStyle','stairs','Normalization','probability');
                hold on
            end
            xlabel(replace(var_selected{k},"_"," "))
            ylim([0 0.75])
        end
        sgtitle(var_sel_method)
        savefig(fullfile(save_dir, "hist_"+var_sel_method+"_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "hist_"+var_sel_method+"_"+join(groups,'_')+".png"));
        save(fullfile(save_dir, "vars_selected_" +var_sel_method+ ".mat"),"var_selected")
    end
end
if options.cluster_observations
    try
        %% unsupervised hierarchical clustering (observations)
        cluster_observations_hierarchical(T, var_selected);
        savefig(fullfile(save_dir, "clustering_obs_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "clustering_obs_"+join(groups,'_')+".png"));
    catch ME
        disp(getReport(ME))
    end
end
if options.cluster_features
    try
        %% unsupervised hierarchical clustering (features)
        cluster_features_hierarchical(T, var_selected);
        savefig(fullfile(save_dir, "clustering_feat_"+join(groups,'_')+".fig"));
        saveas(gcf,fullfile(save_dir, "clustering_obs_"+join(groups,'_')+".png"));
    catch ME
        disp(getReport(ME))
    end
end
end

