function [var_ranked, importance_score] = select_features_SNAP(T,options)
arguments
    T table
    options.method string = "fscmrmr"
    options.plot logical = false
    options.save_dir string = ""
end
%% set up
[T_num,groups] =  prepare_voronoi_table_data(T,"numeric_only",true,"normalize",true); % get only numeric columns
vars = T_num.Properties.VariableNames;
X = table2array(T_num);
%% variable selection
if strcmpi(options.method,"relieff")
    [var_sel_idx,var_sel_scores] = relieff(X,cellstr(groups),10,'method','classification');
else
    [var_sel_idx,var_sel_scores] = feval(lower(options.method),X,groups);
end
var_ranked = vars(var_sel_idx);
importance_score = var_sel_scores(var_sel_idx);
%% plot
if options.plot && options.save_dir ~= ""
    figure('Position', [10 10 1210 910]); 
    plot_var_selected(vars_ranked,importance_score,"n_select",n_vars_sel)
    savefig(fullfile(save_dir, "features_importance_"+var_sel_method+"_"+join(groups,'_')+".fig"));
    saveas(gcf,fullfile(save_dir, "features_importance_"+var_sel_method+"_"+join(groups,'_')+".png"))
end
end

function plot_var_selected(vars,importance_score,options)
    arguments
        vars cell
        importance_score double
        options.n_select double = 4
        options.n_show double = 12
    end
    %% plot
    nexttile
    importance_score = importance_score(1:options.n_show);
    b = bar(importance_score/sum(importance_score),'FaceColor','flat');
    b.CData(1:options.n_select,:) = repmat([0.8500 0.3250 0.0980],options.n_select,1);
    ylabel({'Predictor importance score','(probability)'})
    % ylim([0 1])
    set(gca,'XTick',1:options.n_show, 'XTickLabel',replace(vars(1:options.n_show),"_"," "))
    title('Variable Score')
end