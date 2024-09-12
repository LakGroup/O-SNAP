function var_selected = select_features(T,n_var_sel,method)
%% set up
groups = T.group;
T_num =  prepare_voronoi_table_data(T,"numeric_only",true); % get only numeric columns
vars = T_num.Properties.VariableNames;
X = table2array(T_num);

%% variable selection
if strcmpi(method,"relieff")
    [var_sel_idx,var_sel_scores] = relieff(X,cellstr(groups),10,'method','classification');
else
    [var_sel_idx,var_sel_scores] = feval(lower(method),X,groups);
end
var_selected = vars(var_sel_idx(1:n_var_sel));
importance_score = var_sel_scores(var_sel_idx(1:12));

%% plot
nexttile
b = bar(importance_score/sum(importance_score),'FaceColor','flat');
b.CData(1:n_var_sel,:) = repmat([0.8500 0.3250 0.0980],n_var_sel,1);
ylabel({'Predictor importance score','(probability)'})
ylim([0 1])
set(gca,'XTick',1:12, 'XTickLabel',replace(vars(var_sel_idx(1:12)),"_"," "))
title('Variable Score')
end
