% get table information
T_norm = prepare_voronoi_table_data(T,"numeric_only",true);
vars = replace(T_norm.Properties.VariableNames,'_',' ');
% create correlation matrix
T_corr = abs(corrcoef(table2array(T_norm)));
% remove nans
nan_idx = all(isnan(T_corr));
T_corr = T_corr(~nan_idx,~nan_idx);
vars = vars(~nan_idx);
% create correlation matrix
n_features = length(vars);
% T_corr(T_corr < options.thresh) = 0;
% identify notable features
[feature_rank,score] = fscmrmr(T_norm,T.group);
n_features_selected = length(unique(T.group))+1;
feature_idx = feature_rank(1:n_features_selected);
if all(score == 0)
    return
end

vars(feature_idx) = upper(vars(feature_idx));
for i=1:numel(vars)
    if numel(vars{i}) > 20
        vars{i} = [vars{i}(1:8) '...' vars{i}(end-6:end)];
    end
end

figure('Position',[10 10 1010 910]);
h = heatmap(T_corr,'MissingDataColor','w');
h.XDisplayLabels = vars;
h.YDisplayLabels = vars;
<<<<<<< HEAD:NuCLiPSE/interpret_features_heatmap.m
h.FontSize = 6;
=======
h.FontSize = 6;
>>>>>>> experimental:SNAP/interpret_features_heatmap.m
