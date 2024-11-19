function interpret_features_graph(T, save_dir, options)
arguments
    T table
    save_dir string
    options.thresh double = 0.8
    options.connected_only logical = true;
end
if length(unique(T.group)) < 2
    return
end
% get table information
[T_norm, group_idx] = prepare_voronoi_table_data(T,"numeric_only",true);
vars = replace(T_norm.Properties.VariableNames,'_',' ');
% create correlation matrix
T_corr = abs(corrcoef(table2array(T_norm)));
% % remove nans
% nan_idx = all(isnan(T_corr));
% T_corr = T_corr(~nan_idx,~nan_idx);
% vars = vars(~nan_idx);
% create correlation matrix
n_features = length(vars);
T_corr(T_corr < options.thresh) = 0;
% identify notable features
[feature_rank,score] = fscmrmr(T_norm,group_idx);
n_features_selected = length(unique(T.group))+1;
feature_idx = feature_rank(1:n_features_selected);
if all(score == 0)
    return
end
% create graph object
G = graph(T_corr,vars,'omitselfloops');
% change node sizes, colors, and font weight to correspond to selected vs not selected features
node_size = 5*sqrt(score/max(score)) + 2;
node_size(feature_idx) = node_size(feature_idx);
node_colormap = lines(2);
node_color = repmat(node_colormap(1,:),n_features,1);
node_color(feature_idx,:) = repmat(node_colormap(2,:),n_features_selected,1);
node_font_weight = repmat("Normal",n_features,1);
node_font_weight(feature_idx) = "Bold";
G.Nodes{feature_idx,:} = upper(G.Nodes{feature_idx,:});
node_font_size = repmat(4,n_features,1);
node_font_size(feature_idx) = repmat(6,n_features_selected,1);
% select only for graphs with connections
if options.connected_only
    idx_isolated = find(degree(G) == 0);
    idx_connected = find(degree(G) ~= 0);
    G = rmnode(G,idx_isolated);
    node_color = node_color(idx_connected,:);
    node_size = node_size(idx_connected);
    node_font_weight = node_font_weight(idx_connected);
    node_font_size = node_font_size(idx_connected);
end
figure();
plot(G,...
    "NodeColor",node_color,...
    "MarkerSize",node_size,...
    "NodeFontWeight",node_font_weight, ...
    "NodeFontSize",node_font_size,...
    "NodeLabelMode","auto")
    %,'EdgeLabel',compose("%.2f",G.Edges.Weight))
title(join(unique(T.group), ", "))
savefig(fullfile(save_dir, "features_graph.fig"));
saveas(gcf,fullfile(save_dir, "features_graph.png"));
end
