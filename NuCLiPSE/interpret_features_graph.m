function interpret_features_graph(T, save_dir, options)
arguments
    T table
    save_dir string
    options.thresh double = 0.8
end

% get table information
T_norm = T(:,4:end);
vars = replace(T_norm.Properties.VariableNames,'_',' ');
T_norm = table2array(normalize(T_norm));
% create correlation matrix
T_corr = abs(corrcoef(T_norm));
nan_idx = all(isnan(T_corr));
% filter out NaNs
T_norm = T(:,~nan_idx);
T_corr = T_corr(~nan_idx,~nan_idx);
vars = vars(~nan_idx);
% create correlation matrix
n_features = length(vars);
T_corr(T_corr < options.thresh) = 0;
%% TODO: change so that is compatible with only one group
if length(unique(T.group)) < 2
    return
end
[feature_rank,score] = fscmrmr(T_norm,T.group);

if all(score == 0)
    return
end
% create graph object
G = graph(T_corr,vars,'omitselfloops');
% identify notable features
n_features_selected = length(unique(T.group))+1;
feature_idx = feature_rank(1:n_features_selected);
% change node sizes, colors, and font weight to correspond to selected vs not selected features
node_size = 5*score/max(score) + 2;
node_size(feature_idx) = 2*node_size(feature_idx);
node_colormap = lines(2);
node_color = repmat(node_colormap(1,:),n_features,1);
node_color(feature_idx,:) = repmat(node_colormap(2,:),n_features_selected,1);
node_font_weight = repmat("Normal",n_features,1);
node_font_weight(feature_idx) = "Bold";
G.Nodes{feature_idx,:} = upper(G.Nodes{feature_idx,:});
node_font_size = repmat(8,n_features,1);
node_font_size(feature_idx) = repmat(12,n_features_selected,1);
% % select only for graphs with connections
% idx_isolated = find(degree(G) == 0);
% idx_connected = find(degree(G) ~= 0);
% G = rmnode(G,idx_isolated);
% node_color = node_color(idx_connected,:);
% node_size = node_size(idx_connected);
% node_font_weight = node_font_weight(idx_connected);
% node_font_size = node_font_size(idx_connected);
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