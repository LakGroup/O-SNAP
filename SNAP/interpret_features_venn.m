function interpret_features_venn(T,save_dir,options)
arguments
    T table
    save_dir string
    options.plot = true;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
end
%% find number of sections for venn diagram
groups = unique(T.group);
n_groups = length(groups);
if n_groups == 3
    interpret_features_venn_3(T, save_dir, "plot", options.plot,"alpha",options.alpha,"fold_change_threshold",options.fold_change_threshold)
% elseif n_groups == 4
    % interpret_features_venn_4(T)
else
    disp("   Not enough groups (N>2)")
end
end

function interpret_features_venn_3(T,save_dir, options)
arguments
    T table
    save_dir string
    options.plot = false;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
end
%% data preparation
groups = unique(T.group);
n_groups = length(groups);
group_pairs = nchoosek(groups,2);
group_pairs_string = join(group_pairs, " vs ");
n_group_pairs = size(group_pairs,1);
results = cell(1,n_group_pairs);
venn_data = cell(1,7);
log2_fold_change_threshold = log2(options.fold_change_threshold);

%% get comparison between two groups
for i=1:n_group_pairs
    results{i}.groups = group_pairs(i,:);
    results{i}.feature_table = compare_groups(T,group_pairs(i,1),group_pairs(i,2),"plot",false);
    venn_data{i}.groups = group_pairs_string(i);
    if ~isempty(results{i}.feature_table)
        idx_sig = all([abs(results{i}.feature_table.log2_fold_change)>log2_fold_change_threshold results{i}.feature_table.adj_p_value < options.alpha],2);
        idx_sig_up = all([results{i}.feature_table.log2_fold_change>0 idx_sig],2);
        idx_sig_down = all([results{i}.feature_table.log2_fold_change<0 idx_sig],2);
        venn_data{i}.upregulated_features = results{i}.feature_table{idx_sig_up,"feature"};
        venn_data{i}.downregulated_features = results{i}.feature_table{idx_sig_down,"feature"};
    else
        venn_data{i}.upregulated_features = [];
        venn_data{i}.downregulated_features = [];
    end
        venn_data{i}.n_upregulated_features = size(venn_data{i}.upregulated_features,1);
        venn_data{i}.n_downregulated_features = size(venn_data{i}.downregulated_features,1);
end
%% get shared DEFs between different group comparisons
group_pairs_pairs = nchoosek(group_pairs_string,2);
group_pairs_pairs_string = join(group_pairs_pairs, " VS ");
n_group_pairs_pairs = size(group_pairs_pairs,1);
for i=1:n_group_pairs_pairs
    group_pair_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_pair_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{n_group_pairs+i}.groups = group_pairs_pairs_string(i);
    venn_data{n_group_pairs+i}.upregulated_features = intersect(...
        venn_data{group_pair_idx_1}.upregulated_features,...
        venn_data{group_pair_idx_2}.upregulated_features);
    venn_data{n_group_pairs+i}.downregulated_features = intersect(...
        venn_data{group_pair_idx_1}.downregulated_features,...
        venn_data{group_pair_idx_2}.downregulated_features);
    venn_data{n_group_pairs+i}.n_upregulated_features = size(venn_data{n_group_pairs+i}.upregulated_features,1);
    venn_data{n_group_pairs+i}.n_downregulated_features = size(venn_data{n_group_pairs+i}.downregulated_features,1);
end
%% get shared DEFs between all groups
venn_data{end}.groups = join(groups," vs ");
venn_data{end}.upregulated_features = feature_intersection_3(venn_data{1}.upregulated_features,...
    venn_data{2}.upregulated_features,...
    venn_data{3}.upregulated_features);
venn_data{end}.downregulated_features = feature_intersection_3(venn_data{1}.downregulated_features,...
    venn_data{2}.downregulated_features,...
    venn_data{3}.downregulated_features);
venn_data{end}.n_upregulated_features = length(venn_data{end}.upregulated_features);
venn_data{end}.n_downregulated_features = length(venn_data{end}.downregulated_features);
%% correct counts
for i=1:(length(venn_data)-1)
    venn_data{i}.n_upregulated_features = size(venn_data{i}.upregulated_features,1) - venn_data{end}.n_upregulated_features;
    venn_data{i}.n_downregulated_features = size(venn_data{i}.downregulated_features,1) - venn_data{end}.n_downregulated_features;
end
for i=1:n_group_pairs
    group_pair_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_pair_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{group_pair_idx_1}.n_upregulated_features =  venn_data{group_pair_idx_1}.n_upregulated_features - venn_data{n_group_pairs+i}.n_upregulated_features;
    venn_data{group_pair_idx_2}.n_upregulated_features = venn_data{group_pair_idx_2}.n_upregulated_features - venn_data{n_group_pairs+i}.n_upregulated_features;
    venn_data{group_pair_idx_1}.n_downregulated_features = venn_data{group_pair_idx_1}.n_downregulated_features - venn_data{n_group_pairs+i}.n_downregulated_features;
    venn_data{group_pair_idx_2}.n_downregulated_features = venn_data{group_pair_idx_2}.n_downregulated_features - venn_data{n_group_pairs+i}.n_downregulated_features;
end
%% plot result
if options.plot
    venn_labels_upregulated = cellfun(@(x) x.n_upregulated_features,venn_data,'uni',1);
    venn(n_groups,'sets',group_pairs_string,'labels',venn_labels_upregulated,'alpha',0.7,'colors',autumn(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("UP-REGULATED",'Position', [1.5 -0.3, 0])
    savefig(fullfile(save_dir, "features_venn_upregulated.fig"));
    saveas(gcf,fullfile(save_dir, "features_venn_upregulated.png"));
    venn_labels_downregulated = cellfun(@(x) x.n_downregulated_features,venn_data,'uni',1);
    venn(n_groups,'sets',group_pairs_string,'labels',venn_labels_downregulated,'alpha',0.7,'colors',winter(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("DOWN-REGULATED",'Position', [1.5 -0.3, 0])
    savefig(fullfile(save_dir, "features_venn_downregulated.fig"));
    saveas(gcf,fullfile(save_dir, "features_venn_downregulated.png"));
end
end

function v = feature_intersection_2(v1, v2)
    if (isempty(v1)) || (isempty(v2))
        v = [];
        return
    else
        v = intersect(v1, v2);
    end
    if(isempty(v))
        v = [];
    end
end

function v = feature_intersection_3(v1, v2, v3)
    v_23 = feature_intersection_2(v2, v3);
    if ~isempty(v_23)
        v = feature_intersection_2(v1, v_23);
    else
        v = [];
        return
    end
    if isempty(v)
        v = [];
    end
end
