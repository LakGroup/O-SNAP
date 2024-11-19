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
    venn_data = interpret_features_venn_3(T, save_dir, "plot", options.plot,"alpha",options.alpha,"fold_change_threshold",options.fold_change_threshold);
% elseif n_groups == 4
    % venn_data = interpret_features_venn_4(T);
else
    disp("   Not enough groups (N>2)")
    return
end
save(fullfile(save_dir,"features_venn.mat"),"venn_data");
end

function venn_data = interpret_features_venn_3(T,save_dir, options)
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
        venn_data{i}.increase_features = results{i}.feature_table{idx_sig_up,"feature"};
        venn_data{i}.decrease_features = results{i}.feature_table{idx_sig_down,"feature"};
    else
        venn_data{i}.increase_features = [];
        venn_data{i}.decrease_features = [];
    end
        venn_data{i}.n_increase_features = size(venn_data{i}.increase_features,1);
        venn_data{i}.n_decrease_features = size(venn_data{i}.decrease_features,1);
end
%% get shared DEFs between different group comparisons
group_pairs_pairs = nchoosek(group_pairs_string,2);
group_pairs_pairs_string = join(group_pairs_pairs, " VS ");
n_group_pairs_pairs = size(group_pairs_pairs,1);
for i=1:n_group_pairs_pairs
    group_pair_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_pair_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{n_group_pairs+i}.groups = group_pairs_pairs_string(i);
    venn_data{n_group_pairs+i}.increase_features = intersect(...
        venn_data{group_pair_idx_1}.increase_features,...
        venn_data{group_pair_idx_2}.increase_features);
    venn_data{n_group_pairs+i}.decrease_features = intersect(...
        venn_data{group_pair_idx_1}.decrease_features,...
        venn_data{group_pair_idx_2}.decrease_features);
    venn_data{n_group_pairs+i}.n_increase_features = size(venn_data{n_group_pairs+i}.increase_features,1);
    venn_data{n_group_pairs+i}.n_decrease_features = size(venn_data{n_group_pairs+i}.decrease_features,1);
end
%% get shared DEFs between all groups
venn_data{end}.groups = join(groups," vs ");
venn_data{end}.increase_features = feature_intersection_3(venn_data{1}.increase_features,...
    venn_data{2}.increase_features,...
    venn_data{3}.increase_features);
venn_data{end}.decrease_features = feature_intersection_3(venn_data{1}.decrease_features,...
    venn_data{2}.decrease_features,...
    venn_data{3}.decrease_features);
venn_data{end}.n_increase_features = length(venn_data{end}.increase_features);
venn_data{end}.n_decrease_features = length(venn_data{end}.decrease_features);
%% correct counts
for i=1:(length(venn_data)-1)
    venn_data{i}.n_increase_features = size(venn_data{i}.increase_features,1) - venn_data{end}.n_increase_features;
    venn_data{i}.n_decrease_features = size(venn_data{i}.decrease_features,1) - venn_data{end}.n_decrease_features;
end
for i=1:n_group_pairs
    group_pair_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_pair_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{group_pair_idx_1}.n_increase_features =  venn_data{group_pair_idx_1}.n_increase_features - venn_data{n_group_pairs+i}.n_increase_features;
    venn_data{group_pair_idx_2}.n_increase_features = venn_data{group_pair_idx_2}.n_increase_features - venn_data{n_group_pairs+i}.n_increase_features;
    venn_data{group_pair_idx_1}.n_decrease_features = venn_data{group_pair_idx_1}.n_decrease_features - venn_data{n_group_pairs+i}.n_decrease_features;
    venn_data{group_pair_idx_2}.n_decrease_features = venn_data{group_pair_idx_2}.n_decrease_features - venn_data{n_group_pairs+i}.n_decrease_features;
end
%% plot result
if options.plot
    venn_labels_increase = cellfun(@(x) x.n_increase_features,venn_data,'uni',1);
    venn(n_groups,'sets',group_pairs_string,'labels',venn_labels_increase,'alpha',0.7,'colors',autumn(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE INCREASE",'Position', [1.5 -0.3, 0])
    savefig(fullfile(save_dir, "features_venn_increase.fig"));
    saveas(gcf,fullfile(save_dir, "features_venn_increase.png"));
    venn_labels_decrease = cellfun(@(x) x.n_decrease_features,venn_data,'uni',1);
    venn(n_groups,'sets',group_pairs_string,'labels',venn_labels_decrease,'alpha',0.7,'colors',winter(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE DECREASE",'Position', [1.5 -0.3, 0])
    savefig(fullfile(save_dir, "features_venn_decrease.fig"));
    saveas(gcf,fullfile(save_dir, "features_venn_decrease.png"));
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
