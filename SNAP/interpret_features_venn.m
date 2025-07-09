function venn_data = interpret_features_venn(feature_comparisons,options)
arguments
    feature_comparisons cell
    options.plot = true;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.save_path string = "";
end
%% find number of sections for venn diagram
n_pairs = numel(feature_comparisons);
if n_pairs == 3
    venn_data = interpret_features_venn_3(feature_comparisons,"save_path",options.save_path, "plot", options.plot,"alpha",options.alpha,"fold_change_threshold",options.fold_change_threshold);
elseif n_pairs == 4
    venn_data = interpret_features_venn_4(feature_comparisons,"save_path",options.save_path, "plot", options.plot,"alpha",options.alpha,"fold_change_threshold",options.fold_change_threshold);
else
    disp("   Not enough groups (N<2)")
    venn_data = [];
    return
end
save(fullfile(options.save_path+".mat"),"venn_data");
end

function venn_data = interpret_features_venn_3(feature_comparisons,options)
arguments
    feature_comparisons cell
    options.plot = false;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.save_path string = "";
end
%% data preparation
group_pairs_string = cellfun(@(x) join(x.groups, " vs "),feature_comparisons,'uni',1);
n_group_pairs = numel(feature_comparisons);
venn_data = cell(1,7);
log2_fold_change_threshold = log2(options.fold_change_threshold);
%% get comparison between two groups - level 1 (A, B, C)
for i=1:n_group_pairs
    venn_data{i}.groups = feature_comparisons{i}.groups;
    if ~isempty(feature_comparisons{i}.feature_table)
        idx_sig = all([abs(feature_comparisons{i}.feature_table.log2_fold_change)>log2_fold_change_threshold feature_comparisons{i}.feature_table.adj_p_value < options.alpha],2);
        idx_sig_up = all([feature_comparisons{i}.feature_table.log2_fold_change>0 idx_sig],2);
        idx_sig_down = all([feature_comparisons{i}.feature_table.log2_fold_change<0 idx_sig],2);
        venn_data{i}.increase_features_group_1 = feature_comparisons{i}.feature_table{idx_sig_up,"feature"};
        venn_data{i}.increase_features_group_2 = feature_comparisons{i}.feature_table{idx_sig_down,"feature"};
    else
        venn_data{i}.increase_features_group_1 = [];
        venn_data{i}.increase_features_group_2 = [];
    end
        venn_data{i}.n_increase_features_group_1 = size(venn_data{i}.increase_features_group_1,1);
        venn_data{i}.n_increase_features_group_2 = size(venn_data{i}.increase_features_group_2,1);
end
%% get shared DEFs between different group comparisons - level 2 (AB, AC, BC)
group_pairs_pairs = nchoosek(group_pairs_string,2);
group_pairs_pairs_string = join(group_pairs_pairs, ", ");
n_group_pairs_pairs = size(group_pairs_pairs,1);
for i=1:n_group_pairs_pairs
    idx = n_group_pairs+i;
    group_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{idx}.groups = group_pairs_pairs_string(i);
    venn_data{idx}.increase_features_group_1 = intersect(...
        venn_data{group_idx_1}.increase_features_group_1,...
        venn_data{group_idx_2}.increase_features_group_1);
    venn_data{idx}.increase_features_group_2 = intersect(...
        venn_data{group_idx_1}.increase_features_group_2,...
        venn_data{group_idx_2}.increase_features_group_2);
    venn_data{idx}.n_increase_features_group_1 = size(venn_data{idx}.increase_features_group_1,1);
    venn_data{idx}.n_increase_features_group_2 = size(venn_data{idx}.increase_features_group_2,1);
end
%% get shared DEFs between all groups - level 3 (ABC)
venn_data{end}.groups = join(group_pairs_string,", ");
venn_data{end}.increase_features_group_1 = feature_intersection_3(...
    venn_data{1}.increase_features_group_1,...
    venn_data{2}.increase_features_group_1,...
    venn_data{3}.increase_features_group_1);
venn_data{end}.increase_features_group_2 = feature_intersection_3(...
    venn_data{1}.increase_features_group_2,...
    venn_data{2}.increase_features_group_2,...
    venn_data{3}.increase_features_group_2);
venn_data{end}.n_increase_features_group_1 = length(venn_data{end}.increase_features_group_1);
venn_data{end}.n_increase_features_group_2 = length(venn_data{end}.increase_features_group_2);
%% correct counts
for i=1:(length(venn_data)-1)
    venn_data{i}.n_increase_features_group_1 = size(venn_data{i}.increase_features_group_1,1) - venn_data{end}.n_increase_features_group_1;
    venn_data{i}.n_increase_features_group_2 = size(venn_data{i}.increase_features_group_2,1) - venn_data{end}.n_increase_features_group_2;
end
for i=1:n_group_pairs
    idx = n_group_pairs+i;
    group_pair_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_pair_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{group_pair_idx_1}.n_increase_features_group_1 = venn_data{group_pair_idx_1}.n_increase_features_group_1 - venn_data{idx}.n_increase_features_group_1;
    venn_data{group_pair_idx_2}.n_increase_features_group_1 = venn_data{group_pair_idx_2}.n_increase_features_group_1 - venn_data{idx}.n_increase_features_group_1;
    venn_data{group_pair_idx_1}.n_increase_features_group_2 = venn_data{group_pair_idx_1}.n_increase_features_group_2 - venn_data{idx}.n_increase_features_group_2;
    venn_data{group_pair_idx_2}.n_increase_features_group_2 = venn_data{group_pair_idx_2}.n_increase_features_group_2 - venn_data{idx}.n_increase_features_group_2;
end
%% plot result
if options.plot
    venn_labels_increase_1 = cellfun(@(x) x.n_increase_features_group_1,venn_data,'uni',1);
    venn(3,'sets',group_pairs_string,'labels',venn_labels_increase_1,'alpha',0.7,'colors',autumn(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE INCREASE (GROUP 2)",'Position', [1.5 -0.3, 0])
    if options.save_path~=""
        savefig(options.save_path+"_increase.fig");
        saveas(gcf,options.save_path+"_increase.png");
    end
    venn_labels_increase_2 = cellfun(@(x) x.n_increase_features_group_2,venn_data,'uni',1);
    venn(3,'sets',group_pairs_string,'labels',venn_labels_increase_2,'alpha',0.7,'colors',winter(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE INCREASE (GROUP 1)",'Position', [1.5 -0.3, 0])
    if options.save_path~=""
        savefig(options.save_path+"_decrease.fig");
        saveas(gcf,options.save_path+"_decrease.png");
    end
    venn_labels_increase_all = venn_labels_increase_1+venn_labels_increase_2;
    venn(3,'sets',group_pairs_string,'labels',venn_labels_increase_all,'alpha',0.7,'colors',summer(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE CHANGE",'Position', [1.5 -0.3, 0])
    if options.save_path~=""
        savefig(options.save_path+"_all.fig");
        saveas(gcf,options.save_path+"_all.png");
    end
end
end

function venn_data = interpret_features_venn_4(feature_comparisons,options)
arguments
    feature_comparisons cell
    options.plot = false;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.save_path string = "";
end
groups = cellfun(@(x) x.groups,feature_comparisons,'uni',0);
groups = unique([groups{:}]);
%% data preparation
group_pairs_string = cellfun(@(x) join(x.groups, " vs "),feature_comparisons,'uni',1);
n_group_pairs = numel(feature_comparisons);
venn_data = cell(1,15);
log2_fold_change_threshold = log2(options.fold_change_threshold);
%% get comparison between two groups - level 1 (A, B, C, D)
for i=1:n_group_pairs
    venn_data{i}.groups = feature_comparisons{i}.groups;
    if ~isempty(feature_comparisons{i}.feature_table)
        idx_sig = all([abs(feature_comparisons{i}.feature_table.log2_fold_change)>log2_fold_change_threshold feature_comparisons{i}.feature_table.adj_p_value < options.alpha],2);
        idx_sig_up = all([feature_comparisons{i}.feature_table.log2_fold_change>0 idx_sig],2);
        idx_sig_down = all([feature_comparisons{i}.feature_table.log2_fold_change<0 idx_sig],2);
        venn_data{i}.increase_features_group_1 = feature_comparisons{i}.feature_table{idx_sig_up,"feature"};
        venn_data{i}.increase_features_group_2 = feature_comparisons{i}.feature_table{idx_sig_down,"feature"};
    else
        venn_data{i}.increase_features_group_1 = [];
        venn_data{i}.increase_features_group_2 = [];
    end
        venn_data{i}.n_increase_features_group_1 = size(venn_data{i}.increase_features_group_1,1);
        venn_data{i}.n_increase_features_group_2 = size(venn_data{i}.increase_features_group_2,1);
end
%% get shared DEFs between different group comparisons - level 2 (AB, AC, AD, BC, BC, CD)
group_pairs_pairs = nchoosek(group_pairs_string,2);
group_pairs_pairs_string = join(group_pairs_pairs, ",  ");
n_group_pairs_pairs = size(group_pairs_pairs,1);
for i=1:n_group_pairs_pairs
    idx = n_group_pairs+i;
    group_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{idx}.groups = group_pairs_pairs_string(i);
    venn_data{idx}.increase_features_group_1 = intersect(...
        venn_data{group_idx_1}.increase_features_group_1,...
        venn_data{group_idx_2}.increase_features_group_1);
    venn_data{idx}.increase_features_group_2 = intersect(...
        venn_data{group_idx_1}.increase_features_group_2,...
        venn_data{group_idx_2}.increase_features_group_2);
    venn_data{idx}.n_increase_features_group_1 = size(venn_data{idx}.increase_features_group_1,1);
    venn_data{idx}.n_increase_features_group_2 = size(venn_data{idx}.increase_features_group_2,1);
end
%% get shared DEFs between level 2 - level 3 (ABC, ABD, ACD, ACD, BCD)
group_pairs_triplets_idx = nchoosek([1 2 3 4],3);
group_pairs_triplets = nchoosek(group_pairs_string,3);
group_pairs_triplets_string = join(group_pairs_triplets, ", ");
n_group_pairs_triplets = size(group_pairs_triplets,1);
for i=1:n_group_pairs_triplets
    idx = n_group_pairs+n_group_pairs_pairs+i;
    venn_data{idx}.groups = group_pairs_triplets_string(idx);
    venn_data{idx}.increase_features_group_1 = feature_intersection_3(...
        venn_data{group_pairs_triplets_idx(1)}.increase_features_group_1,...
        venn_data{group_pairs_triplets_idx(2)}.increase_features_group_1,...
        venn_data{group_pairs_triplets_idx(3)}.increase_features_group_1);
    venn_data{idx}.increase_features_group_2 = feature_intersection_3(...
        venn_data{group_pairs_triplets_idx(1)}.increase_features_group_2,...
        venn_data{group_pairs_triplets_idx(2)}.increase_features_group_2,...
        venn_data{group_pairs_triplets_idx(3)}.increase_features_group_2);
    venn_data{idx}.n_increase_features_group_1 = length(venn_data{idx}.increase_features_group_1);
    venn_data{idx}.n_increase_features_group_2 = length(venn_data{idx}.increase_features_group_2);
end
%% get shared DEFs between level 3 - level 4 (ABCD)
venn_data{end}.groups = join(groups," vs ");
venn_data{end}.increase_features_group_1 = feature_intersection_4(...
    venn_data{1}.increase_features_group_1,...
    venn_data{2}.increase_features_group_1,...
    venn_data{3}.increase_features_group_1,...
    venn_data{4}.increase_features_group_1);
venn_data{end}.increase_features_group_2 = feature_intersection_4(...
    venn_data{1}.increase_features_group_2,...
    venn_data{2}.increase_features_group_2,...
    venn_data{3}.increase_features_group_2,...
    venn_data{4}.increase_features_group_2);
venn_data{end}.n_increase_features_group_1 = length(venn_data{end}.increase_features_group_1);
venn_data{end}.n_increase_features_group_2 = length(venn_data{end}.increase_features_group_2);
%% correct counts - level 4
for i=1:(length(venn_data)-1)
    venn_data{i}.n_increase_features_group_1 = size(venn_data{i}.increase_features_group_1,1) - venn_data{end}.n_increase_features_group_1;
    venn_data{i}.n_increase_features_group_2 = size(venn_data{i}.increase_features_group_2,1) - venn_data{end}.n_increase_features_group_2;
end
for i=1:n_group_pairs
    idx = n_group_pairs+i;
    group_pair_idx_1 = find(strcmp(group_pairs_string,group_pairs_pairs(i,1)));
    group_pair_idx_2 = find(strcmp(group_pairs_string,group_pairs_pairs(i,2)));
    venn_data{group_pair_idx_1}.n_increase_features_group_1 = venn_data{group_pair_idx_1}.n_increase_features_group_1 - venn_data{idx}.n_increase_features_group_1;
    venn_data{group_pair_idx_2}.n_increase_features_group_1 = venn_data{group_pair_idx_2}.n_increase_features_group_1 - venn_data{idx}.n_increase_features_group_1;
    venn_data{group_pair_idx_1}.n_increase_features_group_2 = venn_data{group_pair_idx_1}.n_increase_features_group_2 - venn_data{idx}.n_increase_features_group_2;
    venn_data{group_pair_idx_2}.n_increase_features_group_2 = venn_data{group_pair_idx_2}.n_increase_features_group_2 - venn_data{idx}.n_increase_features_group_2;
end
%% plot result
if options.plot
    venn_labels_increase_1 = cellfun(@(x) x.n_increase_features_group_1,venn_data,'uni',1);
    venn(4,'sets',group_pairs_string,'labels',venn_labels_increase_1,'alpha',0.7,'colors',autumn(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE INCREASE",'Position', [1.5 -0.3, 0])
    if options.save_path~=""
        savefig(options.save_path+"_increase.fig");
        saveas(gcf,options.save_path+"_increase.png");
    end
    venn_labels_increase_2 = cellfun(@(x) x.n_increase_features_group_2,venn_data,'uni',1);
    venn(4,'sets',group_pairs_string,'labels',venn_labels_increase_2,'alpha',0.7,'colors',winter(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE DECREASE",'Position', [1.5 -0.3, 0])
    if options.save_path~=""
        savefig(options.save_path+"_decrease.fig");
        saveas(gcf,options.save_path+"_decrease.png");
    end
    venn_labels_increase_all = venn_labels_increase_1+venn_labels_increase_2;
    venn(4,'sets',group_pairs_string,'labels',venn_labels_increase_all,'alpha',0.7,'colors',summer(4),'edgeC',[0 0 0],'labelC',[1 1 1],'edgeW',5);
    title("FEATURE VALUE CHANGE",'Position', [1.5 -0.3, 0])
    if options.save_path~=""
        savefig(options.save_path+"_all.fig");
        saveas(gcf,options.save_path+"_all.png");
    end
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
function v = feature_intersection_4(v1, v2, v3, v4)
    v_234 = feature_intersection_3(v2, v3, v4);
    if ~isempty(v_234)
        v = feature_intersection_3(v1, v_234);
    else
        v = [];
        return
    end
    if isempty(v)
        v = [];
    end
end
