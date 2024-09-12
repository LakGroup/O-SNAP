function [sample_names, group_idx] = sort_sample_to_groups(groups,voronoi_data)
sample_names = cellfun(@(x) x.name, voronoi_data, 'uni',1);
sample_groups = cellfun(@(x) x.group, voronoi_data, 'uni',1);
group_idx = cellfun(@(x) find(sample_groups==x), groups,'uni',0);
end



