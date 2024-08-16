function percentile_data = generate_percentile_data(work_dir, groups, voronoi_data)
n_groups = length(groups);
% create list of percentiles
percentiles = (1:19)*5;
n_percentiles = size(percentiles,2);
percentile_fields = arrayfun(@(elem) num2str(elem), percentiles, 'uni', 0);
percentile_data = cell(1,n_groups);
[sample_names, group_idx] = sort_sample_to_groups(groups,voronoi_data);

%% calculate percentiles to y and the return variable
for g=1:n_groups
    percentile_data{g}.group = groups{g};
    n_samp_in_group = length(group_idx{g});
    if n_samp_in_group > 0
        percentile_data{g}.sample_names = sample_names(group_idx{g});
        for p=1:n_percentiles
            percentile_data{g}.(['prctile_' percentile_fields{p}]) = arrayfun(@(x) prctile(voronoi_data{x}.reduced_log_voronoi_density,percentiles(p)), group_idx{g});
        end
    end
end

save(fullfile(work_dir,'percentile_data'),'percentile_data');
end
