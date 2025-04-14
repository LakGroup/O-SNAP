function plot_ellipse_set_batch(T,work_dir,groups,reps,options)
arguments
    T table
    work_dir string
    groups cell
    reps cell
    options.n_processes double= 12;
end
%% get data
data_info_table = get_valid_SNAP_nucleus_files(work_dir,groups,reps,{'radial_loc_density','radial_dbscan_cluster_density'});
%% split data for parallel processing
[split_data, n_processes] = split_data_to_n(data_info_table,options.n_processes,"shuffle",true);
%% plot
max_radial_density = max(table2array(T(:,contains(T.Properties.VariableNames, 'radial_loc_density'))),[],'all');
max_radial_hetero_density = max(table2array(T(:,contains(T.Properties.VariableNames, 'radial_dbscan_cluster_density'))),[],'all');
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        try
            plot_ellipse_set_locs_and_hetero(filepath, max_radial_density, max_radial_hetero_density);
        catch ME
            disp(getReport(ME))
        end
    end
end
end

