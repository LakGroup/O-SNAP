function T = SNAP_nucleus_to_feature_table_batch(work_dir,groups,reps,options)
arguments
    work_dir string
    groups cell
    reps cell
    options.load logical = true
    options.n_workers double = maxNumCompThreads
end
    vars_string = {'group','biological_replicate','name'};

    %% organize into tables
    % get data
    data_info_table = get_valid_SNAP_nucleus_files(work_dir,groups,reps,{'voronoi_cluster_radius'});
    % split data for parallel processing
    [split_data, options.n_workers] = split_files_for_parallel(data_info_table, options.n_workers, true);
    T_cell = cell(1,options.n_workers);
    % run analysis
    parfor p=1:options.n_workers
        data_info_table_p = split_data{p};
        T_cell{p} = SNAP_nucleus_to_feature_table(data_info_table_p);
    end
    T = sortrows(vertcat(T_cell{:}),vars_string);
end


