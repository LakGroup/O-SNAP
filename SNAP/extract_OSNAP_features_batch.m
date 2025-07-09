function T = extract_OSNAP_features_batch(work_dir,groups,replicates,options)
arguments
    work_dir string
    groups cell
    replicates cell
    options.load logical = true
    options.n_workers double = maxNumCompThreads
end
    vars_string = {'group','biological_replicate','name'};

    %% organize into tables
    % get data
    data_info_table = get_valid_OSNAP_samples(work_dir,groups,replicates,{'voronoi_cluster_radius'});
    % split data for parallel processing
    [split_file_list, options.n_workers] = split_data_to_n_OSNAP(data_info_table, options.n_workers);
    T_cell = cell(1,options.n_workers);
    % run analysis
    parfor p=1:options.n_workers
        data_info_table_p = split_file_list{p};
        T_cell{p} = extract_OSNAP_features_sample(data_info_table_p);
    end
    T = sortrows(vertcat(T_cell{:}),vars_string);
end


