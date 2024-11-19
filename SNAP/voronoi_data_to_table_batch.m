function T = voronoi_data_to_table_batch(work_dir,groups,reps,options)
arguments
    work_dir string
    groups cell
    reps cell
    options.load logical = true
    options.n_workers double = maxNumCompThreads
end

    %% parameters  
    dir_path_split = regexp(work_dir,'\','split');
    analysis_name = dir_path_split(end);
    vars_string = {'group','biological_replicate','name'};

    %% organize into tables
    tic
    % get data
    data_info_table = get_valid_voronoi_data(work_dir,groups,reps,{'cluster_area'});
    % split data for parallel processing
    [split_data, options.n_workers] = split_files_for_parallel(data_info_table, options.n_workers, true);
    T_cell = cell(1,options.n_workers);
    % run analysis
    for p=1:options.n_workers
        data_info_table_p = split_data{p};
        T_cell{p} = voronoi_data_to_table(data_info_table_p);
    end
    T = sortrows(vertcat(T_cell{:}),vars_string);
    disp("      Saving to " + fullfile(work_dir, analysis_name+".mat") + "...")
    save(fullfile(work_dir, analysis_name+".mat"),"T");
    toc
    disp(['      Completed: ' char(datetime)])

end


