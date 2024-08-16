function T = voronoi_data_to_table_batch(work_dir,groups,reps)

    %% parameters  
    dir_path_split = regexp(work_dir,'\','split');
    analysis_name = dir_path_split(end);
    n_processes = 12;
    vars_string = {'group','biological_replicate','name'};

    %% organize into tables
    if(~exist("T","var"))
        if exist(fullfile(work_dir, analysis_name+".mat"),'file')
            disp("Loading " + fullfile(work_dir, analysis_name+".mat"))
            load(fullfile(work_dir, analysis_name+".mat"),"T")
        else
            disp('Generating table data...')
            tic
            data_info_table = get_valid_voronoi_data(work_dir,groups,reps,{'cluster_area'});
            % split data for parallel processing
            [split_data, n_processes] = split_files_for_parallel(data_info_table, n_processes, true);
            % start parallel pool
            p_pool = gcp('nocreate');
            if isempty(p_pool)
                parpool("Processes",n_processes);
            end
            T_cell = cell(1,n_processes);
            % run analysis
            for p=1:n_processes
                data_info_table_p = split_data{p};
                T_cell{p} = voronoi_data_to_table(data_info_table_p);
            end
            T = sortrows(vertcat(T_cell{:}),vars_string);
            disp("Saving to " + fullfile(work_dir, analysis_name+".mat") + "...")
            save(fullfile(work_dir, analysis_name+".mat"),"T");
            toc
            disp(['Completed: ' char(datetime)])
        end
    end

end

