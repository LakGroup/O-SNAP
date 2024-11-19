function T = run_SNAP(root_dir,analysis_name,groups,reps,options)
    arguments
        root_dir char
        analysis_name char
        groups cell
        reps cell
        options.load_table logical = true
        options.generate_features logical = false
        options.plot_PCA logical = true
        options.plot_clustering logical = false
        options.plot_features logical = false
        options.plot_violin logical = false
        options.plot_radial logical = false
        options.run_FSEA logical = false
        options.plot_ripley_k logical = false
    end
    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    disp("Running analysis for: " + analysis_name)
    disp(['RUN START: ' char(datetime)])
    disp("   GROUPS: " + string(join(groups(:),', ')))
    disp("   REPS: " + string(join(reps(:),', ')))
    warning('off','all');
    close all

    %% file management
    work_dir = fullfile(root_dir,analysis_name);
    table_file_path = fullfile(work_dir,analysis_name+".mat");

    %% calculate T
    % try
        if options.generate_features
            disp("  Generating features from scratch...")
            % calculate features
            generate_SNAP_features(work_dir,groups,reps);
            disp("  Creating table from generated features...")
            % coallate features
            T = voronoi_data_to_table_batch(work_dir,groups,reps);
        elseif options.load_table && exist(table_file_path,"file")
            disp("  Loading from "+table_file_path + "...")
            load(table_file_path,"T");
        else
            if ~options.load_table
                disp("  Analysis file not found...")
            end
            try
                disp("  Creating table from generated features...")
                % coallate features
                T = voronoi_data_to_table_batch(work_dir,groups,reps);
            catch ME
                % return
                if (strcmp(ME.identifier,'SNAP:no_valid_files_found'))
                    disp("  Generating features from scratch...")
                    % calculate features
                    generate_SNAP_features(work_dir,groups,reps);
                    disp("  Creating table from generated features...")
                    % coallate features
                    T = voronoi_data_to_table_batch(work_dir,groups,reps);
                else
                    rethrow(ME)
                end
             end
        end
        %% filter T
        % select only desired groups and replicates for analysis
        T = T(all([ismember(T.biological_replicate,reps), ismember(T.group,groups)],2),:);
        % remove rows with only NaNs
        T = T(~all(isnan(T{:,4:end}),2),:);
        %% generate plots
        if options.plot_PCA
            disp('  Plotting variable selection and PCA...')
            if length(unique(T.group)) >= 2
                plot_SNAP(T, groups, work_dir,...
                    'cluster_observations',options.plot_clustering,...
                    'cluster_features',options.plot_features);
            else
                disp('  Not enough groups in table')
            end
            close all
        end
        if options.plot_features
            disp("  Creating feature graph...")
            try
                interpret_features_graph(T, work_dir);
                disp("  Creating feature venn diagram...")
            catch ME
                disp(getReport(ME))
            end
            interpret_features_venn(T, work_dir);
            disp("  Performing feature ontology analysis...")
            interpret_features_ontology(T,work_dir);
            close all
        end
        if options.plot_violin
            disp("  Plotting violin plots...")
            plot_violin(work_dir,T)
        end
        if options.plot_radial
            disp("  Plotting radial densities...")
            plot_ellipse_set_batch(work_dir,groups,reps,T);
            close all
        end
        if options.run_FSEA
            run_SNAP_MrFSEA(root_dir, analysis_name, T, ["feature_universe_1","feature_universe_2","feature_universe_3","feature_universe_4"]);
            close all
        end
        if options.plot_ripley_k
            disp("  Plotting Ripley's K...")
            ripley_k(work_dir,groups,reps);
            close all
        end
    % catch ME
    %     % rethrow(ME)
    %     disp(ME.message)
    % end
    disp(['RUN End: ' char(datetime)])
    disp('Done!')
    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    warning('on','all');
    close all
end