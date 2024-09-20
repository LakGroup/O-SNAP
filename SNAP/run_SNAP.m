function T = run_SNAP(root_dir,analysis_name,groups,reps,options)
    arguments
        root_dir char
        analysis_name char
        groups cell
        reps cell
        options.load logical = false
        options.plot_PCA logical = true
        options.plot_features logical = true
        options.plot_radial logical = false
        options.plot_ripley_k logical = false
        options.run_GSEA logical = true
    end
    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    disp("Running analysis for: " + analysis_name)
    disp(['  RUN START: ' char(datetime)])
    disp("   GROUPS: " + string(join(groups(:),', ')))
    disp("   REPS: " + string(join(reps(:),', ')))
    warning('off','all');
    close all

    %% file management
    work_dir = fullfile(root_dir,analysis_name);
    table_file_path = fullfile(work_dir,analysis_name+".mat");

    %% calculate T
    if options.load && exist(table_file_path,"file")
        disp("Loading from "+table_file_path + "...")
        load(table_file_path,"T");
    elseif options.load && ~exist("table_file_path","file")
        disp("Analysis file not found...")
        try
            disp("Creating table from generated features...")
            % coallate features
            T = voronoi_data_to_table_batch(work_dir,groups,reps);
        catch ME
            if (strcmp(ME.identifier,'SNAP:no_valid_files_found'))
                disp("Generating features from scratch...")
                % calculate features
                generate_SNAP_features(work_dir,groups,reps);
                disp("Creating table from generated features...")
                % coallate features
                T = voronoi_data_to_table_batch(work_dir,groups,reps);
            else
                rethrow(ME)
            end
        end
    else
        disp("Generating features from scratch...")
        % calculate features
        generate_SNAP_features(work_dir,groups,reps);
        disp("Creating table from generated features...")
        % coallate features
        T = voronoi_data_to_table_batch(work_dir,groups,reps,"load",false);
    end
    %% filter T
    % select only desired groups and replicates for analysis
    T = T(all([ismember(T.biological_replicate,reps), ismember(T.group,groups)],2),:);
    % remove rows with only NaNs
    T = T(~all(isnan(T{:,4:end}),2),:);
    %% generate plots
    if options.plot_PCA
        disp('Plotting variable selection and PCA...')
        if length(unique(T.group)) >= 2
            plot_SNAP(T, groups, work_dir);
        else
            disp('Not enough groups in table')
        end
        close all
    end
    if options.plot_features
        % disp("Creating feature graph...")
        % interpret_features_graph(T, work_dir);
        % disp("Creating feature venn diagram...")
        % interpret_features_venn(T, work_dir);
        disp("Performing feature ontology analysis...")
        interpret_features_ontology(T,work_dir);
        close all
    end
    if options.plot_radial
        plot_radial_densities(data_info_table, T);
        close all
    end
    if options.plot_ripley_k
        ripley_k(work_dir,groups,reps);
        close all
    end
    if options.run_GSEA
        run_SNAP_MrGSEA(root_dir, analysis_name, T, ["feature_universe_1","feature_universe_2","feature_universe_3"]);
        close all
    end
    warning('on','all');
    close all
end
