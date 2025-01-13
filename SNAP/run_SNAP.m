function SNAP_data = run_SNAP(root_dir,analysis_name,groups,replicates,options)
    arguments
        root_dir char
        analysis_name char
        groups cell
        replicates cell
        options.generate_features logical = false %
        options.load_feature_table logical = true %
        options.load_feature_selection logical = false %
        options.load_PCA logical = false %
        options.run_classification logical = true %
        options.load_classification logical = false %
        options.suffix = ""
        options.feature_selection_threshold double = 0.5
        options.num_trials double = 3
        options.num_components double = 0.95
        options.compute_validation_accuracy logical = true
        options.k_folds double = 5
        options.verbose logical = false
        options.plot_PCA logical = false
        options.cluster_observations logical = false %
        options.alpha double = 0.05
        options.fold_change_threshold double = 2;
        options.cluster_features logical = false %
        options.graph_features logical = false %
        options.compare_features logical = true %
        options.venn logical = false %
        options.perform_ontology logical = false %
        options.plot_violin logical = false %
        options.plot_radial logical = false %
        options.run_FSEA logical = false %
        options.plot_ripley_k logical = false
    end
    SNAP_data.analysis_name = analysis_name;
    SNAP_data.groups = unique(groups);
    SNAP_data.replicates = unique(replicates);
    SNAP_data.options = options;
    SNAP_data.date = datetime;
    SNAP_classification_data = SNAP_data;
    SNAP_data.starttime = tic;

    %% prepare SNAP run
    work_dir = fullfile(root_dir,analysis_name);
    save_analysis_path = fullfile(work_dir,analysis_name+".mat");
    save_classification_path = fullfile(work_dir,analysis_name+"_classifiers.mat");
    log_file = sprintf('%s_%s.txt',analysis_name,datetime('now','Format','y-MM-dd'));
    logID = fopen(fullfile(work_dir,log_file),'w');
    log_output("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",logID)
    log_output("Running analysis for: " + analysis_name + "\n",logID)
    log_output(sprintf("RUN START: %s\n", string(datetime)),logID)
    log_output(sprintf("   GROUPS: %s\n", string(join(groups(:),', '))),logID)
    log_output(sprintf("   REPS: %s\n", string(join(replicates(:),', '))),logID)
    warning('off','all');
    close all
    %% generate SNAP feature table data
    try
        % generate_features flag takes precedence over load_features flag
        if options.generate_features
            log_output("  Generating features from scratch...\n",logID)
            % calculate features
            starttime_step = tic;
            generate_SNAP_features(work_dir,groups,replicates);
            log_output(sprintf("      Completed %s (%.0f min)...\n",string(datetime),toc(starttime_step)/60),logID);
            log_output("  Creating table from generated features...\n",logID)
            % coallate features
            starttime_step = tic;
            SNAP_data.feature_data = SNAP_nucleus_to_feature_table_batch(work_dir,groups,replicates);
            log_output(sprintf("      Completed %s (%.0f min)...\n",string(datetime),toc(starttime_step)/60),logID);
            save_SNAP_run(save_analysis_path,SNAP_data,'append',false)
        elseif options.load_feature_table
            % load features if file has proper information
            if exist(save_analysis_path,"file\n",logID)...
            && has_variables(save_analysis_path,{'feature_data'})
                load(save_analysis_path,'SNAP_data');
                SNAP_data.feature_data = SNAP_data;
                log_output("  Loaded from "+save_analysis_path + "...\n",logID)
            % if user desires to load features but they are not present,
            % attempt to create table under the assumption that features 
            % are generated
            else
                log_output("  Analysis file not found...\n",logID)
                try
                    log_output("  Creating table from generated features...\n",logID)
                    % coallate features
                    starttime_step = tic;
                    SNAP_data.feature_data = SNAP_nucleus_to_feature_table_batch(work_dir,groups,replicates);
                    log_output(sprintf("      Completed %s (%.0f min)...\n",string(datetime),toc(starttime_step)/60),logID);
                    save_SNAP_run(save_analysis_path,SNAP_data,'append',false)
                % if no files are found with data, then generate features
                % from scratch
                catch ME
                    % return
                    if (strcmp(ME.identifier,'SNAP:no_valid_files_found'))
                        log_output("  Generating features from scratch...\n",logID)
                        % calculate features
                        starttime_step = tic;
                        generate_SNAP_features(work_dir,groups,replicates);
                        log_output(sprintf("      Completed %s (%.0f min)...\n",string(datetime),toc(starttime_step)/60),logID);
                        log_output("  Creating table from generated features...\n",logID)
                        % coallate features
                        starttime_step = tic;
                        SNAP_data.feature_data = SNAP_nucleus_to_feature_table_batch(work_dir,groups,replicates);
                        log_output(sprintf("      Completed %s (%.0f min)...\n",string(datetime),toc(starttime_step)/60),logID);
                        save_SNAP_run(save_analysis_path,SNAP_data,'append',false)
                    else
                        rethrow(ME)
                    end
                 end
            end
        else
            % if user does not want to either load or generate features,
            % then run must end since feature data does not exist in the
            % workspace
            ME = MException('SNAP:no_feature_data_in_RUN', ...
                        sprintf('  No feature data in workspace, ending run...'));
            throw(ME);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID);
        return
    end

    %% filter for groups and replicates
    try
        log_output("  Filtering feature data table for desired groups and replicates...\n",logID)
        feature_data = prepare_voronoi_table_data(SNAP_data.feature_data,...
            "groups",groups,"replicates",replicates,...
            "keep_rep_sample_info",true);
        % end execution if not enough groups after filtering
        if numel(unique(feature_data.group)) < 2
            ME = MException('SNAP:insufficient_group_number', ...
                            sprintf('Only one group found: %s',groups{1}));
            throw(ME)
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end

    %% feature selection
    try
        if options.load_feature_selection && has_variables(save_analysis_path,{'vars_selected','var_select_result'})
                load(save_analysis_path,'vars_selected','var_select_result')
                SNAP_data.vars_selected = vars_selected;
                SNAP_data.var_select_result = var_select_result;
                log_output("  Feature selection results loaded\n",logID)
        else
            log_output("  Selecting features...\n",logID)
            [SNAP_data.vars_selected, SNAP_data.var_select_result] = select_features_SNAP(feature_data,...
                "threshold",options.feature_selection_threshold,...
                "verbose",options.verbose);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    save_SNAP_run(save_analysis_path,SNAP_data)
    %% PCA
    try
        if options.load_PCA && has_variables(save_analysis_path,{'pca_result'})
            load(save_analysis_path,'pca_result')
            SNAP_data.pca_result = pca_result;
            log_output("  PCA results loaded\n",logID)
        else
            if options.suffix ~= ""
                save_path = fullfile(work_dir, join(['PCA',groups,options.suffix],'_'));
            else
                save_path = fullfile(work_dir, join(['PCA',groups],'_'));
            end
            log_output("  Performing PCA...\n",logID)
            SNAP_data.pca_result = run_SNAP_PCA(feature_data,...
                "vars_sel", SNAP_data.vars_selected,...
                "save_path",save_path);
            save_SNAP_run(save_analysis_path,SNAP_data)
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% run classification
    try
        if options.run_classification
            if ~all(isfield(SNAP_classification_data,{'classifiers','classification_summary'}))...
                    || (SNAP_data.options.num_trials ~= options.num_trials)...
                    || ~options.load_classification
                log_output("  Running classification...\n",logID)
                [SNAP_classification_data.classifiers,SNAP_classification_data.classification_summary] = run_SNAP_classification_batch(...
                    feature_data.group,SNAP_data.pca_result,SNAP_data.vars_selected,...
                    "num_trials",options.num_trials,...
                    "num_components",options.num_components,...
                    "compute_validation_accuracy",options.compute_validation_accuracy,...
                    "k_folds", options.k_folds,...
                    "verbose",options.verbose);
                SNAP_data.classification_summary = SNAP_classification_data.classification_summary;
                save_SNAP_run(save_analysis_path,SNAP_data)
                save_SNAP_run(save_classification_path,SNAP_classification_data)
            end
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% cluster observations
    if options.cluster_observations
        log_output("  Clustering observations...\n",logID)
        try
            if options.suffix ~= ""
                save_path = fullfile(work_dir, join(['clustering_obs',groups,options.suffix],'_'));
            else
                save_path = fullfile(work_dir, join(['clustering_obs',groups],'_'));
            end
            if ~all(isfield(SNAP_data,{'cluster_observations_idx_PCA','cluster_observations_idx_tSNE'}))
                [SNAP_data.cluster_observations_idx_PCA,SNAP_data.cluster_observations_idx_tSNE] = cluster_observations_hierarchical(feature_data,SNAP_data.vars_sel,...
                    "save_path",save_path);
            end
        catch ME
            handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
            return
        end
    end
    %% cluster features  
    try
        if options.cluster_features
            log_output("  Clustering features...\n",logID)
            [SNAP_data.cluster_features_idx_PCA,SNAP_data.cluster_features_idx_tSNE] = cluster_features_hierarchical(feature_data, vars_sel,...
                "save_path",fullfile(work_dir,join(['clustering_feat',groups,options.suffix],'_')));
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% create graph of features related by correlation
    try
        if options.graph_features
            log_output("  Creating feature graph...\n",logID)
            interpret_features_graph(feature_data,"save_path",fullfile(work_dir,"feature_graph\n",logID));
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% compare all groups to each other
    try
        if options.compare_features || options.venn || options.perform_ontology
            log_output("  Comparing features...\n",logID)
            if options.suffix ~= ""
                save_path = fullfile(work_dir, "features_volcano_"+suffix);
            else
                save_path = fullfile(work_dir, "features_volcano");
            end
            SNAP_data.feature_comparisons = compare_groups(SNAP_data.feature_data,...
                "alpha",options.alpha,...
                "fold_change_threshold",options.fold_change_threshold,...
                "save_path",save_path);
            close all
            save_SNAP_run(save_analysis_path,SNAP_data)
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% venn diagram
    try
        if options.venn
            log_output("  Creating venn diagrams...\n",logID)
            SNAP_data.venn_data = interpret_features_venn(feature_data,SNAP_data.feature_comparisons,"save_path",fullfile(work_dir, "features_venn"));
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% perform ontology
    try
        if options.perform_ontology
            log_output("  Performing feature ontology analysis...\n",logID)
            interpret_features_ontology(feature_data,"plot",false);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% violin plots
    try
        if options.plot_violin
            log_output("  Plotting violin plots...\n",logID)
            plot_violin(feature_data,work_dir)
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% plot radial
    try
        if options.plot_radial
            log_output("  Plotting radial densities...\n",logID)
            plot_ellipse_set_batch(work_dir,groups,replicates,feature_data);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% FSEA
    try
        if options.run_FSEA
            log_output("  Running feature set enrichment analysis\n",logID)
            run_SNAP_MrFSEA(work_dir, feature_data, ["feature_universe_1","feature_universe_2","feature_universe_3","feature_universe_4"]);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end
    %% ripley K
    try
        if options.plot_ripley_k
            log_output("  Plotting Ripley's K...\n",logID)
            ripley_k(work_dir,groups,replicates);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,logID)
        return
    end

    % Finish run
    conclude_SNAP_run(save_analysis_path, SNAP_data, logID)
end

function conclude_SNAP_run(path,data,logID)
    save_SNAP_run(path,data)
    runtime = toc(data.starttime)/60;
    log_output(sprintf("RUN END: %s\n", string(datetime)),logID)
    log_output(sprintf("TIME ELAPSED: %.2f min\n", runtime),logID)
    log_output("Done!\n",logID)
    log_output("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",logID)
    fclose(logID);
    warning('on','all');
    close all
end

function handle_SNAP_error(ME,save_analysis_path,data,logID)
    if nargin == 3
        log_output("- - - - - - - - - - - SNAP ERROR - - - - - - - - - - -\n")
        log_output(sprintf("%s\n",getReport(ME,'extended','hyperlinks','off')))
        log_output("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
    elseif nargin == 4
        log_output("- - - - - - - - - - - SNAP ERROR - - - - - - - - - - -\n",logID)
        log_output(sprintf("%s\n",getReport(ME,'extended','hyperlinks','off')),logID)
        log_output("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",logID)
    end
    conclude_SNAP_run(save_analysis_path,data,logID);
end

function save_SNAP_run(path,data,options)
arguments
    path string
    data struct
    options.check_overwrite = false
    options.append = true
end
% if structure is empty, return
if isempty(fieldnames(data))
    return
end

% append option
if options.append
    if exist(path,'file')
        save(path,'-struct','data','-append');
    else
        disp("Warning: File does not exist, creating file: "+path)
        save(path,'-struct','data');
    end
% overwrite file
elseif options.check_overwrite
    if exist(path,'file')
        overwrite_dialog();
    else
        save(path,'-struct','data');
    end
end
    
    function overwrite_dialog()
        response = questdlg(sprintf('Warning: %s\nalready exists. Overwrite?',path),...
            'Warning: Overwrite',...
            'Overwrite','New file...','Cancel','Cancel');
        switch response
            case 'Overwrite'
                disp('      Overwriting...')
                save(path,'-struct','data');
            case 'New file...'
                new_path = uigetfile();
                save(new_path,'-struct','data');
            case 'Cancel'
                disp('Save canceled')
        end
    end
end