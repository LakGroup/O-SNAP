function SNAP_data = run_SNAP(root_dir,analysis_name,groups,replicates,options)
    arguments
        root_dir char
        analysis_name char
        groups cell
        replicates cell
        % run options (will overwrite)
        options.run_generate_features logical = 0
        options.run_comparison logical = 0
        options.run_generate_batches logical = 0
        options.run_feature_selection logical = 0
        options.run_PCA logical = 1
        options.run_classification logical = 1
        options.run_cluster_observations logical = 0
        options.run_cluster_features logical = 0 
        options.run_graph_features logical = 0 
        options.run_venn logical = 0 
        options.run_ontology logical = 0 
        options.run_plot_violin logical = 0 
        options.run_plot_radial logical = 0 
        options.run_FSEA logical = 0 
        options.run_ripley_k logical = 0        
        % generate batch options
        options.split_test_train_by string = "replicate"
        options.test_train_ratio double = 0.2;
        options.test_train_k double = 5;
        % feature selection options
        options.feature_selection_threshold double = 0.5
        options.feature_selection_methods cell = {...
            'biPLS',...
            'rPLS',...
            'ChiSquare',...
            'MRMR',...
            'ReliefF',...
            'Boruta',...
            'biPLS_PLSToolbox',... %(Interval PLS from the PLS Toolbox; Eigenvector Inc.)
            'rPLS_PLSToolbox',... %(Recursive PLS from the PLS Toolbox; Eigenvector Inc.)
            'iPLS_PLSToolbox',... %(Interval PLS from the PLS Toolbox; Eigenvector Inc.)
            'GA_PLSToolbox'... %(Genetic Algorithm from the PLS Toolbox; Eigenvector Inc.)
            }
        % PCA options
        options.num_components_explained double = 0.60
        options.plot_PCA logical = 0
        % comparison options
        options.alpha double = 0.05
        options.fold_change_threshold double = 2;
        % other
        options.n_processes double = 12
        options.verbose logical = 0
        % save options
        options.suffix = ""
        options.save = true;
        options.save_if_error = true;
        options.check_overwrite = 0
    end

    %% prepare SNAP run
    work_dir = fullfile(root_dir,analysis_name);
    save_analysis_path = fullfile(work_dir,analysis_name+".mat");
    log_file = fullfile(work_dir,sprintf('%s_%s.txt',analysis_name,datetime('now','Format','y-MM-dd_HH-mm-ss')));
    diary(log_file)
    fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
    fprintf("Running analysis for: " + analysis_name + "\n")
    fprintf("RUN START: %s\n", string(datetime))
    fprintf("   GROUPS: %s\n", string(join(groups(:),', ')))
    fprintf("   REPS: %s\n", string(join(replicates(:),', ')))
    warning('off','all');
    close all

    try
        if exist(save_analysis_path,"file")
            starttime_step = tic;
            fprintf("  Loading data from %s...\n",save_analysis_path)
            try
                SNAP_data = load(save_analysis_path);
            catch ME
                SNAP_data = struct;
                handle_SNAP_error(ME,save_analysis_path,SNAP_data,"stop_run",0,"save",0);
                options.run_generate_features = 1;
            end
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
            SNAP_data.analysis_name = analysis_name;
            SNAP_data.groups = unique(groups);
            SNAP_data.replicates = unique(replicates);
            SNAP_data.date = datetime;
            SNAP_data.starttime = tic;
            SNAP_data.options = options;
            % generate SNAP feature table data if not loaded
            if isfield(SNAP_data,'feature_data')
                feature_data = SNAP_data.feature_data;
            else
                options.run_generate_features = 1;
            end
            % run feature comparisons if not loaded
            if ~isfield(SNAP_data,'feature_comparisons')
                options.run_comparison = 1;
            end
            % run steps prior to classification
            if options.run_classification
                if options.run_generate_features
                    options.run_generate_batches = 1;
                    options.run_feature_selection = 1;
                    options.run_PCA = 1;
                end
                if ~all(isfield(SNAP_data,{'train_idxs','test_idxs'}))
                    options.run_generate_batches = 1;
                elseif isempty(SNAP_data.train_idxs) || isempty(SNAP_data.test_idxs)
                    options.run_generate_batches = 1;
                end
                if ~all(isfield(SNAP_data,{'vars_selected','var_select_result'}))
                    options.run_feature_selection = 1;
                elseif isempty(SNAP_data.vars_selected) || isempty(SNAP_data.var_select_result)
                    options.run_feature_selection = 1;
                elseif ~isfield(SNAP_data,'vars_selected_summary')
                    SNAP_data.vars_selected_summary = summarize_vars_selected(...
                            SNAP_data.feature_data,...
                            SNAP_data.vars_selected,...
                            SNAP_data.var_select_result);
                end
                if ~isfield(SNAP_data,'pca_result')
                    options.run_PCA = 1;
                elseif isempty(SNAP_data.pca_result)
                    options.run_PCA = 1;
                end
            end
            SNAP_data.options = options;
        end
    catch ME
        options.run_generate_features = 1;
        options.run_comparison = 1;
        options.run_generate_batches = 1;
        options.run_feature_selection = 1;
        options.run_PCA = 1;
        options.run_classification = 1;
        SNAP_data.options = options;
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
    end


    %% generate SNAP feature table data
    try
        % generate_features flag takes precedence over load_features flag
        if options.run_generate_features
            fprintf("  Generating features from scratch...\n")
            % calculate features
            starttime_step = tic;
            generate_SNAP_features(work_dir,groups,replicates);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
            fprintf("  Creating table from generated features...\n")
            % coallate features
            starttime_step = tic;
            feature_data = SNAP_nucleus_to_feature_table_batch(work_dir,groups,replicates);
            SNAP_data.feature_data = feature_data;
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
            save_SNAP_run(save_analysis_path,SNAP_data,'append',false)
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    
    %% run classification steps if needed
    % filter data
    try
        feature_data_filtered = filter_SNAP_feature_data(feature_data,groups,replicates);
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    % generate train and test batches
    try
        if options.run_generate_batches
            fprintf("  Creating train/test batches...\n")
            starttime_step = tic;
            [SNAP_data.train_idxs, SNAP_data.test_idxs,options.split_test_train_by] = generate_train_test_sets(...
                feature_data_filtered,...
                "by",options.split_test_train_by,...
                "k",options.test_train_k,...
                "proportion",options.test_train_ratio);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    % feature selection
    try
        if options.run_feature_selection
            % clear variable for new run
            if isfield(SNAP_data,"vars_selected")
                SNAP_data = rmfield(SNAP_data,"vars_selected");
            end
            if isfield(SNAP_data,'var_selected_results')
                SNAP_data = rmfield(SNAP_data,"var_select_result");
            end
            fprintf("  Selecting features...\n")
            starttime_step = tic;
            % split batches for parallel
            train_idx_p = split_data_to_n(SNAP_data.train_idxs,options.n_processes,"shuffle",false);
            n_processes = numel(train_idx_p);
            feature_data_p = cell(1,n_processes);
            vars_selected_p = cell(1,n_processes);
            var_select_result_p = cell(1,n_processes);
            for p=1:n_processes
                feature_data_p{p} = feature_data_filtered;
            end
            parfor p=1:n_processes
                train_idx_b = train_idx_p{p};
                for b=1:numel(train_idx_b)
                    train_data = feature_data_p{p}(train_idx_b{b},:);
                    [vars_selected_p{p}{b}, var_select_result_p{p}{b}] = select_features_SNAP(train_data);
                end
            end
            % bring together batches from each parallel process
            SNAP_data.vars_selected = [vars_selected_p{:}];
            SNAP_data.var_select_result = [var_select_result_p{:}];
            SNAP_data.vars_selected_summary = summarize_vars_selected(...
                    feature_data_filtered,...
                    SNAP_data.vars_selected,...
                    SNAP_data.var_select_result);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% PCA
    try
        if options.run_PCA
            save_path = fullfile(work_dir, join(['PCA',groups],'_'));
            if options.suffix ~= ""
                save_path = join([save_path suffix],'_');
            end
            fprintf("  Performing PCA...\n")
            starttime_step = tic;
            % split batches for parallel
            num_components_explained = options.num_components_explained;
            train_idx_p = split_data_to_n(SNAP_data.train_idxs,options.n_processes,"shuffle",false);
            n_processes = numel(train_idx_p);
            if ~isempty(SNAP_data.vars_selected)
                vars_selected_p = split_data_to_n(SNAP_data.vars_selected,n_processes,"shuffle",false);
            end
            pca_result_p = cell(1,n_processes);
            if ~exist('feature_data_p','var')
                feature_data_p = cell(1,n_processes);
                for p=1:n_processes
                    feature_data_p{p} = feature_data;
                end
            end
            parfor p=1:n_processes
                train_idx_b = train_idx_p{p};
                for b=1:numel(train_idx_b)
                    train_data = feature_data_p{p}(train_idx_b{b},:);
                    pca_result_p{p}{b} = run_SNAP_PCA(train_data,...
                        "vars_sel",vars_selected_p{p}{b},...
                        "save_path",save_path,...
                        "num_components_explained",num_components_explained);
                end
            end
            % bring together batches from each parallel process
            SNAP_data.pca_result = [pca_result_p{:}];
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
    end

    %% run classification
    try
        if options.run_classification
            fprintf("  Running classification...\n")
            starttime_step = tic;
            % cell_dims = [n_batch,numel(SNAP_data.vars_selected),numel(SNAP_data.pca_result)];
            % assert(all(cell_dims == cell_dims(1)),'SNAP:classification_length_mismatch','number of batches not consistent with classification structures (feature selection, PCA)')
            % split batches for parallel
            train_idx_p = split_data_to_n(SNAP_data.train_idxs,options.n_processes,"shuffle",false);
            test_idx_p = split_data_to_n(SNAP_data.test_idxs,options.n_processes,"shuffle",false);
            n_processes = numel(train_idx_p);
            if ~isempty(SNAP_data.vars_selected)
                vars_selected_p = split_data_to_n(SNAP_data.vars_selected,n_processes,"shuffle",false);
            else
                vars_selected_p = split_data_to_n(cell(1,numel(SNAP_data.train_idxs)),n_processes,"shuffle",false);
            end
            if isfield(SNAP_data,'pca_result')
                if ~isempty(SNAP_data.pca_result)
                    pca_result_p = split_data_to_n(SNAP_data.pca_result,n_processes,"shuffle",false);
                else
                    pca_result_p = split_data_to_n(cell(1,numel(SNAP_data.train_idxs)),n_processes,"shuffle",false);
                end
            else
                pca_result_p = split_data_to_n(cell(1,numel(SNAP_data.train_idxs)),n_processes,"shuffle",false);
            end
            if ~exist('feature_data_p','var')
                feature_data_p = cell(1,n_processes);
                for p=1:n_processes
                    feature_data_p{p} = feature_data;
                end
            end
            classifiers_p = cell(1,n_processes);
            for p=1:n_processes
                train_idx_b = train_idx_p{p};
                test_idx_b = test_idx_p{p};
                for b=1:numel(train_idx_b)
                    train_data = feature_data_p{p}(train_idx_b{b},:);
                    test_data = feature_data_p{p}(test_idx_b{b},:);
                    classifiers_p{p}{b} =...
                        run_SNAP_classification_batch(...
                        train_data,...
                        'vars_selected',vars_selected_p{p}{b},...
                        'pca_result',pca_result_p{p}{b},....
                        'test_data',test_data,...
                        'verbose',0);
                end
            end
            % bring together batches from each parallel process
            classifiers = [classifiers_p{:}];
            model_type = reshape(cellfun(@(x) x.ModelType, [classifiers{:}], 'uni', 1),[],1);
            test_accuracy = reshape(cellfun(@(x) x.TestAccuracy, [classifiers{:}], 'uni', 1),[],1);
            classification_summary = table(model_type,test_accuracy);
            classification_summary = sortrows(groupsummary(classification_summary,"model_type",["mean","std"]),"mean_test_accuracy","descend");
            SNAP_data.classifiers = classifiers;
            SNAP_data.classification_summary = classification_summary;
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
    end
   
   %% cluster observations
   try
        if options.run_cluster_observations
            fprintf("  Clustering observations...\n")
            starttime_step = tic;
            save_path = fullfile(work_dir, join(['clustering_observations',groups],'_'));
            if options.suffix ~= ""
                save_path = join([save_path suffix],'_');
            end
            if ~all(isfield(SNAP_data,{'cluster_observations_idx_PCA','cluster_observations_idx_tSNE'}))
                [SNAP_data.cluster_observations_idx_PCA,SNAP_data.cluster_observations_idx_tSNE] = cluster_observations_hierarchical(feature_data,SNAP_data.vars_sel,...
                    "save_path",save_path);
            end
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% cluster features  
    try
        if options.run_cluster_features
            fprintf("  Clustering features...\n")
            starttime_step = tic;
            save_path = fullfile(work_dir, join(['clustering_features',groups],'_'));
            if options.suffix ~= ""
                save_path = join([save_path suffix],'_');
            end
            [SNAP_data.cluster_features_idx_PCA,SNAP_data.cluster_features_idx_tSNE] = cluster_features_hierarchical(feature_data, vars_sel,...
                "save_path",fullfile(work_dir,save_path));
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% create graph of features related by correlation
    try
        if options.run_graph_features
            fprintf("  Creating feature graph...\n")
            starttime_step = tic;
            interpret_features_graph(feature_data,"save_path",fullfile(work_dir,"feature_graph\n"));
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% compare all groups to each other
    try
        if options.run_comparison || options.run_venn || options.run_ontology
            fprintf("  Comparing features...\n")
            starttime_step = tic;
            save_path = fullfile(work_dir, "features_volcano");
            if options.suffix ~= ""
                save_path = join(save_path,options.suffix,"_");
            end
            SNAP_data.feature_comparisons = compare_groups(feature_data,...
                "alpha",options.alpha,...
                "fold_change_threshold",options.fold_change_threshold,...
                "save_path",save_path);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
            save_SNAP_run(save_analysis_path,SNAP_data)
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% venn diagram
    try
        if options.run_venn
            fprintf("  Creating venn diagrams...\n")
            starttime_step = tic;
            SNAP_data.venn_data = interpret_features_venn(feature_data,SNAP_data.feature_comparisons,"save_path",fullfile(work_dir, "features_venn"));
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% perform ontology
    try
        if options.run_ontology
            fprintf("  Performing feature ontology analysis...\n")
            starttime_step = tic;
            interpret_features_ontology(feature_data,"plot",false);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% violin plots
    try
        if options.run_plot_violin
            fprintf("  Plotting violin plots...\n")
            starttime_step = tic;
            plot_violin(feature_data,work_dir)
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% plot radial
    try
        if options.run_plot_radial
            fprintf("  Plotting radial densities...\n")
            plot_ellipse_set_batch(work_dir,groups,replicates,feature_data);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% FSEA
    try
        if options.run_FSEA
            fprintf("  Running feature set enrichment analysis...\n")
            starttime_step = tic;
            run_SNAP_MrFSEA(work_dir, feature_data, ["feature_universe_1","feature_universe_2","feature_universe_3","feature_universe_4"]);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% ripley K
    try
        if options.run_ripley_k
            fprintf("  Plotting Ripley's K...\n")
            ripley_k(work_dir,groups,replicates);
            close all
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end

    % Finish run
    conclude_SNAP_run(save_analysis_path,SNAP_data,"save",options.save)
end

function feature_data_filtered = filter_SNAP_feature_data(feature_data,groups,replicates)
arguments
    feature_data table
    groups cell
    replicates cell
end
% filter for groups and replicates
fprintf("  Filtering feature data table for desired groups and replicates...\n")
feature_data_filtered = preprocess_SNAP_table(feature_data,...
    "groups",groups,"replicates",replicates,...
    "normalize",false,"remove_NaN",true,...
    "keep_rep_sample_info",true);
% end execution if not enough groups after filtering
if numel(unique(feature_data_filtered.group)) < 2
    ME = MException('SNAP:insufficient_group_number', ...
                    sprintf('Only one group found: %s',groups{1}));
    throw(ME)
end
end

function conclude_SNAP_run(path,data,options)
arguments
    path string
    data struct
    options.save logical = true
end
if options.save
    fprintf("  Saving run to: %s...\n",path)
    starttime_step = tic;
    save_SNAP_run(path,data,"append",false)
    fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
end
    runtime = toc(data.starttime)/60;
    fprintf("RUN END: %s\n", string(datetime))
    fprintf("TIME ELAPSED: %.2f min\n", runtime)
    fprintf("Done!\n")
    fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
    fclose('all');
    warning('on','all');
    close all
    diary off
    return
end

function handle_SNAP_error(ME,filepath,data,options)
arguments
    ME
    filepath string
    data struct
    options.stop_run logical = true
    options.save logical = true
end
    if nargin == 3
        fprintf("- - - - - - - - - - - SNAP ERROR - - - - - - - - - - -\n")
        fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
        fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
    elseif nargin == 4
        fprintf("- - - - - - - - - - - SNAP ERROR - - - - - - - - - - -\n")
        fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
        fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
    end
    if options.save
        fprintf("  Saving run to: %s...\n",path)
        starttime_step = tic;
        save_SNAP_run(filepath,data,"append",false)
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
    if options.stop_run
        conclude_SNAP_run(filepath,data,"save",false);
    end
end

function save_SNAP_run(path,data,options)
arguments
    path string
    data struct
    options.check_overwrite logical = 0
    options.append logical = 1
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
        fprintf("Warning: File does not exist, creating file: "+path)
        save(path,'-struct','data');
    end
% overwrite file
elseif options.check_overwrite
    if exist(path,'file')
        overwrite_dialog();
    else
        save(path,'-struct','data');
    end
% overwrite file with no check
else
    save(path,'-struct','data');
end    
    function overwrite_dialog()
        response = questdlg(sprintf('Warning: %s\nalready exists. Overwrite?',path),...
            'Warning: Overwrite',...
            'Overwrite','New file...','Cancel','Cancel');
        switch response
            case 'Overwrite'
                fprintf('      Overwriting...')
                save(path,'-struct','data');
            case 'New file...'
                new_path = uigetfile();
                save(new_path,'-struct','data');
            case 'Cancel'
                fprintf('Save canceled')
        end
    end
end
