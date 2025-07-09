% -------------------------------------------------------------------------
% run_SNAP.m
% -------------------------------------------------------------------------
% Function that runs the suite of analysis available in O-SNAP, excluding
% the pseudotime analysis that implemented in R.
% <FILL>
% Most descriptors/features that are calculated are calculated with respect to
% the point cloud or an accurate representation of a pointcloud (based on
% the alphaShape or polyshape functions of Matlab). 
%
% Example on how to use it:
%   results = run_SNAP("C:\Users\JohnSmith\Documents\O-SNAP_Analysis\,...
%                      "SmithJohn_HeLa_ProteinKO_H2BHistone",...
%                      {'Control','KO'},...
%                      {'20250101','20250108','20250201'},...
%                      "run_generate_features",1,
%                      "run_generate_table",1,...
%                      "run_comparison",1,...
%                      "run_generate_batches",1,
%                      "run_feature_selection",1,
%                      "run_PCA",1, ...
%                      "run_classification_batch",1, ...
%                      "run_plot_violin",1,...
%                      "run_plot_radial",1,...
%                      "run_FSEA",1,...
%                      "save",1);
% -------------------------------------------------------------------------
% Input:
%   root_dir:  The parent directory that contains the directory where the
%              analysis results are saved
%   analysis_name: Analysis identifier string
%   groups: Cell array containing char array of the identifiers of the 
%           phenotypes/cell states *ENSURE THAT EACH SAMPLE FILENAME 
%           CONTAINS EXACTLY ONE IDENTIFIER FROM groups*
%   replicates: Cell array containing char array of the identifiers of the 
%               replicates. *ENSURE THAT EACH REPLCIATE IS STORED IN A
%               SEPARATE FOLDER PREFIXED BY "Rep"*
% Output:
%   SNAP_data:  A struct containing the following fields related to the 
%               SNAP analysis:
%                  - analysis_name: Analysis identifier string
%                  - classification_summary_all: Table of classification
%                                                accuracies using "all"
%                                                approach
%                  - classification_summary_batch_all: Table of
%                                                      classification
%                                                      accuracies using
%                                                      "batch_all" approach
%                  - classification_summary_batch_each: Table of
%                                                       classification 
%                                                       accuracies using 
%                                                       "batch_each"
%                                                       approach
%                  - classifiers_all: Cell array with all instances of 
%                                     trained models using "all" approach
%                  - classifiers_batch_all: Cell array with all instances
%                                           of trained models using
%                                           "batch_all" approach
%                  - classifiers_batch_each: Cell array with all instances
%                                            of trained models using
%                                            "batch_each" approach
%                  - date: The most recent date the analysis was run
%                  - feature_comparisons: A cell array where every cell is
%                                         a pair-wise combination of the
%                                         fold-change analysis used in the
%                                         volcano plots
%                  - feature_data: The table containing the O-SNAP feature
%                                  values where each row represents a
%                                  sample (nucleus) and each column is an
%                                  O-SNAP feature
%                  - groups: Cell array containing char array of the
%                            identifiers of the phenotypes/cell states
%                            *ENSURE THAT EACH SAMPLE FILENAME CONTAINS 
%                            EXACTLY ONE IDENTIFIER FROM groups*
%                  - options: Struct array of parameter values for 
%                             O-SNAP analysis
%                  - pca_result_all: Info on PCA transformation in the 
%                                    classification pipeline using the
%                                    "all" approach
%                  - pca_result_batch_all: Info on the PCA transformation
%                                          in the classification pipeline
%                                          using the "batch_all" approach
%                  - pca_result_batch_each: Info on the PCA transformation
%                                          in the classification pipeline
%                                          using the "batch_each" approach"
%                  - replicates: Cell array containing char array of the 
%                                identifiers of the replicates.
%                                *ENSURE THAT EACH REPLCIATE IS STORED IN
%                                 A SEPARATE FOLDER PREFIXED BY "Rep"*
%                  - starttime: The starttime identifier for the analysis
%                               run
%                  - test_idxs: A cell array where each cell contains a
%                               logical array indicating samples for the
%                               test data of each fold
%                  - train_idxs: A cell array where each cell contains a
%                                logical array indicating samples for the
%                                training data of each fold
%                  - vars_select_result_all: MRMR scores of features for
%                                            ranking performed on entire
%                                            feature data
%                  - vars_select_result_batch: Table with MRMR scores of
%                                              each fold and the sum of the
%                                              scores across all folds
%                  - vars_selected_all: Features selected from MRMR
%                                       performed on entire dataset
%                  - vars_selected_batch: Table of boolean values
%                                         indicating whether a feature is
%                                         selected within a batch and from
%                                         the aggregated selection
%                  - venn_data: Data for Venn diagram to compare changes
%                               between 3-4 comparisons of phenotype pairs
% Options:
%   run_generate_features:  
%   run_generate_table:  
%   run_comparison:  
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
function SNAP_data = run_SNAP(root_dir,analysis_name,groups,replicates,options)
    arguments
        root_dir char
        analysis_name char
        groups cell
        replicates cell
        % run options (will overwrite)
        options.run_generate_features logical = 1
        options.run_generate_table logical = 1
        options.run_comparison logical = 1
        options.run_generate_batches logical = 1
        options.run_feature_selection logical = 1
        options.run_PCA logical = 1
        options.run_classification_all logical = 1
        options.run_classification_batch logical = 1
        options.run_venn logical = 0
        options.run_plot_violin logical = 0
        options.run_plot_radial logical = 0
        options.run_FSEA logical = 0
        % generate batch options
        options.split_method string = "k-fold"
        options.test_train_ratio double = 0.2;
        options.test_train_k double = 5;
        % feature select options
        options.max_idx = 12;
        % PCA options
        options.num_components_explained double = 0.75
        % classification options
        options.n_models_per_type = 5;
        % comparison options
        options.alpha double = 0.05
        options.fold_change_threshold double = 2;
        % FSEA options
        options.feature_universe_names = ["universe_1"];
        options.FSEA_rank_type string = "S2N";
        % other
        options.n_processes double = 12
        % save options
        options.suffix = ""
        options.save = 1;
        options.save_if_error = 1;
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
    close('all');

    rng('default');

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
            % generate SNAP feature table data if not loaded
            if isfield(SNAP_data,'feature_data')
                feature_data = SNAP_data.feature_data;
            else
                options.run_generate_table = 1;
            end
            % run feature comparisons if not loaded
            if ~isfield(SNAP_data,'feature_comparisons')
                options.run_comparison = 1;
            end
            % run steps prior to classification
            if options.run_classification_batch
                if options.run_generate_features
                    options.run_generate_table = 1;
                    options.run_generate_batches = 1;
                    options.run_feature_selection = 1;
                    options.run_PCA = 1;
                end
                if ~all(isfield(SNAP_data,{'train_idxs','test_idxs'}))
                    options.run_generate_batches = 1;
                elseif isempty(SNAP_data.train_idxs) || isempty(SNAP_data.test_idxs)
                    options.run_generate_batches = 1;
                end
                if ~all(isfield(SNAP_data,{'vars_select_result_batch','vars_selected_batch'}))
                    options.run_feature_selection = 1;
                elseif isempty(SNAP_data.vars_select_result_batch) || isempty(SNAP_data.vars_selected_batch)
                    options.run_feature_selection = 1;
                end
                if ~isfield(SNAP_data,'pca_result_batch_each')
                    options.run_PCA = 1;
                elseif isempty(SNAP_data.pca_result_batch_each)
                    options.run_PCA = 1;
                end
                if ~isfield(SNAP_data,'pca_result_batch_all')
                    options.run_PCA = 1;
                elseif isempty(SNAP_data.pca_result_batch_all)
                    options.run_PCA = 1;
                end
            end
            SNAP_data.options = options;
        end
    catch ME
        options.run_generate_table = 1;
        options.run_comparison = 1;
        options.run_generate_batches = 1;
        options.run_feature_selection = 1;
        options.run_PCA = 1;
        options.run_classification_batch = 1;
        SNAP_data.options = options;
        if exist('SNAP_data','var')
            handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        end
    end

    SNAP_data.date = datetime;
    SNAP_data.starttime = tic;

    %% generate SNAP features per nucleus
    try
        % generate_features flag takes precedence over load_features flag
        if options.run_generate_features
            fprintf("  Generating features from scratch...\n")
            % calculate features
            starttime_step = tic;
            generate_SNAP_features(work_dir,groups,replicates,'overwrite',true,'n_processes',options.n_processes);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        fprintf("%s\n",getReport(ME));
        return
    end

    %% generate SNAP feature table data
    try
        if options.run_generate_table
            fprintf("  Creating table from generated features...\n");
            % coallate features
            starttime_step = tic;
            feature_data = SNAP_nucleus_to_feature_table_batch(work_dir,groups,replicates);
            SNAP_data.feature_data = feature_data;
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end

    %% filter data
    if exist("feature_data","var")
        try
            feature_data_filtered = filter_SNAP_feature_data(feature_data,groups,replicates);
            writetable(feature_data_filtered,replace(save_analysis_path,".mat",".csv"));
        catch ME
            handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
            return
        end

    else
        fprintf("No feature data, ending run...\n");
        return
    end

    %% violin plots
    try
        if options.run_plot_violin
            fprintf("  Plotting violin plots...\n")
            starttime_step = tic;
            plot_SNAP_violin(feature_data_filtered,work_dir)
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end

    %% compare all groups to each other
    try
        if options.run_comparison || (options.run_venn && ~isfield(SNAP_data,"feature_comparisons"))
            fprintf("  Comparing features...\n")
            starttime_step = tic;
            save_path = fullfile(work_dir, "features_volcano");
            if options.suffix ~= ""
                save_path = join(save_path,options.suffix,"_");
            end
            SNAP_data.feature_comparisons = compare_groups(feature_data_filtered,...
                "alpha",options.alpha,...
                "fold_change_threshold",options.fold_change_threshold,...
                "save_path",save_path);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    if isfield(SNAP_data,"feature_comparisons")
        for i=1:numel(SNAP_data.feature_comparisons)
            disp(join(SNAP_data.feature_comparisons{i}.groups," vs "))
            fprintf("UP: %3.0f\n",sum(all([SNAP_data.feature_comparisons{i}.feature_table{:,5}>1, SNAP_data.feature_comparisons{i}.feature_table{:,7} < 0.05],2)))
            fprintf("DOWN: %3.0f\n",sum(all([SNAP_data.feature_comparisons{i}.feature_table{:,5}<-1, SNAP_data.feature_comparisons{i}.feature_table{:,7} < 0.05],2)))
        end
    end
    %% venn diagram
    try
        if options.run_venn
            fprintf("  Creating venn diagrams...\n")
            starttime_step = tic;
            SNAP_data.venn_data = interpret_features_venn(SNAP_data.feature_comparisons,"save_path",fullfile(work_dir, "features_venn"));
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    %% run classification steps if needed
    % generate train and test batches
    try
        if options.run_generate_batches
            fprintf("  Creating train/test batches...\n")
            starttime_step = tic;
            [SNAP_data.train_idxs, SNAP_data.test_idxs,by] = generate_train_test_sets(...
                feature_data_filtered,...
                "split_method",options.split_method,...
                "k",options.test_train_k,...
                "proportion",options.test_train_ratio);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
            options.split_method = by;
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    % feature selection
    try
        if options.run_feature_selection
            % clear variable for new run
            if isfield(SNAP_data,'vars_select_result_batch')
                SNAP_data = rmfield(SNAP_data,"vars_select_result_batch");
            end
            fprintf("  Selecting features (batch-wise)...\n")
            starttime_step = tic;
            save_path = fullfile(work_dir, join(['feature_selection_batch',groups],'_'));
            % split batches for parallel
            train_idx_p = split_data_to_n(SNAP_data.train_idxs,options.n_processes,"shuffle",false);
            n_processes = numel(train_idx_p);
            feature_data_p = cell(1,n_processes);
            var_scores_p = cell(1,n_processes);
            for p=1:n_processes
                feature_data_p{p} = feature_data_filtered;
            end
            parfor p=1:n_processes
                train_idx_b = train_idx_p{p};
                for b=1:numel(train_idx_b)
                    train_data = feature_data_p{p}(train_idx_b{b},:);
                    [~, var_scores_p{p}{b}] = fscmrmr(train_data(:,vartype('numeric')),train_data.group);
                end
            end
            % bring together batches from each parallel process
            scores = [var_scores_p{:}];
            scores = vertcat(scores{:})';
            % store values
            [SNAP_data.vars_select_result_batch,SNAP_data.vars_selected_batch] = select_features_SNAP( ...
                scores, ...
                feature_data_filtered(:,vartype('numeric')).Properties.VariableNames, ...
                "max_idx",options.max_idx, ...
                "save_path",save_path);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end
    if isfield(SNAP_data,'vars_selected_batch')
        vars_display = SNAP_data.vars_selected_batch.Properties.RowNames(logical(SNAP_data.vars_selected_batch.batch_sum));
        fprintf("      Selected features (n = %.0f):\n",numel(vars_display))
        for i=1:numel(vars_display)
            fprintf("        - %s\n", vars_display{i})
        end
    end
    %% PCA
    try
        if options.run_PCA
            save_path = fullfile(work_dir, string(join(['PCA',groups],'_')));
            if options.suffix ~= ""
                save_path = string(join([save_path options.suffix],'_'));
            end
            fprintf("  Performing PCA (batch-wise)...\n")
            starttime_step = tic;
            % split batches for parallel
            num_components_explained = options.num_components_explained;
            train_idx_p = split_data_to_n(SNAP_data.train_idxs,options.n_processes,"shuffle",false);
            n_processes = numel(train_idx_p);
            n_batch = numel(SNAP_data.train_idxs);
            if ~exist('feature_data_p','var')
                feature_data_p = cell(1,n_processes);
                for p=1:n_processes
                    feature_data_p{p} = feature_data_filtered;
                end
            end
            % pca based on feature selection by EACH batch
            if isfield(SNAP_data,'vars_select_result_batch') 
                vars_selected_b = split_data_to_n(arrayfun(@(x) SNAP_data.vars_selected_batch.Properties.RowNames(logical(SNAP_data.vars_selected_batch{:,x})), 1:n_batch,'uni',0),n_processes,"shuffle",false);
            else
                vars_selected_b = split_data_to_n(cell(1,n_batch),n_processes,"shuffle",false);
            end
            pca_result_b = cell(1,n_processes);
            parfor p=1:n_processes
                train_idx_b = train_idx_p{p};
                for b=1:numel(train_idx_b)
                    train_data = feature_data_p{p}(train_idx_b{b},:);
                    pca_result_b{p}{b} = run_SNAP_PCA(train_data,...
                        "vars_sel",vars_selected_b{p}{b},...
                        "save_path",sprintf("%s_%02.0f_%02.0f",save_path,p,b),...
                        "num_components_explained",num_components_explained);
                end
            end
            SNAP_data.pca_result_batch_each = [pca_result_b{:}];
            % pca based on feature selection on ALL batches
            if isfield(SNAP_data,'vars_select_result_batch')
                vars_selected_b = split_data_to_n(repmat({SNAP_data.vars_selected_batch.Properties.RowNames(logical(SNAP_data.vars_selected_batch{:,end}))},1,n_batch),n_processes,"shuffle",false);
            else
                vars_selected_b = split_data_to_n(cell(1,n_batch),n_processes,"shuffle",false);
            end
            pca_result_b = cell(1,n_processes);
            parfor p=1:n_processes
                train_idx_b = train_idx_p{p};
                for b=1:numel(train_idx_b)
                    train_data = feature_data_p{p}(train_idx_b{b},:);
                    pca_result_b{p}{b} = run_SNAP_PCA(train_data,...
                        "vars_sel",vars_selected_b{p}{b},...
                        "save_path",save_path+"_all",...
                        "num_components_explained",num_components_explained);
                end
            end
            SNAP_data.pca_result_batch_all = [pca_result_b{:}];
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
    end
    %% run classification batch-wise
    try
        if options.run_classification_batch
            n_models_per_type = options.n_models_per_type;
            fprintf("  Running classification (batch-wise)...\n")
            starttime_step = tic;
            % cell_dims = [numel(SNAP_data.train_idxs),numel(SNAP_data.vars_selected_batch),numel(SNAP_data.pca_result)];
            % assert(all(cell_dims == cell_dims(1)),'SNAP:classification_length_mismatch','number of batches not consistent with classification structures (feature selection, PCA)')
            % split batches for parallel
            n_batch = numel(SNAP_data.train_idxs);
            % classification based on feature selection by EACH batch
            if isfield(SNAP_data,'vars_select_result_batch')
                vars_selected_b = arrayfun(@(x) SNAP_data.vars_selected_batch.Properties.RowNames(logical(SNAP_data.vars_selected_batch{:,x}))',1:n_batch,'uni',0);
            else
                vars_selected_b = cell(1,n_batch);
            end
            if isfield(SNAP_data,'pca_result_batch_each')
                if ~isempty(SNAP_data.pca_result_batch_each)
                    pca_result_b = SNAP_data.pca_result_batch_each;
                else
                    pca_result_b = cell(1,n_batch);
                end
            else
                pca_result_b = cell(1,n_batch);
            end
            classifiers = cell(1,numel(SNAP_data.train_idxs));
            for b=1:n_batch
                train_data = feature_data_filtered(SNAP_data.train_idxs{b},:);
                test_data = feature_data_filtered(SNAP_data.test_idxs{b},:);
                classifiers{b} = run_SNAP_classification_batch(...
                    train_data,...
                    'n_models_per_type',n_models_per_type,...
                    'vars_selected',vars_selected_b{b},...
                    'pca_result',pca_result_b{b},....
                    'test_data',test_data,...
                    'verbose',0);
            end
            SNAP_data.classifiers_batch_each = classifiers;
            % plot figures
            if options.save
                SNAP_data.classification_summary_batch_each = evaluate_classifiers(SNAP_data.classifiers_batch_each,"save_path",fullfile(work_dir,"batch_each"));
            else
                SNAP_data.classification_summary_batch_each = evaluate_classifiers(SNAP_data.classifiers_batch_each);
            end
            
            % classification based on feature selection on ALL batches
            if isfield(SNAP_data,'vars_select_result_batch')
                vars_selected_b = SNAP_data.vars_selected_batch.Properties.RowNames(logical(SNAP_data.vars_selected_batch{:,end}))';
            else
                vars_selected_b = [];
            end
            if isfield(SNAP_data,'pca_result_batch_all')
                if ~isempty(SNAP_data.pca_result_batch_all)
                    pca_result_b = SNAP_data.pca_result_batch_all;
                else
                    pca_result_b = cell(1,n_batch);
                end
            else
                pca_result_b = cell(1,n_batch);
            end
            classifiers = cell(1,n_batch);
            for b=1:n_batch
                train_data = feature_data_filtered(SNAP_data.train_idxs{b},:);
                test_data = feature_data_filtered(SNAP_data.test_idxs{b},:);
                classifiers{b} =...
                    run_SNAP_classification_batch(...
                    train_data,...
                    'n_models_per_type',n_models_per_type,...
                    'vars_selected',vars_selected_b,...
                    'pca_result',pca_result_b{b},....
                    'test_data',test_data,...
                    'verbose',0);
            end
            SNAP_data.classifiers_batch_all = classifiers;
            if options.save
                SNAP_data.classification_summary_batch_all = evaluate_classifiers(SNAP_data.classifiers_batch_all,"save_path",fullfile(work_dir,"batch_all"));
            else
                SNAP_data.classification_summary_batch_all = evaluate_classifiers(SNAP_data.classifiers_batch_allSNAP_data.classifiers_batch_all);
            end
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
    end
    if isfield(SNAP_data,"classification_summary_batch_each")
        fprintf('    Classification summary (batch - fully independent feature selection):\n');
        head(SNAP_data.classification_summary_batch_each,3)
    end
    if isfield(SNAP_data,"classification_summary_batch_all")
        fprintf('    Classification summary (batch - aggregated feature selection):\n');
        head(SNAP_data.classification_summary_batch_all,3)
    end

    %% run classification on whole data
    if options.run_classification_all
        n_models_per_type = options.n_models_per_type;
        fprintf("  Running classification (whole dataset)...\n")
        starttime_step = tic;
        % cell_dims = [numel(SNAP_data.train_idxs),numel(SNAP_data.vars_selected_batch),numel(SNAP_data.pca_result)];
        % assert(all(cell_dims == cell_dims(1)),'SNAP:classification_length_mismatch','number of batches not consistent with classification structures (feature selection, PCA)')
        % split batches for parallel
        % variable selection
        save_path = fullfile(work_dir, join(['feature_selection_all',groups],'_'));
        if options.suffix ~= ""
            save_path = string(join([save_path options.suffix],'_'));
        end
        [~, SNAP_data.vars_select_result_all] = fscmrmr(feature_data_filtered(:,vartype('numeric')),feature_data_filtered.group);
        [~,SNAP_data.vars_selected_all] = select_features_SNAP(SNAP_data.vars_select_result_all',feature_data_filtered(:,vartype('numeric')).Properties.VariableNames,"max_idx",options.max_idx,"save_path",save_path);
        fprintf("      Selected features - ALL (n = %.0f):\n",numel(SNAP_data.vars_selected_all))
        for i=1:numel(SNAP_data.vars_selected_all)
            fprintf("        - %s\n", SNAP_data.vars_selected_all{i})
        end
        % pca
        save_path = fullfile(work_dir, string(join(['PCA',groups],'_')));
        if options.suffix ~= ""
            save_path = string(join([save_path options.suffix],'_'));
        end
        SNAP_data.pca_result_all = run_SNAP_PCA(feature_data_filtered,...
                                        "vars_sel",SNAP_data.vars_selected_all,...
                                        "save_path",save_path,...
                                        "num_components_explained",options.num_components_explained);
        % classification
        classifiers =...
            run_SNAP_classification_batch(...
                feature_data_filtered,...
                'n_models_per_type',n_models_per_type,...
                'vars_selected',SNAP_data.vars_selected_all,...
                'pca_result',SNAP_data.pca_result_all,....
                'verbose',0);
        SNAP_data.classifiers_all = {classifiers};
        if options.save
            SNAP_data.classification_summary_all = evaluate_classifiers(SNAP_data.classifiers_all,"save_path",fullfile(work_dir,"all"));
        else
            SNAP_data.classification_summary_all = evaluate_classifiers(SNAP_data.classifiers_allSNAP_data.classifiers_all);
        end
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
    if isfield(SNAP_data,"classification_summary_all")
        fprintf('    Classification summary (batch - aggregated feature selection):\n');
        head(SNAP_data.classification_summary_all,3)
    end
    %% plot radial
    try
        if options.run_plot_radial
            fprintf("  Plotting radial densities...\n")
            plot_ellipse_set_batch(feature_data_filtered,work_dir,groups,replicates);
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
            calculate_feature_set_coverage(SNAP_data.feature_comparisons,work_dir,"plot",false,"feature_universe_names",options.feature_universe_names);
            run_SNAP_MrFSEA(work_dir, feature_data_filtered, options.feature_universe_names, "FSEA_rank_type",options.FSEA_rank_type);
            fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        end
    catch ME
        handle_SNAP_error(ME,save_analysis_path,SNAP_data,"save",options.save_if_error);
        return
    end

    % Finish run
    SNAP_data.options = options;
    conclude_SNAP_run(save_analysis_path,SNAP_data,"save",options.save)
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
    set(findall(groot,'Type','Figure'),'visible','on');
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