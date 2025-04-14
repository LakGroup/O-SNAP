% setup for MrFSEA analysis
function run_SNAP_MrFSEA(work_dir, T, feature_universe_names, options)
arguments
    work_dir string
    T table
    feature_universe_names string = ["universe_1"];
    options.alpha double = 0.05;
    options.plot logical = true;
    options.FSEA_rank_type string = "S2N";
end
groups = unique(T.group);
group_pairs = nchoosek(groups,2);
for g=1:size(group_pairs,1)
    group_ctrl = group_pairs(g,1);
    group_case = group_pairs(g,2);
    %% directory management
    for i=1:length(feature_universe_names)
        current_dir = cd(work_dir);
        %% file management
        feature_universe_name = feature_universe_names(i);
        save_dir = "FSEA_plots_"+feature_universe_name + "_" + group_ctrl + "_" + group_case;
        save_path = char(fullfile(work_dir,"FSEA_" + feature_universe_name + "_" + group_ctrl + "_" + group_case));
        % if exist(save_path+".mat",'file') && exist(save_dir,'dir') && ~options.plot
            % load(save_path,"T_FSEA")
            % display_FSEA_results(T_FSEA,group_pairs(g,:),"alpha",options.alpha);
            % continue
        % elseif ~exist("save_dir","dir")
            mkdir(save_dir)
        % end
        %% prepare data
        % create data array (pt 1)
        T_norm = preprocess_SNAP_table(T,"groups",cellstr([group_ctrl group_case]),"replicates",cellstr(unique(T.biological_replicate)),"normalize",true);
        % get feature names (feature_IDs)
        feature_IDs = T_norm.Properties.VariableNames(2:end);
        % define groups
        group_IDs = strcmpi(T_norm.group,group_case);
        % create data array (pt 2)
        T_norm = T_norm{:,2:end}';
        % get gene set data collection for appropriate feature universes and feature IDs (FS)
        feature_universe = get_feature_universes(feature_universe_names(i));
        file_path_feature_set = feature_universe_name;
        % if ~exist(file_path_feature_set+ ".xlsx","file")
            generate_SNAP_feature_set(feature_IDs',feature_universe{i})
        % end
        %% set options
        opts = default_FSEA_opts();
        opts.show = true;       %if plot results
        opts.save =  false;      %if save results
        opts.perm_nb = 1000;    %number of permutations
        % Default parameters for FSEA method.
        opts.FS_name = file_path_feature_set; %FeatureSet database name
        opts.FS_filt = [1,500]; %minimum and maximum number of genes in FS
        opts.sort_type = 'descend'; %type of ranks sorting
        opts.p = 1; %0 - Kolmogorov-Smirnov statistic, 1-weighting genes by ranking metric
        opts.rank_type = options.FSEA_rank_type;%ranking method:  {'S2N','ttest','ratio','diff',...
        % 'log2_ratio','WAD','FCROS','MWT','external'}...
        % (dedicated for own rank implementations)
        opts.abs = false; %if use absolute value of ranking statistic
        opts.tied_rank = false; %if create tied ranks of ranking statistic
        opts.perm_type = 'pheno'; %type of permutation {'entrez','pheno'}
        opts.perm_nb = 1000; %number of permutations
        opts.perm_ext = false; %if extend number of permutation till one significant result occur(may be time consuming)
        %% run MrFSEA
        % POS – results for increased value feature sets (Positive value of Enrichment Score),
        % NEG – results for decreased value feature sets (Negative value of Enrichment Score),
        % DESCR – variable with description of each column in POS and NEG matrices.,
        % FeaturePval – P-value for each feature (useful for ranking metrics where permutations are needed to assess significance of expression differentiation).
        [res_pos,res_neg,res_descr,~] = MrFSEA(work_dir,T_norm,group_IDs,1:length(feature_IDs),char(save_path),opts);
        %% create output
        T_FSEA = cell2table([res_pos; res_neg],'VariableNames',res_descr(1:9));
        save([save_path '.mat'], "T_FSEA");
        display_FSEA_results(T_FSEA,group_pairs(g,:),"alpha",options.alpha);
        %% plot
        plot_SNAP_FSEA(T_FSEA,group_ctrl,group_case,fullfile(save_dir,"FSEA_plots_"+feature_universe_name),"alpha",options.alpha);
        %% wrap up
        close all
        cd(current_dir)
    end
end
end