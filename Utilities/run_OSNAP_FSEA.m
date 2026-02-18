% -------------------------------------------------------------------------
% run_OSNAP_FSEA.m
% -------------------------------------------------------------------------
% Runs the Feature Set Enrichment analysis, finding information on how a
% feature universe is defined and then determining how feature sets change
% between all combination of group pairs.
%
% Example on how to use it:
%   run_OSNAP_FSEA("D:\Analysis\ExperimentA",...
%                   OSNAP_feature_data,... 
%                   ["universe_1"])
% -------------------------------------------------------------------------
% Input:
%   work_dir: O-SNAP analysis directory. Should contain one or more folders
%             with the MAT O-SNAP feature files for each sample (nucleus),
%             divided by replicates with each folder prefixed by "Rep"
%   feature_data: The table containing the O-SNAP feature
%                 values where each row represents a sample (nucleus) and 
%                 each column is an O-SNAP feature
%   feature_universe_names: Identifier for the feature universe to analyze.
%                           Must be a valid universe defined in 
%                           "get_OSNAP_feature_universes.m"
% Options:
%   alpha: The p-value cutoff for significance testing
%   plot: Flag on whether to plot the FSEA results
%   FSEA_rank_type: Method to rank features in the FSEA algorithm
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
function run_OSNAP_FSEA(work_dir, feature_data, feature_universe_names, options)
arguments
    work_dir string
    feature_data table
    feature_universe_names string = ["universe_1"];
    options.alpha double = 0.05;
    options.plot logical = true;
    options.FSEA_rank_type string = "S2N";
end
groups = unique(feature_data.group);
group_pairs = nchoosek(groups,2);
for g=1:size(group_pairs,1)
    group_ctrl = group_pairs(g,1);
    group_case = group_pairs(g,2);
    for i=1:length(feature_universe_names)
        current_dir = cd(work_dir);
        %% Directory management
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
        %% Prepare feature set information
        % 1. Preprocess feature data array
        T_norm = preprocess_OSNAP_feature_data(feature_data,"groups",cellstr([group_ctrl group_case]),"replicates",cellstr(unique(feature_data.biological_replicate)),"normalize",true);
        % 2. Get feature names (feature_IDs)
        feature_IDs = T_norm.Properties.VariableNames(2:end);
        % 3. Define groups
        group_IDs = strcmpi(T_norm.group,group_case);
        % 4. Create data array
        T_norm = T_norm{:,2:end}';
        % 5. Get feature set data collection for appropriate feature universes and feature IDs (FS)
        feature_universe = get_OSNAP_feature_universes(feature_universe_names(i));
        file_path_feature_set = feature_universe_name;
        % 6. Create feature set from the feature universe
        % if ~exist(file_path_feature_set+ ".xlsx","file")
        generate_OSNAP_feature_set(feature_IDs',feature_universe{i})
        % end
        %% Set options
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
        %% Run MrFSEA
        % POS – results for increased value feature sets (Positive value of Enrichment Score),
        % NEG – results for decreased value feature sets (Negative value of Enrichment Score),
        % DESCR – variable with description of each column in POS and NEG matrices.,
        % FeaturePval – P-value for each feature (useful for ranking metrics where permutations are needed to assess significance of expression differentiation).
        [res_pos,res_neg,res_descr,~] = MrFSEA(work_dir,T_norm,group_IDs,1:length(feature_IDs),char(save_path),opts);
        %% Create output
        T_FSEA = cell2table([res_pos; res_neg],'VariableNames',res_descr(1:9));
        save([save_path '.mat'], "T_FSEA");
        display_FSEA_results(T_FSEA,group_pairs(g,:),"alpha",options.alpha);
        %% Plot
        plot_OSNAP_FSEA(T_FSEA,group_ctrl,group_case,fullfile(save_dir,"FSEA_plots_"+feature_universe_name),"alpha",options.alpha);
        %% Wrap up
        close all
        cd(current_dir)
    end
end
end
