% setup for MrGSEA analysis
function run_NuCLiPSE_MrGSEA(root_dir, analysis_name, T, feature_universe_names, options)
arguments
    root_dir string
    analysis_name string
    T table
    feature_universe_names string = ["feature_universe_1","feature_universe_2","feature_universe_3"];
    options.alpha double = 0.05;
    options.plot logical = true;
end
groups = unique(T.group);
group_pairs = nchoosek(groups,2);
for g=1:size(group_pairs,1)
    group_ctrl = group_pairs(g,1);
    group_case = group_pairs(g,2);
    %% directory management
    work_dir = fullfile(root_dir,analysis_name);
    for i=1:length(feature_universe_names)
        current_dir = cd(work_dir);
        %% file management
        feature_universe_name = feature_universe_names(i);
        save_name = char(fullfile(work_dir,analysis_name + "_GSEA_" + feature_universe_name + "_" + group_ctrl + "_" + group_case));
        plot_dir = "GSEA_plots_"+feature_universe_name + "_" + group_ctrl + "_" + group_case;
        if exist(save_name+".mat",'file') && exist(plot_dir,'dir') && ~options.plot
            load(save_name,"T_GSEA")
            display_GSEA_results(T_GSEA,group_pairs(g,:),"alpha",options.alpha);
            continue
        end
        %% prepare data
        % create data array (pt 1)
        T_norm = prepare_voronoi_table_data(T,"groups",cellstr([group_ctrl group_case]),"replicates",cellstr(unique(T.biological_replicate)),"normalize",true);
        % get feature names (feature_IDs)
        feature_IDs = T_norm.Properties.VariableNames(4:end);
        % define groups
        group_IDs = strcmpi(T_norm.group,group_ctrl);
        % create data array (pt 2)
        X = table2array(T_norm(:,4:end))';
        % get gene set data collection for appropriate feature universes and feature IDs (GS)
        feature_universe = get_feature_universes(feature_universe_names{i});
        GS_name = ['NuCLiPSE_' char(feature_universe_name) '_GS'];
        file_path_GS = [GS_name '.xlsx'];
        if ~exist(file_path_GS,"file")
            generate_NuCLiPSE_GS(feature_IDs',feature_universe{1},"filepath",file_path_GS)
        end
        %% set options
        opts = default_GSEA_opts();
        opts.show = true;       %if plot results
        opts.save =  false;      %if save results
        opts.perm_nb = 1000;    %number of permutations
        % Default parameters for GSEA method.
        opts.GS_name = GS_name; %GeneSet database name
        opts.GS_filt = [1,500]; %minimum and maximum number of genes in GS
        opts.sort_type = 'descend'; %type of ranks sorting
        opts.p = 1; %0 - Kolmogorov-Smirnov statistic, 1-weighting genes by ranking metric
        opts.rank_type = 'S2N';%ranking method:
        % {'S2N','ttest','ratio','diff','log2_ratio','SoR','BWS','ReliefF','WAD',...
        % 'FCROS','MWT','MSD','external'} (dedicated for own rank implementations)
        opts.abs = true; %if use absolute value of ranking statistic
        opts.tied_rank = false; %if create tied ranks of ranking statistic
        opts.perm_type = 'pheno'; %type of permutation {'entrez','pheno'}
        opts.perm_nb = 1000; %number of permutations
        opts.perm_ext = false; %if extend number of permutation till one significant result occur(may be time consuming)
        %% run MrGSEA
        % POS – results for up-regulated enriched gene sets (Positive value of Enrichment Score),
        % NEG – results for down-regulated enriched gene sets (Negative value of Enrichment Score),
        % DESCR – variable with description of each column in POS and NEG matrices.,
        % GenePval – P-value for each gene (useful for ranking metrics where permutations are needed to assess significance of expression differentiation).
        [res_pos,res_neg,res_descr,~] = MrGSEA(X,group_IDs,1:length(feature_IDs),char(save_name),opts);
        %% create output
        T_GSEA = cell2table([res_pos; res_neg],'VariableNames',res_descr(1:9));
        save([save_name '.mat'], "T_GSEA");
        display_GSEA_results(T_GSEA,group_pairs(g,:),"alpha",options.alpha);
        %% plot
        plot_GSEA(work_dir,T_GSEA,feature_universe_name,"alpha",options.alpha);
        %% wrap up
        movefile(fullfile(work_dir,"GSEA_plots","GSEA_plots_"+feature_universe_name +".fig"),fullfile(work_dir,"GSEA_plots","GSEA_plots_"+feature_universe_name + "_" + group_ctrl + "_" + group_case+".fig"));
        movefile(fullfile(work_dir,"GSEA_plots","GSEA_plots_"+feature_universe_name +".png"),fullfile(work_dir,"GSEA_plots","GSEA_plots_"+feature_universe_name + "_" + group_ctrl + "_" + group_case+".png"));
        movefile("GSEA_plots",plot_dir);
        %% wrap up
        close all
        cd(current_dir)
    end
end
end

function display_GSEA_results(T_GSEA, group_pair, options)
arguments
    T_GSEA table
    group_pair string
    options.alpha double = 0.05
end
    families_significant = T_GSEA.("NES_q-val") < options.alpha;
    if any(families_significant)
        disp(group_pair(1) + " VS " + group_pair(2))
        disp("- - - - - - - - - - - - - - - - -")
        T_GSEA_filt = sortrows(T_GSEA(families_significant,:),"NES_q-val");
        for j=1:sum(families_significant)
            fprintf("%s:\t%4.2f\t%4.2e\n\n",string(T_GSEA_filt{j,"GS.NAME"}),T_GSEA_filt{j,"NES"},T_GSEA_filt{j,"NES_q-val"});
        end
        disp("- - - - - - - - - - - - - - - - -")
    end
end
