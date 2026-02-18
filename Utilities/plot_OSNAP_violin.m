% -------------------------------------------------------------------------
% plot_OSNAP_violin.m
% -------------------------------------------------------------------------
% Creates a violin plot for each O-SNAP feature across each phenotype
% present in the feature data.
%
% Example on how to use it:
%   plot_OSNAP_violin(feature_data,...
%                     "D:\Analysis\ExperimentA")
% 
% -------------------------------------------------------------------------
% Input:
%   feature_data: The feature data output table, where each row represents 
%                 a sample(nucleus). The first three columns represent (1) 
%                 Group/Phenotype, (2) replicate, and (3) Sample 
%                 Identifier. Each subsequent column is  an O-SNAP feature.
%   work_dir: O-SNAP analysis directory where resulting figures will be 
%             saved in a sub-directory called "violin_plots"
% Options:
%   n_processes: Number of cores to run parallel processes on 
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
function plot_OSNAP_violin(feature_data,work_dir,options)
arguments
    feature_data table
    work_dir string
    options.n_processes double = 12
end
%% Directory management
save_dir = fullfile(work_dir,"violin_plots");
if ~exist(save_dir,"dir")
    mkdir(save_dir)
else
    delete(fullfile(save_dir,'*.png'))
end
%% Get info on groups
group_cat_split = repelem({categorical(replace(feature_data.group,"_"," "))},options.n_processes);
[groups,~,group_idx] = unique(feature_data.group);
n_groups = length(groups);
%% Split for parallel
feature_idx = split_data_to_n_OSNAP(num2cell(4:(size(feature_data,2))),options.n_processes);
T_split = cellfun(@(x) feature_data{:,cell2mat(x)},feature_idx,'uni',0);
feature_split = cellfun(@(x) feature_data.Properties.VariableNames(cell2mat(x)),feature_idx,'uni',0);
%% Create violinplots using multiple workers
parfor p=1:options.n_processes
    feature_idx_p = feature_idx{p};
    c_map = lines(numel(groups));
    for i=1:numel(feature_idx_p)
        feature = feature_split{p}(i);
        feature_data = T_split{p}(:,i);
        group_cat = group_cat_split{p};
        figure('visible','off');
        l = [];
        for g=1:n_groups
            idx = group_idx==g;
            group_idx_g = group_cat(idx);
            feature_data_g = feature_data(idx);
            [feature_data_g_no_outliers, outliers_idx] = rmoutliers(feature_data_g);
            group_idx_g = group_idx_g(~outliers_idx);
            % Violin
            v = violinplot(group_idx_g,feature_data_g_no_outliers);
            v.FaceColor = c_map(g,:);
            hold on
            % Included points
            swarmchart(group_idx_g,feature_data_g_no_outliers,20,c_map(g,:),'filled');
            hold on
            % Median
            median_val = median(feature_data_g,"omitmissing");
            plot([g-0.4 g+0.4], [median_val median_val],"-k","LineWidth",2);
            hold on
            % Mean - NOTE: OUTLIERS ARE NOT REMOVED WHEN CALCULATING MEAN
            scatter(g,mean(feature_data_g,"omitmissing"),50,'pentagramk','filled')
            hold on
            l = [l, string(groups(g)), "", "", ""];
            % Outliers
            if sum(outliers_idx > 0)
                scatter(repmat(g,1,sum(outliers_idx)),feature_data_g(outliers_idx),"xk");
                l = [l ""];
            end
        end
        scatter(nan,nan,50,'pentagramk','filled');
        plot(nan, nan,"-k","LineWidth",2);
        l = [l "Mean" "Median"];
        legend(l,"Location","eastoutside","interpreter","none");
        xlabel("Groups",'FontSize',2)
        xtickangle(30)
        set(gca,'FontSize',16)
        % title(replace(feature,"_"," "),'FontSize',10)
        saveas(gcf,fullfile(save_dir,"violin_plot_"+feature+".png"));
    end
end
