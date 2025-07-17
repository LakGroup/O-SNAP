% -------------------------------------------------------------------------
% plot_OSNAP_ellipses_batch.m
% -------------------------------------------------------------------------
% Performs the ellipse plot rendering in parallel for multiple samples
%
% Example on how to use it:
%   plot_OSNAP_ellipses_batch(feature_data,...
%                             "D:\Analysis\ExperimentA",...
%                             {'Control','KO'},...
%                             {'20250101','20250108','20250201'})
% -------------------------------------------------------------------------
% Input:
%   feature_data: The feature data output table, where each row represents 
%                 a sample(nucleus). The first three columns represent (1) 
%                 Group/Phenotype, (2) replicate, and (3) Sample 
%                 Identifier. Each subsequent column is  an O-SNAP feature. 
%   work_dir: O-SNAP analysis directory. Should contain one or more folders
%             with the MAT O_SNAP feature files for each sample (nucleus),
%             divided by replicates with each folder prefixed by "Rep"
%   groups: Cell array containing char array of the identifiers of the 
%           phenotypes/cell states
%           *ENSURE THAT EACH SAMPLE FILENAME CONTAINS EXACTLY ONE 
%            IDENTIFIER FROM groups*
%   replicates: Cell array containing char array of the identifiers of the 
%               replicates.
%               *ENSURE THAT EACH REPLCIATE IS STORED IN A SEPARATE FOLDER 
%                PREFIXED BY "Rep"* 
% Options:
%   n_processes: Number of cores to run parallel processes on 
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function plot_OSNAP_ellipses_batch(feature_data,work_dir,groups,replicates,options)
arguments
    feature_data table
    work_dir string
    groups cell
    replicates cell
    options.n_processes double= 12;
end
%% get data
data_info_table = get_valid_OSNAP_samples(work_dir,groups,replicates,{'radial_loc_density','radial_dbscan_cluster_density'});
%% split data for parallel processing
[split_data, n_processes] = split_data_to_n_OSNAP(data_info_table,options.n_processes,"shuffle",true);
%% plot
max_radial_density = max(table2array(feature_data(:,contains(feature_data.Properties.VariableNames, 'radial_loc_density'))),[],'all');
max_radial_hetero_density = max(table2array(feature_data(:,contains(feature_data.Properties.VariableNames, 'radial_dbscan_cluster_density'))),[],'all');
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        try
            plot_OSNAP_ellipses_locs_and_clusters(filepath, max_radial_density, max_radial_hetero_density);
        catch ME
            disp(getReport(ME))
        end
    end
end
end


