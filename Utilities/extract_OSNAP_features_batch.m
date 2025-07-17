% -------------------------------------------------------------------------
% extract_OSNAP_features_batch.m
% -------------------------------------------------------------------------
% Coordinates the parallel processing extracting features from sample files
% into a single feature table. Each process creates its own table, and
% every table across the multiple processes are then coallated into a final
% table.
%
% Example on how to use it:
%   feature_data = extract_OSNAP_features_batch("D:\Analysis\ExperimentA",...
%                                 {'Control','KO'},...
%                                 {'20250101','20250108','20250201'})
% -------------------------------------------------------------------------
% Input:
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
% Output:
%   T: The feature data output table, where each row represents a sample
%      (nucleus). The first three columns represent (1) Group/Phenotype, 
%      (2) replicate, and (3) Sample Identifier. Each subsequent column is  
%      an O-SNAP feature. 
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
function T = extract_OSNAP_features_batch(work_dir,groups,replicates,options)
arguments
    work_dir string
    groups cell
    replicates cell
    options.n_processes double = maxNumCompThreads
end
    vars_string = {'group','biological_replicate','name'};

    %% organize into tables
    % get data
    data_info_table = get_valid_OSNAP_samples(work_dir,groups,replicates,{'voronoi_cluster_radius'});
    % split data for parallel processing
    [split_file_list, options.n_processes] = split_data_to_n_OSNAP(data_info_table, options.n_processes);
    T_cell = cell(1,options.n_processes);
    % run analysis
    parfor p=1:options.n_processes
        data_info_table_p = split_file_list{p};
        T_cell{p} = extract_OSNAP_features(data_info_table_p);
    end
    T = sortrows(vertcat(T_cell{:}),vars_string);
end



