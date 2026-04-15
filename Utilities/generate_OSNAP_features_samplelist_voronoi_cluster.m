% -------------------------------------------------------------------------
% generate_OSNAP_features_samplelist_voronoi_cluster.m
% -------------------------------------------------------------------------
% Coordinates calculation of features related to Voronoi clusters
% for multiple samples (nuclei). Saves information to the MAT analysis 
% files.
%
% Example on how to use it:
%   extract_OSNAP_features_samples(OSNAP_sample_file_list)
% -------------------------------------------------------------------------
% Input:
%   OSNAP_sample_file_list: Table with details on the samples of interest
%                           to extract feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
% Options:
%   n_processes: Number of cores to run parallel processes on
%   area_threshold_arr: Area thresholds for multi-scale Voronoi clustering.
%                       Act as upper bound for Voronoi cell areas. Locs
%                       with areas under this value are filtered.
%   min_number_of_localizations: The minimum number of locs in a Voronoi 
%                                cluster. Clusters with fewer locs are
%                                filtered.
%   overwrite: Flag on whether to ask user before overwriting data
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
function generate_OSNAP_features_samplelist_voronoi_cluster(OSNAP_sample_file_list, options)
arguments
    OSNAP_sample_file_list table
    options.n_processes double = maxNumCompThreads;
    options.area_threshold_arr double = 10.^(1.25:0.25:2);
    options.min_number_of_localizations double = 5;
    options.overwrite logical = false;
end
area_threshold_arr = options.area_threshold_arr;
min_number_of_localizations = options.min_number_of_localizations;
overwrite = options.overwrite;
data_vars = {...
    'area_thresholds',...
    'min_number_of_localizations',...
    'voronoi_clusters',...
    'voronoi_cluster_n_locs',...
    'voronoi_cluster_radius',...
    'voronoi_cluster_density',...
    'voronoi_cluster_gyration_radius'};
% Identify files to process
idx = ~logical(arrayfun(@(x) has_variables_OSNAP(x,data_vars,"verbose",0), OSNAP_sample_file_list{:,'filepath'},'uni',1));
[split_file_list, n_processes] = split_data_to_n_OSNAP(OSNAP_sample_file_list(idx,:),options.n_processes);
fprintf("      Performing Voronoi Clustering Analysis...\n")
starttime_step = tic;
if any(cellfun('length',split_file_list))
    parfor p=1:n_processes
        sample_file_list_p = split_file_list{p};
        filepaths = sample_file_list_p{:,'filepath'};
        if numel(filepaths)==0
            continue;
        end
        for s=1:length(filepaths) 
            filepath = filepaths(s);
            fprintf("Worker %2.0f: File %3.0f/%3.0f...\t%s\n",p,s,numel(filepaths),filepath);
            if exist(filepath,'file') && ~has_variables_OSNAP(filepath,data_vars,"verbose",0)
                try
                    generate_OSNAP_features_sample_voronoi_cluster(filepath,...
                        area_threshold_arr,...
                        min_number_of_localizations,...
                        "overwrite",overwrite);
                catch ME
                    fprintf("- - Skipping: %s - -\n",filepath)
                    fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                end
            end 
        end
        fprintf("Worker %2.0f: COMPLETE\n",p);
    end
end
fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
end