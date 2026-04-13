% -------------------------------------------------------------------------
% generate_OSNAP_features_samplelist.m
% -------------------------------------------------------------------------
% Obtains the O-SNAP features for every sample under the specified
% phenotype (groups) and replicates. The feature information is stored in
% the file that corresponds to the samples of the given groups and
% replicates.
%
% Example on how to use it:
%   generate_OSNAP_features_samples("D:\Analysis\ExperimentA",...
%                                 {'Control','KO'},...
%                                 {'20250101','20250108','20250201'})
% -------------------------------------------------------------------------
% Input:
%   OSNAP_sample_file_list: Table with details on the samples of interest
%                           to generate feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
% Options:
%   n_processes: Number of cores to run parallel processes on
%   min_log_vor_density: Lower bound for filtering and visualizing reduced 
%                        Voronoi density values
%   max_log_vor_density: Upper bound for filtering and visualizing reduced
%                        Voronoi density values
%   min_log_vor_area: Lower bound for filtering and visualizing Voronoi
%                     area values
%   max_log_vor_area: Upper bound for filtering and visualizing Voronoi
%                     area values
%   dbscan_thresh: Lower percentile cutoff for localizations to include in
%                  the DBSCAN clustering analysis for heterochromatin-
%                  associated domains
%   eps: The search radius of the DBSCAN clustering algorithm
%   min_num: The minimum number of localizations to include in a DBSCAN
%            domain
%   ellipse_inc: For radial analysis of localization/cluster density.
%                Percentage of the major/minor axis length to subdivide the
%                ellipses that define the bounds of the rings. The number 
%                of rings is N_r = floor(1/inc)
%   periphery_thresh: Cutoff for what percentage of the polygon to consider
%                     as the periphery. Ex. 0.15 results in a periphery
%                     region with an area equivalent to 15% of the original
%                     polygon shape
%   area_threshold_arr: Area thresholds for multi-scale Voronoi clustering.
%                       Act as upper bound for Voronoi cell areas. Locs
%                       with areas under this value are filtered.
%   min_number_of_localizations: The minimum number of locs in a Voronoi 
%                                cluster. Clusters with fewer locs are
%                                filtered.
%   plot: Flag to plot results
%   overwrite: Flag on whether to ask user before overwriting data
%   filter: Flag to retain samples where names do not satisfy a condition
%           (group name in information table is left as an empty string)
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
function generate_OSNAP_features_samplelist(OSNAP_sample_file_list, options)
arguments
    OSNAP_sample_file_list table
    options.n_processes double = maxNumCompThreads;
    options.min_log_reduced_vor_density double = 0.5;
    options.max_log_reduced_vor_density double = 3;
    options.min_log_vor_area double = 1;
    options.max_log_vor_area double = 3;
    options.dbscan_thresh double = 70;
    options.eps double = 20;
    options.min_num double = 3;
    options.ellipse_inc double = 0.1;
    options.periphery_thresh double = 0.15;
    options.area_threshold_arr double = 10.^(1.25:0.25:2);
    options.min_number_of_localizations double = 5;
    options.plot logical = true;
    options.overwrite logical = false;
    options.filter logical = false;
end
warning('off','all')

%% Calculate voronoi density
generate_OSNAP_features_samplelist_voronoi_segmentations(OSNAP_sample_file_list,...
    "n_processes",options.n_processes,...
    "min_log_reduced_vor_density",options.min_log_reduced_vor_density,...
    "max_log_reduced_vor_density",options.max_log_reduced_vor_density,...
    "min_log_vor_area",options.min_log_vor_area,...
    "max_log_vor_area",options.max_log_vor_area,...
    "plot",options.plot_flag,...
    "overwrite",options.overwrite...
    );
%% Perform morphometric, dbscan clustering, and radial analysis
generate_OSNAP_features_samplelist_morphometrics_DBSCAN_radial(OSNAP_sample_file_list,...
    "n_processes",options.n_processes,...
    "dbscan_thresh",options.dbscan_thresh,... 
    "eps",options.eps,... 
    "min_num",options.min_num,... 
    "ellipse_inc",options.ellipse_inc,...
    "periphery_thresh",options.periphery_thresh,... 
    "plot",options.plot,...
    "overwrite",options.overwrite... 
    );
%% Perform voronoi clustering analysis
generate_OSNAP_features_samplelist_voronoi_segmentations(OSNAP_sample_file_list,...
    "n_processes",options.n_processes,... 
    "area_threshold_arr",options.area_threshold_arr,... 
    "min_number_of_localizations",options.min_number_of_localizations,... 
    "overwrite",options.overwrite...
    );
%% Conclude output
fprintf("  RUN END: %s\n",datetime);
fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
warning('on','all')
end


