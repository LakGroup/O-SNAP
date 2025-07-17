% -------------------------------------------------------------------------
% extract_OSNAP_features_samples.m
% -------------------------------------------------------------------------
% Obtains the O-SNAP features for every sample under the specified
% phenotype (groups) and replicates. The feature information is stored in
% the file that corresponds to the samples of the given groups and
% replicates.
%
% Example on how to use it:
%   extract_OSNAP_features_samples("D:\Analysis\ExperimentA",...
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
% Options:
%   n_processes: Number of cores to run parallel processes on
%   min_log_vor_density: Lower bound for filtering and visualizing voronoi
%                        density values
%   max_log_vor_density: Upper bound for filtering and visualizing voronoi
%                        density values
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
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
function extract_OSNAP_features_samples(work_dir, groups, replicates, options)
arguments
    work_dir string
    groups cell
    replicates cell
    options.n_processes double = maxNumCompThreads;
    options.min_log_vor_density double = 0;
    options.max_log_vor_density double = 3;
    options.dbscan_thresh double = 70;
    options.eps double = 20;
    options.min_num double = 3;
    options.ellipse_inc double = 0.1;
    options.periphery_thresh double = 0.15;
    options.area_threshold_arr double = 10.^(1.25:0.25:2);
    options.min_number_of_localizations double = 5;
    options.plot logical = true;
    options.overwrite logical = false;
end
min_log_vor_density = options.min_log_vor_density;
max_log_vor_density = options.max_log_vor_density;
area_threshold_arr = options.area_threshold_arr;
min_number_of_localizations_arr = options.min_number_of_localizations;
plot_flag = options.plot;
overwrite = options.overwrite;
warning('off','all')

%% load data and split for processing 
% smaller n_processes for voronoi density due to high memory demand
n_processes = ceil(options.n_processes/3);
OSNAP_sample_file_list = get_valid_OSNAP_samples(work_dir,groups,replicates,{'x','y'});
[split_file_list, n_processes] = split_data_to_n_OSNAP(OSNAP_sample_file_list,n_processes);

%% calculate voronoi density
data_vars = {...
            'voronoi_areas_all',...
            'voronoi_neighbors',...
            'voronoi_areas',...
            'reduced_log_voronoi_density',...
            'faces',...
            'vertices'};
fprintf("      Calculating Voronoi densities...\n");
starttime_step = tic;
parfor p=1:n_processes
    OSNAP_sample_file_list_p = split_file_list{p};
    filepaths = OSNAP_sample_file_list_p{:,'filepath'};
    samps_to_remove = [];
    for s=1:length(filepaths)
        filepath = filepaths(s);
        fprintf("Worker %2.0f: File %3.0f/%3.0f...\n",p,s,numel(filepaths));
        if exist(filepath,'file') && ~has_variables_OSNAP(filepath,data_vars,"verbose",0)
            try
                generate_OSNAP_voronoi_segmentations(filepath, ...
                    "min_log_vor_density", min_log_vor_density, ...
                    "max_log_vor_density", max_log_vor_density,...
                    "plot", plot_flag,"overwrite",overwrite);
            catch ME
                fprintf("- - Removing: %s - -\n",filepath)
                fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                samps_to_remove = [samps_to_remove s];
                fprintf("- - - - - - - - - - -\n")
            end
        end
    end
end
fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);

%% load data and split for processing 
OSNAP_sample_file_list = get_valid_OSNAP_samples(work_dir,groups,replicates,{'x','y'});
[split_file_list, n_processes] = split_data_to_n_OSNAP(OSNAP_sample_file_list,options.n_processes);
%% perform morphometric, dbscan clustering, and radial analysis
data_vars = {'nucleus_radius',...
            'locs_number',...
            'locs_density',...
            'boundary',... 
            'locs_dbscan_cluster_labels',...
            'locs_dbscan_cluster_periphery_labels',...
            'locs_dbscan_cluster_interior_labels',...
            'domain_boundary',...
            'interior_dbscan_cluster_center',...
            'interior_dbscan_cluster_radius',...
            'interior_dbscan_cluster_gyration_radius',...
            'interior_dbscan_cluster_n_locs',...
            'interior_dbscan_cluster_density',...
            'interior_dbscan_cluster_spacing',...
            'periphery_dbscan_cluster_center',...
            'periphery_dbscan_cluster_radius',...
            'periphery_dbscan_cluster_gyration_radius',...
            'periphery_dbscan_cluster_n_locs',...
            'periphery_dbscan_cluster_density',...
            'model_lad_thickness',...
            'model_lad_segment_locs',...
            'model_lad_segment_length',...
            'model_lad_segment_thickness',...
            'model_lad_segment_normals',...
            'components',...
            'eigenvalues',...
            'locs_norm',...
            'x_length',...
            'y_length',...
            'polygon',...
            'elastic_energy',...
            'bending_energy',...
            'border_curvature',...
            'radial_loc_density',...
            'radial_dbscan_cluster_density',...
            'interior_loc_density',...
            'periphery_loc_density',...
            'interior_cluster_density',...
            'periphery_cluster_density'};
%% calculate density threshold
if ~all(arrayfun(@(x) has_variables_OSNAP(x,data_vars,"verbose",0),OSNAP_sample_file_list{:,'filepath'},'uni',1))
    eps = options.eps;
    min_num = options.min_num;
    ellipse_inc = options.ellipse_inc;
    periphery_thresh = options.periphery_thresh;
    density_data=cell(1,n_processes);
    fprintf("      Calculating density threshold...\n");
    starttime_step = tic;
    parfor p=1:n_processes
        OSNAP_sample_file_list_p = split_file_list{p};
        filepaths = OSNAP_sample_file_list_p{:,'filepath'};
        density_data{p} = cell(1,length(filepaths));
        samps_to_remove = [];
        for s=1:numel(filepaths)
            try
                sample = load(filepaths(s),'voronoi_areas');
                density_data{p}{s} = 1./sample.voronoi_areas;
            catch ME
                fprintf("- - Removing: %s - -\n",filepaths(s))
                fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                samps_to_remove = [samps_to_remove s];
            end
        end
        density_data{p} = vertcat(density_data{p}{:});
        split_file_list{p}(samps_to_remove,:) = [];
    end
    density_threshold = prctile(vertcat(density_data{:}),options.dbscan_thresh);
    clearvars density_data
    fprintf("      Calculating morphometric, dbscan clustering, and radial features...\n")
    parfor p=1:n_processes
        OSNAP_sample_file_list_p = split_file_list{p};
        filepaths = OSNAP_sample_file_list_p{:,'filepath'};
        samps_to_remove = [];
        for s=1:length(filepaths)
            filepath = filepaths(s);
            fprintf("Worker %2.0f: File %3.0f/%3.0f...\n",p,s,numel(filepaths));
            if exist(filepath,'file') && ~has_variables_OSNAP(filepath,data_vars,"verbose",0)
                try
                    extract_OSNAP_morphometrics_dbscan_cluster_radial_features(filepath,density_threshold, ...
                        "eps",eps,...
                        "min_num",min_num, ...
                        "ellipse_inc",ellipse_inc,...
                        "periphery_thresh",periphery_thresh,...
                        "overwrite",overwrite);
                catch ME
                    fprintf("- - Removing: %s - -\n",filepath)
                    fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                    samps_to_remove = [samps_to_remove s];
                end
            end
        end
        split_file_list{p}(samps_to_remove,:) = [];
    end
    fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
end
%% perform voronoi clustering analysis
data_vars = {...
    'area_thresholds',...
    'min_number_of_localizations',...
    'voronoi_clusters',...
    'voronoi_cluster_n_locs',...
    'voronoi_cluster_radius',...
    'voronoi_cluster_density',...
    'voronoi_cluster_gyration_radius'};
fprintf("      Performing Voronoi Clustering Analysis...\n")
starttime_step = tic;
parfor p=1:n_processes
    OSNAP_sample_file_list_p = split_file_list{p};
    filepaths = OSNAP_sample_file_list_p{:,'filepath'};
    for s=1:length(filepaths) 
        filepath = filepaths(s);
        fprintf("Worker %2.0f: File %3.0f/%3.0f...\n",p,s,numel(filepaths));
        if exist(filepath,'file') && ~has_variables_OSNAP(filepath,data_vars,"verbose",0)
            try
                extract_OSNAP_voronoi_cluster_features(filepath,...
                    area_threshold_arr,...
                    min_number_of_localizations_arr,...
                    "overwrite",overwrite);
            catch ME
                fprintf("- - Removing: %s - -\n",filepath)
                fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
            end
        end 
    end
end
fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
%%
fprintf("  RUN END: %s\n",datetime);
fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
warning('on','all')
end


