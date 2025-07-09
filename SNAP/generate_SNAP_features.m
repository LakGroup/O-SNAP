% -------------------------------------------------------------------------
% generate_SNAP_features.m
% -------------------------------------------------------------------------
% <FILL>
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
%   work_dir:  The parent directory that contains the directory where the
%              analysis results are saved
%   analysis_name: Analysis identifier string
%   groups: Cell array containing char array of the identifiers of the 
%           phenotypes/cell states *ENSURE THAT EACH SAMPLE FILENAME 
%           CONTAINS EXACTLY ONE IDENTIFIER FROM groups*
%   replicates: Cell array containing char array of the identifiers of the 
%               replicates. *ENSURE THAT EACH REPLCIATE IS STORED IN A
%               SEPARATE FOLDER PREFIXED BY "Rep"*
% Options:
%   n_processes:  
%   min_log_vor_density:  
%   max_log_vor_density:  
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
function generate_SNAP_features(work_dir, groups, replicates, options)
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
    options.logID = [];
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
SNAP_nucleus_file_list = get_valid_SNAP_nucleus_files(work_dir,groups,replicates,{'x','y'});
[split_file_list, n_processes] = split_data_to_n(SNAP_nucleus_file_list,n_processes);

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
    data_info_table_p = split_file_list{p};
    filepaths = data_info_table_p{:,'filepath'};
    samps_to_remove = [];
    for s=1:length(filepaths)
        filepath = filepaths(s);
        fprintf("Worker %2.0f: File %3.0f/%3.0f...\n",p,s,numel(filepaths));
        if exist(filepath,'file') && ~has_variables(filepath,data_vars,"verbose",0)
            try
                generate_voronoi_data(filepath, ...
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
SNAP_nucleus_file_list = get_valid_SNAP_nucleus_files(work_dir,groups,replicates,{'x','y'});
[split_file_list, n_processes] = split_data_to_n(SNAP_nucleus_file_list,options.n_processes);
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
if ~all(arrayfun(@(x) has_variables(x,data_vars,"verbose",0),SNAP_nucleus_file_list{:,'filepath'},'uni',1))
    eps = options.eps;
    min_num = options.min_num;
    ellipse_inc = options.ellipse_inc;
    periphery_thresh = options.periphery_thresh;
    density_data=cell(1,n_processes);
    fprintf("      Calculating density threshold...\n");
    starttime_step = tic;
    parfor p=1:n_processes
        data_info_table_p = split_file_list{p};
        filepaths = data_info_table_p{:,'filepath'};
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
        data_info_table_p = split_file_list{p};
        filepaths = data_info_table_p{:,'filepath'};
        samps_to_remove = [];
        for s=1:length(filepaths)
            filepath = filepaths(s);
            fprintf("Worker %2.0f: File %3.0f/%3.0f...\n",p,s,numel(filepaths));
            if exist(filepath,'file') && ~has_variables(filepath,data_vars,"verbose",0)
                try
                    generate_morphometrics_dbscan_radial_data(filepath,density_threshold, ...
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
    data_info_table_p = split_file_list{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths) 
        filepath = filepaths(s);
        fprintf("Worker %2.0f: File %3.0f/%3.0f...\n",p,s,numel(filepaths));
        if exist(filepath,'file') && ~has_variables(filepath,data_vars,"verbose",0)
            try
                generate_voronoi_clusters(filepath,...
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

