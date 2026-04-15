% -------------------------------------------------------------------------
% extract_OSNAP_features_samples_morphometrics_dbscan_cluster_radial.m
% -------------------------------------------------------------------------
% Coordinates calculation of features related to morphometrics, DBSCAN 
% heterochromatin-associated localization clustering, and radial features 
% for multiple samples (nuclei). Saves information to the MAT analysis 
% file specific to that sample.
%
% Example on how to use it:
%   generate_OSNAP_features_samplelist_morphometrics_DBSCAN_radial(OSNAP_sample_file_list)
% -------------------------------------------------------------------------
% Input:
%   OSNAP_sample_file_list: Table with details on the samples of interest
%                           to extract feature information from, including 
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
function generate_OSNAP_features_samplelist_morphometrics_DBSCAN_radial(OSNAP_sample_file_list, options)
arguments
    OSNAP_sample_file_list table
    options.n_processes double = maxNumCompThreads;
    options.dbscan_thresh double = 70;
    options.eps double = 20;
    options.min_num double = 3;
    options.ellipse_inc double = 0.1;
    options.periphery_thresh double = 0.15;
    options.plot logical = true;
    options.overwrite logical = false;
end

n_processes = options.n_processes;
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
% Identify files to process
idx = ~logical(arrayfun(@(x) has_variables_OSNAP(x,data_vars,"verbose",0), OSNAP_sample_file_list{:,'filepath'},'uni',1));
[split_file_list, n_processes] = split_data_to_n_OSNAP(OSNAP_sample_file_list(idx,:),n_processes);
starttime_step = tic;
%% Calculate density threshold
if any(cellfun('length',split_file_list))
    eps = options.eps;
    min_num = options.min_num;
    ellipse_inc = options.ellipse_inc;
    periphery_thresh = options.periphery_thresh;
    plot_flag = options.plot;
    overwrite = options.overwrite;
    density_data=cell(1,n_processes);
    fprintf("      Calculating density threshold...\n");
    parfor p=1:n_processes
        sample_file_list_p = split_file_list{p};
        filepaths = sample_file_list_p{:,'filepath'};
        if numel(filepaths)==0
            continue;
        end
        density_data{p} = cell(1,length(filepaths));
        samps_to_remove = [];
        for s=1:numel(filepaths)
            try
                sample = load(filepaths(s),'voronoi_areas');
                density_data{p}{s} = 1./sample.voronoi_areas;
            catch ME
                fprintf("- - Skipping: %s - -\n",filepaths(s))
                fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                samps_to_remove = [samps_to_remove s];
            end
        end
        % Remove samples that give errors
        density_data{p} = vertcat(density_data{p}{:});
        split_file_list{p}(samps_to_remove,:) = [];
    end
    density_threshold = prctile(vertcat(density_data{:}),options.dbscan_thresh);
    fprintf("          Density threshold: %.4f...\n",density_threshold);
    clearvars density_data
    fprintf("      Calculating morphometric, dbscan clustering, and radial features...\n")
    parfor p=1:n_processes
        sample_file_list_p = split_file_list{p};
        filepaths = sample_file_list_p{:,'filepath'};
        if numel(filepaths)==0
            continue;
        end
        % samps_to_remove = [];
        for s=1:length(filepaths)
            filepath = filepaths(s);
            fprintf("Worker %2.0f: File %3.0f/%3.0f...\t%s\n",p,s,numel(filepaths),filepath);
            if exist(filepath,'file') && ~has_variables_OSNAP(filepath,data_vars,"verbose",0)
                try
                    generate_OSNAP_features_sample_morphometrics_DBSCAN_radial(filepath,density_threshold, ...
                        "eps",eps,...
                        "min_num",min_num, ...
                        "ellipse_inc",ellipse_inc,...
                        "periphery_thresh",periphery_thresh,...
                        "plot",plot_flag,...
                        "overwrite",overwrite);
                catch ME
                    fprintf("- - Skipping: %s - -\n",filepath)
                    fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                    % samps_to_remove = [samps_to_remove s];
                end
            end
        end
        % split_file_list{p}(samps_to_remove,:) = [];
        fprintf("Worker %2.0f: COMPLETE\n",p);
    end
end
fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);