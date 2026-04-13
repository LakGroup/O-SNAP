% -------------------------------------------------------------------------
% extract_OSNAP_features_samples_voronoi_density.m
% -------------------------------------------------------------------------
% Coordinates Voronoi segmentation and plotting for multiple samples 
% (nuclei). Saves information to the MAT analysis files.
%
% Example on how to use it:
%   extract_OSNAP_features_samples_voronoi_density(OSNAP_sample_file_list)
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
%   plot: Flag to plot results
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
function generate_OSNAP_features_samplelist_voronoi_segmentations(OSNAP_sample_file_list, options)
arguments
    OSNAP_sample_file_list table
    options.n_processes double = maxNumCompThreads;
    options.min_log_reduced_vor_density double = 0.5;
    options.max_log_reduced_vor_density double = 3;
    options.min_log_vor_area double = 1;
    options.max_log_vor_area double = 3;
    options.plot logical = true;
    options.overwrite logical = false;
    options.filter logical = false;
end
min_log_reduced_vor_density = options.min_log_reduced_vor_density;
max_log_reduced_vor_density = options.max_log_reduced_vor_density;
min_log_vor_area = options.min_log_vor_area;
max_log_vor_area = options.max_log_vor_area;
plot_flag = options.plot;
overwrite = options.overwrite;
warning('off','all')

% Smaller n_processes for voronoi density due to high memory demand
n_processes = ceil(options.n_processes/3);
fprintf("      Using %.0f workers (fewer due to high memory demand)...\n",n_processes);
% Variables to calculate
data_vars = {...
            'x','y',...
            'voronoi_areas_all',...
            'voronoi_neighbors',...
            'voronoi_areas',...
            'reduced_log_voronoi_density',...
            'faces',...
            'vertices'};
% Identify files to process
if ~overwrite
    idx = ~logical(arrayfun(@(x) has_variables_OSNAP(x,data_vars,"verbose",0), OSNAP_sample_file_list{:,'filepath'},'uni',1));
else
    idx = true(1,size(OSNAP_sample_file_list,1));
end
[split_file_list, n_processes] = split_data_to_n_OSNAP(OSNAP_sample_file_list(idx,:),n_processes);
fprintf("      Calculating Voronoi densities...\n");
starttime_step = tic;
if any(cellfun('length',split_file_list))
    parfor p=1:n_processes
        sample_file_list_p = split_file_list{p};
        filepaths = sample_file_list_p{:,'filepath'};
        % samps_to_remove = [];
        for s=1:length(filepaths)
            filepath = filepaths(s);
            fprintf("Worker %2.0f: File %3.0f/%3.0f...\t%s\n",p,s,numel(filepaths),filepath);
            try
                if exist(filepath,'file') && (~has_variables_OSNAP(filepath,data_vars,"verbose",0) || overwrite)
                    generate_OSNAP_features_sample_voronoi_segmentations(filepath, ...
                        "min_log_vor_density", min_log_reduced_vor_density, ...
                        "max_log_vor_density", max_log_reduced_vor_density,...
                        "min_log_vor_area", min_log_vor_area, ...
                        "max_log_vor_area", max_log_vor_area,...
                        "plot", plot_flag,"overwrite",overwrite);
                end
            catch ME
                fprintf("- - Skipping: %s - -\n",filepath)
                fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
                % samps_to_remove = [samps_to_remove s];
                fprintf("- - - - - - - - - - -\n")
            end
        end
        % split_file_list{p}(samps_to_remove,:) = [];
        fprintf("Worker %2.0f: COMPLETE\n",p);
    end
end
fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);

warning('on','all')
