function generate_SNAP_features(work_dir, groups, reps, options)
arguments
    work_dir string
    groups cell
    reps cell
    options.n_workers double = maxNumCompThreads;
    options.min_log_vor_density double = 0;
    options.max_log_vor_density double = 3;
    options.dbscan_thresh double = 70;
    options.eps double = 20;
    options.min_num double = 3;
    options.ellipse_inc double = 0.1;
    options.periphery_thresh double = 0.15;
    options.area_threshold_arr double = 10.^(1:0.25:2);
    options.min_number_of_localizations double = 5;
    options.plot logical = true;
    options.logID = [];
end
plot_flag = options.plot;
warning('off','all')

%% load data
SNAP_nucleus_file_list = get_valid_SNAP_nucleus_files(work_dir,groups,reps,{'x','y'});

%% split full data for parallel processing
[split_file_list, options.n_workers] = split_data_to_n(SNAP_nucleus_file_list,options.n_workers);

%% calculate voronoi density
min_log_vor_density = options.min_log_vor_density;
max_log_vor_density = options.max_log_vor_density;
fprintf("      Calculating Voronoi densities...\n");
starttime_step = tic;
parfor p=1:options.n_workers
    data_info_table_p = split_file_list{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file')
            try
                generate_voronoi_data(filepath, ...
                    "min_log_vor_density", min_log_vor_density, ...
                    "max_log_vor_density", max_log_vor_density,...
                    "plot", plot_flag);
            catch ME
                split_file_list{p} = remove_sample(ME,split_file_list{p},s)
            end
        end
    end
end
fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
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
            'interior_dbscan_cluster_n_locs',...
            'interior_dbscan_cluster_density',...
            'periphery_dbscan_cluster_center',...
            'periphery_dbscan_cluster_radius',...
            'periphery_dbscan_cluster_n_locs',...
            'periphery_dbscan_cluster_density',...
            'model_lad_thickness',...
            'model_lad_segment_locs',...
            'model_lad_segment_length',...
            'model_lad_segment_thickness',...
            'model_lad_segment_spacing',...
            'model_lad_segment_normals',...
            'components',...
            'eigenvalues',...
            'locs_norm',...
            'x_length',...
            'y_length',...
            'polygon',...
            'elasstarttime_step = tic;_energy',...
            'bending_energy',...
            'border_curvature',...
            'radial_loc_density',...
            'radial_dbscan_cluster_density',...
            'interior_loc_density',...
            'periphery_loc_density',...
            'interior_cluster_density',...
            'periphery_cluster_density'};
% calculate density threshold
if ~all(arrayfun(@(x) has_variables(x,data_vars,"verbose",0),SNAP_nucleus_file_list{:,'filepath'},'uni',1))
    eps = options.eps;
    min_num = options.min_num;
    ellipse_inc = options.ellipse_inc;
    periphery_thresh = options.periphery_thresh;
    density_data=cell(1,options.n_workers);
    fprintf("      Calculating density threshold...\n");
    starttime_step = tic;
    parfor p=1:options.n_workers
        data_info_table_p = split_file_list{p};
        filepaths = data_info_table_p{:,'filepath'};
        density_data{p} = cell(1,length(filepaths));
        for s=1:length(filepaths)
            try
                sample = load(filepaths(s),'voronoi_areas');
                density_data{p}{s} = 1./sample.voronoi_areas;
            catch ME
                split_file_list{p} = remove_sample(ME,split_file_list{p},s)
            end
        end
        density_data{p} = vertcat(density_data{p}{:});
    end
    density_threshold = prctile(vertcat(density_data{:}),options.dbscan_thresh);
    clearvars density_data
    fprintf("      Calculating morphometric, dbscan clustering, and radial features...\n")
    parfor p=1:options.n_workers
        data_info_table_p = split_file_list{p};
        filepaths = data_info_table_p{:,'filepath'};
        for s=1:length(filepaths)
            filepath = filepaths(s);
            if exist(filepath,'file') && ~has_variables(filepath,data_vars,"verbose",0)
                try
                    generate_morphometrics_dbscan_radial_data(filepath,density_threshold, ...
                        "eps",eps,...
                        "min_num",min_num, ...
                        "ellipse_inc",ellipse_inc,...
                        "periphery_thresh",periphery_thresh);
                catch ME
                    split_file_list{p} = remove_sample(ME,split_file_list{p},s)
                end
            end
        end
    end
    fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
end
%% perform voronoi clustering analysis
fprintf("      Performing Voronoi Clustering Analysis...\n")
area_threshold_arr = options.area_threshold_arr;
min_number_of_localizations_arr = options.min_number_of_localizations;
starttime_step = tic;
parfor p=1:options.n_workers
    data_info_table_p = split_file_list{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file')
            try
                generate_voronoi_clusters(filepath,area_threshold_arr,min_number_of_localizations_arr);
            catch ME
                split_file_list{p} = remove_sample(ME,split_file_list{p},s)
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

function split_file_list = remove_sample(ME,split_file_list,s)
arguments
    ME
    split_file_list table
    s double
end
    fprintf("- - Removing: %s - -\n",filepath)
    fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
    fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
    split_file_list(s,:) = [];
end