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
    options.min_number_of_localizations_arr double = [5];
    options.plot logical = true;
    options.logID = [];
end
plot_flag = options.plot;
warning('off','all')

%% load data
SNAP_nucleus_file_list = get_valid_SNAP_nucleus_files(work_dir,groups,reps,{'x','y'});

%% split full data for parallel processing
[split_data, options.n_workers] = split_files_for_parallel(SNAP_nucleus_file_list, options.n_workers, true);

%% calculate voronoi density
if plot_flag
    vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density','faces','vertices'};
else
    vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density'};
end
min_log_vor_density = options.min_log_vor_density;
max_log_vor_density = options.max_log_vor_density;
disp("      Calculating Voronoi densities...");
tic
parfor p=1:options.n_workers
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file') && ~has_variables(filepath,vor_data_vars,"verbose",0)
            try
                generate_voronoi_data(filepath, ...
                    "min_log_vor_density", min_log_vor_density, ...
                    "max_log_vor_density", max_log_vor_density,...
                    "plot", plot_flag);
            catch ME
                disp(getReport(ME));
                disp("        Removing from analysis: " + filepath)
                split_data{p}(s,:) = [];
            end
        end
    end
end
toc
disp(['      Completed: ' char(datetime)])
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
            'elastic_energy',...
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
    disp("      Calculating density threshold...");
    tic
    parfor p=1:options.n_workers
        data_info_table_p = split_data{p};
        filepaths = data_info_table_p{:,'filepath'};
        density_data{p} = cell(1,length(filepaths));
        for s=1:length(filepaths)
            try
                sample = load(filepaths(s),'voronoi_areas');
                density_data{p}{s} = 1./sample.voronoi_areas;
            catch ME
                disp(getReport(ME));
                disp("        Removing from analysis: " + filepaths(s))
                split_data{p}(s,:) = [];
            end
        end
        density_data{p} = vertcat(density_data{p}{:});
    end
    density_threshold = prctile(vertcat(density_data{:}),options.dbscan_thresh);
    clearvars density_data
    disp("      Calculating morphometric, dbscan clustering, and radial features...")
    parfor p=1:options.n_workers
        data_info_table_p = split_data{p};
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
                    disp(getReport(ME));
                    disp("        Removing from analysis: " + filepath)
                    split_data{p}(s,:) = [];
                end
            end
        end
    end
    toc
    disp(['      Completed: ' char(datetime)])
end
%% perform voronoi clustering analysis
voronoi_cluster_data_vars1 = {'area_thresholds','min_number_of_localizations','voronoi_clusters'};
voronoi_cluster_data_vars2 = {'voronoi_cluster_n_locs','voronoi_cluster_radius','voronoi_cluster_density','voronoi_cluster_gyration_radius'};
if ~all(arrayfun(@(x) has_variables(x,[voronoi_cluster_data_vars1 voronoi_cluster_data_vars2],"verbose",0),SNAP_nucleus_file_list{:,'filepath'},'uni',1))
    disp("      Performing Voronoi Clustering Analysis...")
    area_threshold_arr = options.area_threshold_arr;
    min_number_of_localizations_arr = options.min_number_of_localizations_arr;
    tic
    parfor p=1:options.n_workers
        data_info_table_p = split_data{p};
        filepaths = data_info_table_p{:,'filepath'};
        for s=1:length(filepaths)
            filepath = filepaths(s);
            if exist(filepath,'file')
                if ~has_variables(filepath,voronoi_cluster_data_vars1,"verbose",0)
                    try
                        generate_voronoi_clusters(filepath, ...
                            "area_threshold_arr",area_threshold_arr,...
                            "min_number_of_localizations_arr",min_number_of_localizations_arr);
                    catch ME
                        disp(getReport(ME));
                        disp("        Removing from analysis: " + filepath)
                        split_data{p}(s,:) = [];
                    end
                end
                if ~has_variables(filepath,voronoi_cluster_data_vars2,"verbose",0) && has_variables(filepath,voronoi_cluster_data_vars1,"verbose",0)
                    try
                        extract_cluster_features(filepath);
                    catch ME
                        disp(getReport(ME));
                        disp("        Removing from analysis: " + filepath)
                        split_data{p}(s,:) = [];
                    end
                end
            end 
        end
    end
    toc
end
%%
disp(['  RUN END: ' char(datetime)])
disp("- - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
warning('on','all')
end

