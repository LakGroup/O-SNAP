function generate_SNAP_features(work_dir, groups, reps, options)
arguments
    work_dir string
    groups cell
    reps cell
    options.n_processes double = 12;
    options.min_log_vor_density double = 0;
    options.max_log_vor_density double = 3;
    options.heterochromatin_thresh double = 70;
    options.eps double = 20;
    options.min_num double = 3;
    options.ellipse_inc double = 0.1;
    options.periphery_thresh double = 0.15;
    options.area_threshold_arr double = 10.^(1:0.25:2);
    options.min_number_of_localizations_arr double = [5];
    options.plot logical = true;
end
plot_flag = options.plot;
warning('off','all')

%% load data
data_info_table = get_valid_voronoi_data(work_dir,groups,reps,{'x','y'});

%% split full data for parallel processing
[split_data, options.n_processes] = split_files_for_parallel(data_info_table, options.n_processes, true);

%% start parallel pool
p_pool = gcp('nocreate');
if isempty(p_pool)
    parpool("Processes",options.n_processes);
end

%% calculate voronoi density
disp('Calculating Voronoi densities...');
if plot_flag
    vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density','faces','vertices'};
else
    vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density'};
end
min_log_vor_density = options.min_log_vor_density;
max_log_vor_density = options.max_log_vor_density;
tic
parfor p=1:options.n_processes
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
                disp("Error in " + filepath + ", removing from analysis")
                split_data{p}(s,:) = [];
            end
        end
    end
end
toc
disp(['Completed: ' char(datetime)])
%% calculate density threshold
density_data=cell(1,options.n_processes);
if ~exist('density_threshold','var')
    disp('Calculating density threshold...');
    tic
    parfor p=1:options.n_processes
        data_info_table_p = split_data{p};
        filepaths = data_info_table_p{:,'filepath'};
        density_data{p} = cell(1,length(filepaths));
        for s=1:length(filepaths)
            sample = load(filepaths(s),'voronoi_areas');
            density_data{p}{s} = 1./sample.voronoi_areas;
        end
        density_data{p} = vertcat(density_data{p}{:});
    end
    density_threshold = prctile(vertcat(density_data{:}),options.heterochromatin_thresh);
    clearvars density_data
    toc
end
disp(['Completed: ' char(datetime)])
%% perform heterochromatin analysis
disp("Performing Heterochromatin Analysis...")
hetero_data_vars = {'nucleus_radius',...
            'locs_number',...
            'locs_density',...
            'boundary',... 
            'hetero_flt_with_label',...
            'lads',...
            'non_lads',...
            'lads2total',...
            'domain_boundary',...
            'hetero_center',...
            'hetero_radius',...
            'hetero_n_locs',...
            'hetero_density',...
            'lads_area',...
            'lads_n_locs',...
            'lads_density',...
            'spacing',...
            'seg_normals',...
            'seg_locs',...
            'lads_seg_len',...
            'lads_seg_thickness',...
            'components',...
            'eigenvalues',...
            'locs_norm',...
            'x_length',...
            'y_length',...
            'polygon',...
            'elastic_energy',...
            'bending_energy',...
            'border_curvature',...
            'radial_density',...
            'radial_hetero_density',...
            'interior_density',...
            'periphery_density',...
            'interior_hetero_density',...
            'periphery_hetero_density'};
eps = options.eps;
min_num = options.min_num;
ellipse_inc = options.ellipse_inc;
periphery_thresh = options.periphery_thresh;
tic
parfor p=1:options.n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file') && ~has_variables(filepath,hetero_data_vars,"verbose",0)
            try
                analyze_heterochromatin(filepath,density_threshold, ...
                    "eps",eps,...
                    "min_num",min_num, ...
                    "ellipse_inc",ellipse_inc,...
                    "periphery_thresh",periphery_thresh);
            catch ME
                disp(getReport(ME));
                disp("Removing from analysis: " + filepath)
                split_data{p}(s,:) = [];
            end
        end
    end
end
toc
disp(['Completed: ' char(datetime)])
%% perform voronoi clustering analysis
disp("Performing Voronoi Clustering Analysis...")
cluster_data_vars1 = {'area_thresholds','min_number_of_localizations','clusters'};
cluster_data_vars2 = {'cluster_n_locs','cluster_area','cluster_density','cluster_gyration_R'};
area_threshold_arr = options.area_threshold_arr;
min_number_of_localizations_arr = options.min_number_of_localizations_arr;
tic
parfor p=1:options.n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file')
            if ~has_variables(filepath,cluster_data_vars1,"verbose",0)
                try
                    generate_voronoi_clusters(filepath, ...
                        "area_threshold_arr",area_threshold_arr,...
                        "min_number_of_localizations_arr",min_number_of_localizations_arr);
                catch ME
                    disp(getReport(ME));
                    disp("Removing from analysis: " + filepath)
                    split_data{p}(s,:) = [];
                end
            end
            if ~has_variables(filepath,cluster_data_vars2,"verbose",0) && has_variables(filepath,cluster_data_vars1,"verbose",0)
                try
                    extract_cluster_features(filepath);
                catch ME
                    disp(getReport(ME));
                    disp("Removing from analysis: " + filepath)
                    split_data{p}(s,:) = [];
                end
            end
        end 
    end
end
toc
%%
disp(['  RUN END: ' char(datetime)])
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
warning('on','all')
end

