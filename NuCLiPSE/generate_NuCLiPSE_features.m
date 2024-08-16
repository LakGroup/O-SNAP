function generate_NuCLiPSE_features(work_dir, groups, reps)
%% analysis parameters
n_processes = 12;
min_log_vor_density = 0;
max_log_vor_density = 3;
heterochromatin_thresh = 70;
eps = 20;
min_num = 3;
ellipse_inc = 0.1;
periphery_thresh = 0.15;
area_threshold_arr = 10.^(1:0.25:2);
min_number_of_localizations_arr = [5];

%% load data
data_info_table = get_valid_voronoi_data(work_dir,groups,reps,{'x','y'});

%% split full data for parallel processing
[split_data, n_processes] = split_files_for_parallel(data_info_table, n_processes, true);

%% start parallel pool
p_pool = gcp('nocreate');
if isempty(p_pool)
    parpool("Processes",n_processes);
end

%% calculate voronoi density
disp('Calculating Voronoi densities...');
vor_data_vars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density'};
tic
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file') && ~has_variables(filepath,vor_data_vars,"verbose",0)
            try
                generate_voronoi_data(filepath, min_log_vor_density, max_log_vor_density);
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
density_data=cell(1,n_processes);
if ~exist('density_threshold','var')
    disp('Calculating density threshold...');
    tic
    parfor p=1:n_processes
        data_info_table_p = split_data{p};
        filepaths = data_info_table_p{:,'filepath'};
        density_data{p} = cell(1,length(filepaths));
        for s=1:length(filepaths)
            sample = load(filepaths(s),'voronoi_areas');
            density_data{p}{s} = 1./sample.voronoi_areas;
        end
        density_data{p} = vertcat(density_data{p}{:});
    end
    density_threshold = prctile(vertcat(density_data{:}),heterochromatin_thresh);
    clearvars density_data
    toc
end
disp(['Completed: ' char(datetime)])
%% perform Heterochromatin analysis
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
tic
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file') && ~has_variables(filepath,hetero_data_vars,"verbose",0)
            try
                analyze_heterochromatin(filepath,density_threshold,eps,min_num,ellipse_inc,periphery_thresh);
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
tic
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file')
            if ~has_variables(filepath,cluster_data_vars1,"verbose",0)
                try
                    generate_voronoi_clusters(filepath,area_threshold_arr,min_number_of_localizations_arr);
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
end
