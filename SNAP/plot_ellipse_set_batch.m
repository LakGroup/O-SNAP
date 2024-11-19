function plot_ellipse_set_batch(work_dir,groups,reps,T,options)
arguments
    work_dir string
    groups cell
    reps cell
    T table
    options.n_processes double= 12;
end
%% get data
data_info_table = get_valid_voronoi_data(work_dir,groups,reps,{'radial_density','radial_hetero_density'});
%% split data for parallel processing
[split_data, n_processes] = split_files_for_parallel(data_info_table, options.n_processes, true);
% start parallel pool
p_pool = gcp('nocreate');
if isempty(p_pool)
    parpool("Processes",options.n_processes);
end
%% plot
max_radial_density = max(table2array(T(:,contains(T.Properties.VariableNames, 'radial_density'))),[],'all');
max_radial_hetero_density = max(table2array(T(:,contains(T.Properties.VariableNames, 'radial_heterochromatin_density'))),[],'all');
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        try
            plot_ellipse_set_locs_and_hetero(filepath, max_radial_density, max_radial_hetero_density);
        catch ME
            disp(getReport(ME))
        end
    end
end
end

