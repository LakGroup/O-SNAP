function plot_radial_densities(data_info_table, T)

%% parameters
n_processes = 12;

%% split data for parallel processing
[split_data, n_processes] = split_files_for_parallel(data_info_table, n_processes, true);
% start parallel pool
p_pool = gcp('nocreate');
disp('Generating table data...')
tic
if isempty(p_pool)
    parpool("Processes",n_processes);
end

%% plot
max_radial_density = max(table2array(T(:,contains(T.Properties.VariableNames, 'radial_density'))),[],'all');
max_radial_hetero_density = max(table2array(T(:,contains(T.Properties.VariableNames, 'radial_hetero_density'))),[],'all');
parfor p=1:n_processes
    data_info_table_p = split_data{p};
    filepaths = data_info_table_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        batchplot_radial_density(filepath, max_radial_density, max_radial_hetero_density);
    end
end

end
