function generate_localizations_batch(data, work_dir)
%% SET UP FOR ANALYSIS
% create save directory
if work_dir == 0
    error('Error: Directory selection canceled');
else
    save_dir = fullfile(work_dir,'voronoi_data');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end
end

% parameters for voronoi density plots
input_values = inputdlg({
    'Pixel Size (nm per pixel):',...
    'Number of processes:'},...
    'Input', ...
    [1 50],...
    {'117',...
    '12'});
pixel_size = str2double(input_values{1});
n_processes = str2double(input_values{2});

% split data
[split_data, n_processes] = split_data_for_parallel(data, n_processes);

% save scaled localizations
disp('Saving localizations...');
tic
parfor p=1:n_processes
    data_p = split_data{p};
    for s=1:length(data_p)
        sample = data_p{s};
        %% GENERATE VORONOI DENSITY MAP
        sample = generate_localizations(sample, pixel_size);
        save_voronoi_data(fullfile(save_dir,[sample.name '.mat']), sample);
        disp("   " + sample.name + ": " + string(datetime))
    end
end
toc

disp('Done!')
end 

