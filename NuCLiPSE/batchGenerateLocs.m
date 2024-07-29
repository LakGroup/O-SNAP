function batchGenerateLocs(data, workDir)
%% SET UP FOR ANALYSIS
% create save directory
if workDir == 0
    error('Error: Directory selection canceled');
else
    saveDir = fullfile(workDir,'VoronoiAreaData');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir)
    end
    mapDir = fullfile(workDir,'ReducedVoronoiDensityMaps');
    if ~exist(mapDir, 'dir')
        mkdir(mapDir)
    end
end

% parameters for voronoi density plots
inputValues = inputdlg({
    'Pixel Size (nm per pixel):',...
    'Number of processes:'},...
    'Input', ...
    [1 50],...
    {'117',...
    '12'});
pixel_size = str2double(inputValues{1});
nProcesses = str2double(inputValues{2});

% split data
[splitData, nProcesses] = splitDataForParallel(data, nProcesses);

% save scaled localizations
disp('Saving localizations...');
tic
parfor p=1:nProcesses
    data_p = splitData{p};
    for s=1:length(data_p)
        sample = data_p{s};
        %% GENERATE VORONOI DENSITY MAP
        sample = generateLocData(sample, pixel_size);
        saveVoronoiAnalysisData(fullfile(saveDir,[sample.name '.mat']), sample);
        disp("   " + sample.name + ": " + string(datetime))
    end
end
toc

disp('Done!')
end 