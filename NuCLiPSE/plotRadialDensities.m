function plotRadialDensities(dataInfoTable, T)

%% parameters
nProcesses = 12;

%% split data for parallel processing
[splitData, nProcesses] = splitFilesForParallel(dataInfoTable, nProcesses, true);
% start parallel pool
pPool = gcp('nocreate');
disp('Generating table data...')
tic
if isempty(pPool)
    parpool("Processes",nProcesses);
end

%% plot
maxRadialDensity = max(table2array(T(:,contains(T.Properties.VariableNames, 'radialDensity'))),[],'all');
maxRadialHeteroDensity = max(table2array(T(:,contains(T.Properties.VariableNames, 'radialHeteroDensity'))),[],'all');
parfor p=1:nProcesses
    dataInfoTable_p = splitData{p};
    filepaths = dataInfoTable_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        batchPlotRadialDensity(filepath, maxRadialDensity, maxRadialHeteroDensity);
    end
end

end