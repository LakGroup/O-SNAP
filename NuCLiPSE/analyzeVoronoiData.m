%% setup
% workDir = uigetdir('C:\Users\HannahKim\Documents\Sessions','Select folder with data');

% Rep 0
workDir = "C:\Users\HannahKim\Documents\Sessions\ANALYSIS_HK_BovChon_Rep0";
colIVDataDirList = {
'F:\20230712_HK-BovChon-P0-CMNeg-ColVIH3-Rep0\',...
'F:\20230712_HK-BovChon-P0-CMPos-ColVIH3-Rep0\',...
'F:\20230714_HK-BovChon-P2-CMNeg-ColVIH3-Rep0\',...
'F:\20230715_HK-BovChon-P2-CMPos-ColVIH3-Rep0\',...
'F:\20230715_HK-BovChon-P6-CMNeg-ColVIH3-Rep0\',...
'F:\20230716_HK-BovChon-P6-CMPos-ColVIH3-Rep0\',...
};

% % Rep 1
% workDir = "C:\Users\HannahKim\Documents\Sessions\ANALYSIS_HK_BovChon_Rep1";
% colIVDataDirList = {
% 'G:\20230921_HK-BovChon-P0-CMNeg-ColVIH3-Rep1\',...
% 'G:\20230921_HK-BovChon-P0-CMPos-ColVIH3-Rep1\',...
% 'G:\20230924_HK-BovChon-P2-CMNeg-ColVIH3-Rep1\',...
% 'G:\20230924_HK-BovChon-P2-CMPos-ColVIH3-Rep1\',...
% 'G:\20231007_HK-BovChon-P6-CMNeg-ColVIH3-Rep1\',...
% 'G:\20231007_HK-BovChon-P6-CMPos-ColVIH3-Rep1\',...
% };

% % Rep 2
% workDir = "C:\Users\HannahKim\Documents\Sessions\ANALYSIS_HK_BovChon_Rep2";
% colIVDataDirList = {
% 'H:\20231001_HK-BovChon-P0-CMNeg-ColVIH3-Rep2',...
% 'H:\20231001_HK-BovChon-P0-CMPos-ColVIH3-Rep2',...
% 'H:\20231021_HK-BovChon-P2-CMNeg-ColVIH3-Rep2',...
% 'H:\20231021_HK-BovChon-P2-CMPos-ColVIH3-Rep2',...
% 'H:\20231017_HK-BovChon-P6-CMNeg-ColVIH3-Rep2',...
% 'H:\20231017_HK-BovChon-P6-CMPos-ColVIH3-Rep2',...
% 'H:\20231025_HK-BovChon-P0-CMNeg-ColVIH3-Rep2'
% };

% % BACKGROUND Rep 1
% workDir = "C:\Users\HannahKim\Documents\Sessions\ANALYSIS_HK_BovChon_Rep1";
% colIVDataDirList = {
% 'H:\20231012_HK-BovChon-P0-CMNeg-ColVINegative-Rep1',...
% 'H:\20231012_HK-BovChon-P0-CMPos-ColVINegative-Rep1',...
% 'H:\20231012_HK-BovChon-P2-CMNeg-ColVINegative-Rep1',...
% 'H:\20231012_HK-BovChon-P2-CMPos-ColVINegative-Rep1',...
% 'H:\20231022_HK-BovChon-P6-CMNeg-ColVINegative-Rep1',...
% 'H:\20231022_HK-BovChon-P6-CMPos-ColVINegative-Rep1'
% };

% % BACKGROUND Rep 2
% workDir = "C:\Users\HannahKim\Documents\Sessions\ANALYSIS_HK_BovChon_Rep2";
% colIVDataDirList = {
% 'H:\20231022_HK-BovChon-P0-CMNeg-ColVINegative-Rep2',...
% 'H:\20231022_HK-BovChon-P0-CMPos-ColVINegative-Rep2',...
% 'H:\20231022_HK-BovChon-P2-CMNeg-ColVINegative-Rep2',...
% 'H:\20231022_HK-BovChon-P2-CMPos-ColVINegative-Rep2',...
% 'H:\20231025_HK-BovChon-P6-CMNeg-ColVNegative-Rep2',...
% 'H:\20231025_HK-BovChon-P6-CMPos-ColVNegative-Rep2'
% };

% % ATDC5 Test
% workDir = "C:\Users\HannahKim\Documents\Sessions\ANALYSIS_HK_ATDC5_H3_BM_HighVsLow";

if ~exist('y','var')
    load handel.mat
end

if ~exist(workDir, 'dir')
    error(['Error: Work directory ' char(workDir) ' is invalid']);
elseif exist(fullfile(workDir,'VoronoiAreaData'),'dir')
    voronoiDir = fullfile(workDir,'VoronoiAreaData');
else
    uigetdir('C:\Users\HannahKim\Documents\Sessions','Select folder with voronoi area data');
end
workDir = fullfile(workDir, 'Analysis');
% workDir = fullfile(workDir, 'Background');
if ~exist(workDir, 'dir')
    mkdir(workDir)
end

disp("Analyzing data...")
disp("Saving data to " + workDir)

plotPercentile = 0;
plotCDF = 0;
processColVI = 0;
plotColVI = 0;
plotLocPerArea = 0;
plotColVIvsVoronoi = 0;

nbins = 500;

% define groups
% groups = {'P0-CMNeg','P1-CMNeg','P2-CMNeg','P6-CMNeg',...
% 'P0-CMPos','P1-CMPos','P2-CMPos','P6-CMPos'};
% groups = {'P0-CMNeg','P2-CMNeg','P6-CMNeg','P0-CMPos','P2-CMPos','P6-CMPos'};
% group_styles = {'--','--','--',...
%     '-','-','-'};
groups = {'LowColVI','HighColVI'};
group_styles = {'-','-'};
group_colors = repmat(lines(length(groups)/2),2,1);

% load data
% voronoi data
if ~exist('voronoiData','var')
    disp("Loading Voronoi data...")
    voronoiData = loadVoronoiData(voronoiDir, groups);
end

%% loc per area plot
if plotLocPerArea
    disp("Plotting localizations vs nuclear area...")
    plotLocPerNuclearArea(workDir, groups, voronoiData)
end

%% plot percentiles
if plotPercentile
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if ~exist('percentileData','var')
        if exist(fullfile(workDir,'percentileData.mat'),'file')
            disp("Loading percentile data...")
            load(fullfile(workDir,'percentileData.mat'));
        else
            disp('Calculating percentile data...')
            percentileData = generatePercentileData(workDir, groups, voronoiData);
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    disp('Plotting percentiles...')
    plotVoronoiDensityPercentiles(workDir, groups, voronoiData, percentileData);
end

%% ColVI plot
if processColVI 
    if ~exist('outputTables','var')
        disp("Loading ColVI data...")
        if exist(fullfile(workDir,'ColVIzStack_outputTable.mat'),'file')
            load(fullfile(workDir,'ColVIzStack_outputTable.mat'),'outputTables');
        else
            disp('Analyzing ColVI data...')
            [outputTables, ~,] = ProcessResults_ColVIZstack(workDir,colIVDataDirList,plotColVI);
        end
    end
end

if plotColVIvsVoronoi
    if ~exist('percentileData','var')
        if exist(fullfile(workDir,'percentileData.mat'),'file')
            disp("Loading percentile data...")
            load(fullfile(workDir,'percentileData.mat'));
        else
            disp('Calculating percentile data...')
            percentileData = generatePercentileData(workDir, groups, voronoiData);
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    disp('Plotting ColIV vs voronoi densities...')
    % all
    logVoronoiDensity_vs_ColVI(workDir, {'P0-CMNeg','P2-CMNeg','P6-CMNeg',...
        'P0-CMPos','P2-CMPos','P6-CMPos'},...
        voronoiData, outputTables, percentileData,'')
    % grouped by treatment
    logVoronoiDensity_vs_ColVI(workDir, {'P0-CMNeg','P2-CMNeg','P6-CMNeg'},...
        voronoiData, outputTables, percentileData,'CMNeg')
    logVoronoiDensity_vs_ColVI(workDir, {'P0-CMPos','P2-CMPos','P6-CMPos'},...
        voronoiData, outputTables, percentileData, 'CMPos')
    % grouped by passage
    logVoronoiDensity_vs_ColVI(workDir, {'P0-CMNeg','P0-CMPos'}, voronoiData, outputTables, percentileData, 'P0')
    % logVoronoiDensity_vs_ColVI(workDir, {'P1-CMNeg','P1-CMPos'}, voronoiData, outputTables, percentileData, 'P1')
    logVoronoiDensity_vs_ColVI(workDir, {'P2-CMNeg','P2-CMPos'}, voronoiData, outputTables, percentileData, 'P2')
    logVoronoiDensity_vs_ColVI(workDir, {'P6-CMNeg','P6-CMPos'}, voronoiData, outputTables, percentileData, 'P6')
    % compare difference between P0-CMPos and P6-CMNeg
    logVoronoiDensity_vs_ColVI(workDir, {'P0-CMPos','P6-CMNeg'}, voronoiData, outputTables, percentileData, 'P0P6')
end
% 
% %% cdf plot
% if plotCDF
%     disp("Plotting Voronoi density CDFs...")
%     groups = {'LowColVI','HighColVI'};
%     group_styles = {'-','-'};
%     group_colors = repmat(lines(length(groups)/2),2,1);
%     plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'',nbins)
%     % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 
%     % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%     % % all
%     % groups = {'P0-CMNeg','P2-CMNeg','P6-CMNeg',...
%     %     'P0-CMPos','P2-CMPos','P6-CMPos'};
%     % group_styles = {'--','--','--',...
%     %     '-','-','-'};
%     % group_colors = repmat(lines(length(groups)/2),2,1);
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'',nbins)
%     % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%     % % grouped by treatment
%     % group_styles = {'-','-','-','-'};
%     % group_colors = lines(length(groups));
%     % groups = {'P0-CMNeg','P2-CMNeg','P6-CMNeg'};
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'CMNeg',nbins)
%     % groups = {'P0-CMPos','P2-CMPos','P6-CMPos'};
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'CMPos',nbins)
%     % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%     % % grouped by passage
%     % group_styles = {'-','-'};
%     % group_colors = lines(length(groups));
%     % groups = {'P0-CMNeg','P0-CMPos'};
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'P0',nbins)
%     % % groups = {'P1-CMNeg','P1-CMPos'};
%     % % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'P1',nbins)
%     % groups = {'P2-CMNeg','P2-CMPos'};
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'P2',nbins)
%     % groups = {'P6-CMNeg','P6-CMPos'};
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'P6',nbins)
%     % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%     % % compare difference between P0-CMPos and P6-CMNeg
%     % groups = {'P0-CMPos','P6-CMNeg'};
%     % group_colors = lines(4);
%     % group_colors = group_colors([1 4],:);
%     % plotVoronoiDensityCDFs(workDir, groups, group_styles, group_colors, voronoiData,'P0P6',nbins)
% end

sound(y)
disp("Done!")
