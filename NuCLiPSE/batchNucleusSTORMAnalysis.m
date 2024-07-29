function batchNucleusSTORMAnalysis(workDir, groups, reps)
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
disp("Running analysis for: " + workDir)
disp(['  RUN START: ' char(datetime)])
disp("   GROUPS: " + string(join(groups(:),', ')))
disp("   REPS: " + string(join(reps(:),', ')))
%% analysis parameters
nProcesses = 12;
min_logVorDensity = 0;
max_logVorDensity = 3;
heterochromatin_thresh = 70;
eps = 20;
min_num = 3;
ellipse_inc = 0.1;
peripheryThresh = 0.15;
area_threshold_arr = 10.^(1:0.25:2);
min_number_of_localizations_arr = [5];

%% load data
dataInfoTable = getValidVoronoiAreaData(workDir,groups,reps,{'x','y'});

%% split full data for parallel processing
[splitData, nProcesses] = splitFilesForParallel(dataInfoTable, nProcesses, true);

%% start parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    parpool("Processes",nProcesses);
end

%% calculate voronoi density
disp('Calculating Voronoi densities...');
vorDataVars = {'voronoi_areas_all','voronoi_neighbors','voronoi_areas','reduced_log_voronoi_density'};
tic
parfor p=1:nProcesses
    dataInfoTable_p = splitData{p};
    filepaths = dataInfoTable_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file') && ~hasVars(filepath,vorDataVars)
            try
                generateVoronoiData(filepath, min_logVorDensity, max_logVorDensity);
            catch ME
                disp(getReport(ME));
                disp("Error in " + filepath + ", removing from analysis")
                splitData{p}(s,:) = [];
            end
        end
    end
end
toc
disp(['Completed: ' char(datetime)])
%% calculate density threshold
densityData=cell(1,nProcesses);
if ~exist('density_threshold','var')
    disp('Calculating density threshold...');
    tic
    parfor p=1:nProcesses
        dataInfoTable_p = splitData{p};
        filepaths = dataInfoTable_p{:,'filepath'};
        densityData{p} = cell(1,length(filepaths));
        for s=1:length(filepaths)
            sample = load(filepaths(s),'voronoi_areas');
            densityData{p}{s} = 1./sample.voronoi_areas;
        end
        densityData{p} = vertcat(densityData{p}{:});
    end
    density_threshold = prctile(vertcat(densityData{:}),heterochromatin_thresh);
    clearvars densityData
    toc
end
disp(['Completed: ' char(datetime)])
%% perform heterochromatin analysis
disp("Performing Heterochromatin Analysis...")
heteroDataVars = {'nucleus_radius',...
            'locs_number',...
            'locs_density',...
            'boundary',... 
            'Hetero_flt_with_label',...
            'lads',...
            'non_lads',...
            'lads2total',...
            'DomainBoundary',...
            'hetero_center',...
            'hetero_radius',...
            'hetero_nLocs',...
            'hetero_density',...
            'lads_area',...
            'lads_nLocs',...
            'lads_density',...
            'spacing',...
            'seg_normals',...
            'seg_locs',...
            'lads_seg_len',...
            'lads_seg_thickness',...
            'components',...
            'Eigenvalues',...
            'locs_norm',...
            'x_length',...
            'y_length',...
            'polygon',...
            'ElasticEnergy',...
            'BendingEnergy',...
            'BorderCurvature',...
            'radialDensity',...
            'radialHeteroDensity',...
            'interiorDensity',...
            'peripheryDensity',...
            'interiorHeteroDensity',...
            'peripheryHeteroDensity'};
tic
parfor p=1:nProcesses
    dataInfoTable_p = splitData{p};
    filepaths = dataInfoTable_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file') && ~hasVars(filepath,heteroDataVars)
            try
                Nucleus_STORM_Analysis_MATLAB_v4_HK(filepath,density_threshold,eps,min_num,ellipse_inc,peripheryThresh);
            catch ME
                disp(getReport(ME));
                disp("Removing from analysis: " + filepath)
                splitData{p}(s,:) = [];
            end
        end
    end
end
toc
disp(['Completed: ' char(datetime)])
%% perform voronoi clustering analysis
disp("Performing Voronoi Clustering Analysis...")
clusterDataVars1 = {'area_thresholds','min_number_of_localizations','clusters'};
clusterDataVars2 = {'clusterNLocs','clusterArea','clusterDensity','clusterGyrationR'};
tic
parfor p=1:nProcesses
    dataInfoTable_p = splitData{p};
    filepaths = dataInfoTable_p{:,'filepath'};
    for s=1:length(filepaths)
        filepath = filepaths(s);
        if exist(filepath,'file')
            if ~hasVars(filepath,clusterDataVars1)
                try
                    generateVoronoiClusters(filepath,area_threshold_arr,min_number_of_localizations_arr);
                catch ME
                    disp(getReport(ME));
                    disp("Removing from analysis: " + filepath)
                    splitData{p}(s,:) = [];
                end
            end
            if ~hasVars(filepath,clusterDataVars2) && hasVars(filepath,clusterDataVars1)
                try
                    extractClusterFeatures(filepath);
                catch ME
                    disp(getReport(ME));
                    disp("Removing from analysis: " + filepath)
                    splitData{p}(s,:) = [];
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