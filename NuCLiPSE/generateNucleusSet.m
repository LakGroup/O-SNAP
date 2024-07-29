% Generate multiple nuclei
function simData = generateNucleusSet(saveDir,n,options)
    arguments
        saveDir string;
        n double = 1;
        options.majorAxis (1,2) double = [2e4,2.5e3];
        options.minorAxis (1,2) double = [1.3e4,2e3];
        options.intPeripheryCutoff (1,2) double = [0.85,0.05];
        options.densityInterior (1,2) double = [7.4e-3,2e-3];
        options.densityPeriphery (1,2) double = [7.5e-3,2.3e-3];
        options.densityHeteroCenter (1,2) double = [1e-5,3e-6];
        options.densityHetero (1,2) double = [2.3e-2,4.8e-3];
        options.radiusHetero (1,2) double = [30,4.5];
        options.plotFlag logical = false;
        options.nProcesses double = 12;
    end
    if options.plotFlag
        figure();
        tiledlayout('flow');
    end
    % reduce broadcast overhead
    majorAxis = repmat(options.majorAxis,options.nProcesses,1);
    majorAxis_mean = majorAxis(:,1);
    majorAxis_std = majorAxis(:,2);
    minorAxis = repmat(options.minorAxis,options.nProcesses,1);
    minorAxis_mean = minorAxis(:,1);
    minorAxis_std = minorAxis(:,2);
    intPeripheryCutoff = repmat(options.intPeripheryCutoff,options.nProcesses,1);
    intPeripheryCutoff_mean = intPeripheryCutoff(:,1);
    intPeripheryCutoff_std = intPeripheryCutoff(:,2);
    densityInterior = repmat(options.densityInterior,options.nProcesses,1);
    densityInterior_mean = densityInterior(:,1);
    densityInterior_std = densityInterior(:,2);
    densityPeriphery = repmat(options.densityPeriphery,options.nProcesses,1);
    densityPeriphery_mean = densityPeriphery(:,1);
    densityPeriphery_std = densityPeriphery(:,2);
    densityHeteroCenter = repmat(options.densityHeteroCenter,options.nProcesses,1);
    densityHeteroCenter_mean = densityHeteroCenter(:,1);
    densityHeteroCenter_std = densityHeteroCenter(:,2);
    densityHetero = repmat(options.densityHetero,options.nProcesses,1);
    densityHetero_mean = densityHetero(:,1);
    densityHetero_std = densityHetero(:,2);
    radiusHetero = repmat(options.radiusHetero,options.nProcesses,1);
    radiusHetero_mean = radiusHetero(:,1);
    radiusHetero_std = radiusHetero(:,2);
    plotFlag = repmat(options.plotFlag,options.nProcesses,1);
    % split data for parallel
    [simData, nProcesses] = splitDataForParallel(cell(1,n),...
        options.nProcesses, false);
    idx = splitDataForParallel(num2cell(1:n), options.nProcesses,false);
    % create nuclei
    tic
    parfor p=1:nProcesses
        for s=1:length(simData{p})
            % create random values based on options
            while 1
                majorAxis_val = random('Normal',majorAxis_mean(p),majorAxis_std(p));
                minorAxis_val = random('Normal',minorAxis_mean(p),minorAxis_std(p));
                intPeripheryCutoff_val = random('Normal',intPeripheryCutoff_mean(p),intPeripheryCutoff_std(p));
                densityInterior_val = random('Normal',densityInterior_mean(p),densityInterior_std(p));
                densityPeriphery_val = random('Normal',densityPeriphery_mean(p),densityPeriphery_std(p));
                densityHeteroCenter_val = random('Normal',densityHeteroCenter_mean(p),densityHeteroCenter_std(p));
                densityHetero_val = random('Normal',densityHetero_mean(p),densityHetero_std(p)) - densityInterior_mean(p);
                radiusHetero_val = random('Normal',radiusHetero_mean(p),radiusHetero_std(p)); % random avg cluster size per nucleus
                if any([majorAxis_val, minorAxis_val,...
                    intPeripheryCutoff_val, densityInterior_val,...
                    densityPeriphery_val, densityHeteroCenter_val,...
                    densityHetero_val, radiusHetero_val])
                    break
                end
            end
            angle_val = 360*rand();
            if plotFlag(p)
                nexttile
            end
            points = generateNucleus([majorAxis_val minorAxis_val],[0 0],angle_val,intPeripheryCutoff_val,...
            densityInterior_val,densityPeriphery_val,densityHeteroCenter_val,densityHetero_val,radiusHetero_val,...
            "plotFlag",plotFlag(p));
            simData{p}{s}.x_data = points(:,1);
            simData{p}{s}.y_data = points(:,2);
            %% save localization data for subsequent analysis
            simData{p}{s} = generateLocData(simData{p}{s}, 1);
            simData{p}{s}.name = ['simulatedNucleus-Group1-storm-' sprintf('%02.0f',idx{p}{s})];
            saveVoronoiAnalysisData(fullfile(saveDir,[simData{p}{s}.name '.mat']), simData{p}{s});
            % disp("   " + simData{p}{s}.name + ": " + string(datetime))
        end
    end
    simData = [simData{:}];
    toc
end