function T = voronoiDataToTable_parallel(dataInfoTable)

    %% parameters   
    dirPathSplit = regexp(dataInfoTable{1,1},'\','split');
    rootDir = join(dirPathSplit(1:end-3),'\');
    analysisName = dirPathSplit(end-2);
    nProcesses = 12;
    workDir = fullfile(rootDir,analysisName);
    varsString = {'Group','BioRep','name'};

    %% organize into tables
    if(~exist("T","var"))
        if exist(fullfile(workDir, analysisName+".mat"),'file')
            disp("Loading " + fullfile(workDir, analysisName+".mat"))
            load(fullfile(workDir, analysisName+".mat"),"T")
        else
            tic
            % split data for parallel processing
            [splitData, nProcesses] = splitFilesForParallel(dataInfoTable, nProcesses, true);
            % start parallel pool
            pPool = gcp('nocreate');
            disp('Generating table data...')
            tic
            if isempty(pPool)
                parpool("Processes",nProcesses);
            end
            T_cell = cell(1,nProcesses);
            % run analysis
            for p=1:nProcesses
                dataInfoTable_p = splitData{p};
                T_cell{p} = voronoiDataToTable(dataInfoTable_p);
            end
            T = sortrows(vertcat(T_cell{:}),varsString);
            disp("Saving to " + fullfile(workDir, analysisName+".mat") + "...")
            save(fullfile(workDir, analysisName+".mat"),"T");
        end
    end

end
