function T = getValidVoronoiAreaData(workDir,groups,reps,varsToLoad)
    if(exist(workDir,"dir"))
        folderInfo = dir(workDir);
        folderInfo = folderInfo([folderInfo(:).isdir]);
        dirNames = {folderInfo(3:end).name};
        repsName = cellfun(@(x) ['Rep' x], reps,'uni',0);
        validDirsIdx = cellfun(@(x) contains(x, repsName),dirNames,'uni',1);
        validDirs = dirNames(validDirsIdx);
        %% get preliminary information on valid files to load based on input
        n_samps = 0;
        dirNames = [];
        sampleNames = [];
        repNames = [];
        groupNames = [];
        for d=1:length(validDirs)
            validDir = string(fullfile(workDir,validDirs(d),'VoronoiAreaData'));
            folderInfo = dir(validDir);
            folderInfo = folderInfo(~[folderInfo(:).isdir]);
            fileNames = {folderInfo(:).name};
            validIdx = all([endsWith(fileNames,".mat"); contains(fileNames, 'storm','IgnoreCase',true); contains(fileNames, groups,'IgnoreCase',true)]);
            nValidFiles = sum(validIdx);
            sampleNames_d = string(extractBefore(fileNames(validIdx),".mat"))';
            dirNames = [dirNames; repmat(validDir,nValidFiles,1)];
            sampleNames = [sampleNames; sampleNames_d];
            repNames = [repNames; repmat(string(reps{d}),nValidFiles,1)];
            groupNames = [groupNames; arrayfun(@(x) findGroup(x,groups), sampleNames_d,'uni',1)];
            n_samps = n_samps + nValidFiles;
        end
        filePaths = arrayfun(@(x,y) fullfile(x,y+".mat"),dirNames,sampleNames,'uni',1);
        T = table(dirNames,sampleNames,repNames,groupNames,filePaths);
        T.Properties.VariableNames={'dir','name','rep','group','filepath'};
        T = T(idxVoronoiAreaDataContainingVar(T,varsToLoad),:);
    else
        disp("No valid directory: " + workDir)
        T = [];
    end
end

function group = findGroup(sampleName,groups)
    groupIdx = cell2mat(cellfun(@(x) contains(sampleName, ['-' x]), groups,'uni',0));
    if any(groupIdx)
        group = string(groups{groupIdx});
    else
        groupIdx = cell2mat(cellfun(@(x) contains(sampleName, x), groups,'uni',0));
        if any(groupIdx)
            group = string(groups{groupIdx});
        else
            error('Error: At least one group specified could not be found in the data set')
        end
    end
end