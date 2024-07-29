function [splitData, nProcesses] = splitFilesForParallel(dataInfoTable, nProcesses, shuffle)
    arguments
        dataInfoTable table
        nProcesses double
        shuffle logical = true
    end
    nSamps = size(dataInfoTable,1);
    dataInfoTable_shuffled = dataInfoTable;
    if shuffle
        dataInfoTable_shuffled = dataInfoTable_shuffled(randperm(nSamps),:);
    end
    % if the number of samples is less than parallel processing
    if nSamps < nProcesses
        disp('Number of samples is less than number of processes; not very efficient for parallel...')
        nProcesses = nSamps;
        splitData = cell(1,nSamps);
        for i=1:nSamps
            splitData{i} = dataInfoTable_shuffled(i,:);
        end
    end
    % else, divide into groups for parallel processing
    nSampPerProcess = getNSampPerProcess(nSamps,nProcesses);
    splitData = cell(1,nProcesses);
    idx = 1;
    for i=1:nProcesses
        splitData{i} = dataInfoTable_shuffled(idx:idx+nSampPerProcess(i)-1,:);
        idx = idx + nSampPerProcess(i);
    end
end