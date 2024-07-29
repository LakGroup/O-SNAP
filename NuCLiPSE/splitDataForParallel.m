function [splitData, nProcesses] = splitDataForParallel(data, nProcesses, shuffle)
    arguments
        data cell
        nProcesses double = 12
        shuffle logical = true
    end
    nSamps = size(data,2);
    data_shuffled = data;
    if shuffle
        data_shuffled = data_shuffled(randperm(nSamps));
    end
    % if the number of samples is less than parallel processing
    if nSamps < nProcesses
        disp('Number of samples is less than number of processes; not very efficient for parallel...')
        nProcesses = nSamps;
        splitData = cell(1,nSamps);
        for i=1:nSamps
            splitData{i} = data(i);
        end
    end
    % else, divide into groups for parallel processing
    nSampPerProcess = getNSampPerProcess(nSamps,nProcesses);
    splitData = cell(1,nProcesses);
    idx = 1;
    for i=1:nProcesses
        splitData{i} = data_shuffled(idx:idx+nSampPerProcess(i)-1);
        idx = idx + nSampPerProcess(i);
    end
end