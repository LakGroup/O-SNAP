function nSampPerProcess = getNSampPerProcess(nSamps,nProcesses)
    nSampPerProcess_nMin1 = floor(nSamps/nProcesses);
    nSampPerProcess_n = nSamps - nSampPerProcess_nMin1*(nProcesses);
    nSampPerProcess = [repmat(nSampPerProcess_nMin1, 1, nProcesses)];
    nSampPerProcess(1:nSampPerProcess_n) = nSampPerProcess(1:nSampPerProcess_n) + 1;
end