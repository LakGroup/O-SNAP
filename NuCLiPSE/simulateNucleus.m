saveDir = uigetdir('H:\','Select analysis directory');

if ~exist(fullfile(saveDir,"SimData1"),'dir')
    mkdir(fullfile(saveDir,"SimData1"))
end
if ~exist(fullfile(saveDir,"SimData2"),'dir')
    mkdir(fullfile(saveDir,"SimData2"))
end

generateNucleusSet(fullfile(saveDir,"SimData1"),80);
% generateNucleusSet(fullfile(saveDir,"SimData2"),80,"majorAxis",[2e3,2.4e2],"minorAxis",[1.3e3,2e2]);

% batchGenerateLocs(data, workDir)