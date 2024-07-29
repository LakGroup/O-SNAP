% function S = preprocessForClassification(T)
%%
rootDir = 'F:';
analysisName = '\ANALYSIS_AD_BJ_H2B';
reps = {'0'};
groups = {'Ctrl','TSA'};
% %%
% rootDir = 'F:';
% analysisName = 'ANALYSIS_AM_H3-1_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% %%
% rootDir = 'F:';
% analysisName = 'ANALYSIS_AM_H3-2_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% %%
% rootDir = 'F:';
% analysisName = 'ANALYSIS_AM_H3-3_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% %%
% rootDir = 'F:';
% analysisName = 'ANALYSIS_AM_H3-1_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% %%
% rootDir = 'F:';
% analysisName = 'ANALYSIS_AM_H3-2_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% %%
% rootDir = 'F:';
% analysisName = 'ANALYSIS_AM_H3-3_H3';
% reps = {'0'};
% groups = {'Ctrl','D_Ala'};

workDir = fullfile(rootDir,analysisName);
load(fullfile(workDir, analysisName+".mat"),"T");

% find cells from groups of interest
% groups = {''};
rowIdx = cellfun(@(x) strcmpi(T.Group,x), groups,'uni',0);
rowIdx = any([rowIdx{:}],2);
T_norm = T(rowIdx,:);

% normalize and remove unnecessary columns
T_norm(:,4:end) = normalize(T_norm(:,4:end));
T_norm = T_norm(:,[1 4:end]);

save(fullfile(workDir, analysisName+".mat"),"T","T_norm");

% end