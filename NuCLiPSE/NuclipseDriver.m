% prepare compiler for faster run
% if ~exist('firstRun','var')
%     nuclei_classification_setup
%     firstRun = true;
% end
% %%

rootDir = 'F:\';

% analysisName = 'ANALYSIS_Simulation';
% reps = {'0'};
% groups = {'Group1','Group2'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AD_BJ_H2B';
% reps = {'0'};
% groups = {'Ctrl','TSA'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AD_BJ_EdC';
% reps = {'0'};
% groups = {'Ctrl','TSA'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_PC_BJ_H2B';
% reps = {'0'};
% groups = {'Ctrl','TSA'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_H3-1_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_H3-2_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_H3-3_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_H3-1_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_H3-2_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_H3-3_H3';
% reps = {'0'};
% groups = {'Ctrl','D_Ala'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_EK_HCT116_H2B';
% reps = {'0'};
% groups = {'ntc2','HNRNPK2'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_HK_BJ_H2B';
% % reps = {'0-0','0-1','1','2'};
% reps = {'0-1','1','2'};
% groups = {'Ctrl','LMNAKD'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_AM_Heterokaryon_H2B';
% reps = {'0','1','2','3'};
% % groups = {'CtrlhFb','CtrlmESC','HK06h','HK12h','HK24h','HK36h','HK48h'};
% groups = {'CtrlhFb','CtrlmESC','HK48h'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_SH_hTen_H2B';
% reps = {'0'};
% groups = {'Healthy','Tendinosis'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_SH_hMSC_H2B_substrate';
% reps = {'0','1'};
% groups = {'SoftMeHA','StiffMeHA','Glass'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_SH_hMSC_H2B_stiffeninggel';
% reps = {'0'};
% groups = {'00h','02h','04h','06h','18h','24h','48h'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_SH_hTC_H2B_substrate';
% reps = {'0'};
% groups = {'SoftMeHA','StiffMeHA','Glass'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_SH_hTC_H2B_stiffeninggel';
% reps = {'0'};
% groups = {'Healthy-04h','Healthy-24h','Old-04h','Old-24h','Tendinosis-04h','Tendinosis-24h'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_KS_CM_H3K9me2';
% reps = {'0','1'};
% % groups = {'Control1','LMNAKD','LMNADNK','DNK'};
% groups = {'Control1','LMNAKD'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_JZ_Heterokaryon_H1.0';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_JZ_Heterokaryon_H1.4';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_JZ_Heterokaryon_H1.5';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_JZ_Heterokaryon_H2B';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos','HK7d-SOX2Neg','HK7d-SOX2Pos'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_JZ_Heterokaryon_H2B';
% reps = {'0','1'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos','HK7d-SOX2Neg','HK7d-SOX2Pos'};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %%
% analysisName = 'ANALYSIS_HK_ATDC5_H3_CMPos_YAPTAZ';
% reps = {'0','1'};
% groups = {''};
% runNuclipseAnalysis(rootDir,analysisName,groups,reps);
% %% wrap up
% delete(pPool)
disp('Done!')
%%
function T = runNuclipseAnalysis(rootDir,analysisName,groups,reps)
    warning('off','all');
    workDir = fullfile(rootDir,analysisName);
    tableFilePath = fullfile(workDir,analysisName+".mat");
    if ~exist(tableFilePath,"file")
        %% calculate features
        batchNucleusSTORMAnalysis(workDir,groups,reps);
        %% coallate features
        dataInfoTable = getValidVoronoiAreaData(workDir,groups,reps,{'nucleus_radius','clusterNLocs'});
        T = voronoiDataToTable_parallel(dataInfoTable);
    else
        load(tableFilePath,"T");
    end
    T = T(all([ismember(T.BioRep,reps), ismember(T.Group,groups)],2),:);
    %% generate plots
    disp('Plotting variable selection and PCA...')
    if length(unique(T.Group)) >= 2
        plotNuclipse(T, groups, workDir);
    else
        disp('Not enough groups in table')
    end
    % plotRadialDensities(dataInfoTable, T)
    % ripleyK(workDir,groups,reps);
    warning('on','all');
end