% %% prepare compiler for faster run
% if ~exist('first_run','var')
%     setup_SNAP
%     first_run = true;
% end
%% start parallel pool
p_pool = gcp('nocreate');
if isempty(p_pool)
    parpool("Processes");
elseif class(p_pool) == "parallel.ThreadPool"
    delete(p_pool)
    parpool("Processes");
end

root_dir = 'F:\';
% %%
% analysis_name = 'ANALYSIS_AD_BJ_EdC';
% reps = {'0'};
% groups = {'Ctrl','TSA'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %% 
% analysis_name = 'ANALYSIS_AD_BJ_H2B';
% reps = {'0'};
% groups = {'Ctrl','TSA'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_PC_BJ_H2B';
% reps = {'0'};
% groups = {'Ctrl','TSA'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %
% analysis_name = 'ANALYSIS_YZ_hChon_H2B';
% reps = {'A46yoM','B55yoM'};
% groups = {'P0','P3','P6'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
analysis_name = 'ANALYSIS_HK_BJ_TSA_H2B';
% reps = {'20240920','20240930','20241010','20241011','20241018'};
reps = {'20240930','20241010','20241011','20241018'};
groups = {'Ctrl','TSA'};
run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_HK_BJ_LMNAKD_H2B';
% % reps = {'20240524_0-0','20240604_0 -1','20240605_1','20240612_2','20240925_3','20241002_4','20241011_5'};
% % reps = {'20240604_0-1','20240605_1','20240612_2','20240925_3','20241002_4','20241011_5','20241023_P9_6','20241023_P10_7'};
% reps = {'20240925_3','20241002_4','20241011_5','20241023_P9_6','20241023_P10_7'};
% groups = {'Ctrl','LMNAKD'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H2B';
% reps = {'0','1','2','3'};
% groups = {'CtrlhFb','CtrlmESC','HK06h','HK12h','HK24h','HK36h','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K4me3';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'CtrlhFb','CtrlmESC','HK06h','HK24h','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K9ac';
% reps = {'20201118','20201203','20201209'};
% groups = {'CtrlhFb','CtrlmESC','HK06h','HK24h','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K9me3';
% reps = {'20201118','20201203','20201209','20210520'};
% groups = {'CtrlhFb','CtrlmESC','HK06h','HK24h','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K27me3';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'CtrlhFb','CtrlmESC','HK06h','HK24h','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-1_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-2_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-3_H2B';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-1_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-2_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-3_H3';
% reps = {'0'};
% groups = {'Control','D_Ala'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% HistoneMarksGroup()
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H2B_CtrlhFb_CtrlmESC_HK48h';
% reps = {'0','1','2','3'};
% groups = {'CtrlhFb','CtrlmESC','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K4me3_CtrlhFb_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'CtrlhFb','CtrlmESC','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K9ac_CtrlhFb_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'CtrlhFb','CtrlmESC','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K9me3_CtrlhFb_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'CtrlhFb','CtrlmESC','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K27me3_CtrlhFb_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'CtrlhFb','CtrlmESC','HK48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K9RepressiveActive_CtrlhFb';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'H3K9ac','H3K9me3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_H3K9RepressiveActive_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'H3K9ac','H3K9me3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_Active_CtrlhFb';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'H3K4me3','H3K9ac'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_Active_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'H3K4me3','H3K9ac'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_Repressive_CtrlhFb';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'H3K9me3','H3K27me3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_Heterokaryon_Repressive_HK48h';
% reps = {'20201118','20201203','20201209','20210824'};
% groups = {'H3K9me3','H3K27me3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% HistoneVariantRename()
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H2B_Control';
% reps = {'0'};
% groups = {'H3-1','H3-2','H3-3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H2B_D_Ala';
% reps = {'0'};
% groups = {'H3-1','H3-2','H3-3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3_Control';
% reps = {'0'};
% groups = {'H3-1','H3-2','H3-3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3_D_Ala';
% reps = {'0'};
% groups = {'H3-1','H3-2','H3-3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-1_Control';
% reps = {'0'};
% groups = {'H2B','H3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-2_Control';
% reps = {'0'};
% groups = {'H2B','H3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-3_Control';
% reps = {'0'};
% groups = {'H2B','H3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-1_D_Ala';
% reps = {'0'};
% groups = {'H2B','H3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-2_D_Ala';
% reps = {'0'};
% groups = {'H2B','H3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_AM_H3Variant_H3-3_D_Ala';
% reps = {'0'};
% groups = {'H2B','H3'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_SH_hTC_H2B';
% reps = {'0'};
% groups = {'Healthy','Tendinosis'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_SH_hMSC_H2B_substrate';
% reps = {'0','1'};
% groups = {'Soft_MeHA','Stiff_MeHA','Glass'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_SH_hMSC_H2B_stiffeninggel';
% reps = {'0'};
% groups = {'00h','02h','04h','06h','18h','24h','48h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_SH_hTC_H2B_substrate';
% reps = {'0'};
% groups = {'Soft_MeHA','Stiff_MeHA','Glass'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_SH_hTC_H2B_stiffeninggel';
% reps = {'0'};
% groups = {'Healthy-04h','Healthy-24h','Old-04h','Old-24h','Tendinosis-04h','Tendinosis-24h'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_EK_HCT116_H2B';
% reps = {'0'};
% groups = {'ntc2','HNRNPK2'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_KS_CM_H3K9me2';
% reps = {'0','1'};
% % groups = {'Control1','LMNAKD','LMNADNK','DNK'};
% groups = {'Control1','LMNAKD'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_JZ_Heterokaryon_H1.0';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_JZ_Heterokaryon_H1.4';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_JZ_Heterokaryon_H1.5';
% reps = {'0'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_JZ_Heterokaryon_H2B';
% reps = {'0','1'};
% groups = {'hFbCtrl','HK4d-SOX2Neg','HK4d-SOX2Pos','HK7d-SOX2Neg','HK7d-SOX2Pos'};
% run_SNAP(root_dir,analysis_name,groups,reps);
% %%
% analysis_name = 'ANALYSIS_HK_ATDC5_H3_CMPos_YAPTAZ';
% % reps = {'0','1'};
% reps = {'0'};
% groups = {''};
% run_SNAP(root_dir,analysis_name,groups,reps);


% %% wrap up
% delete(p_pool)
% disp('Done!')
% diary off 