analysis_name = '';

pca_path = fullfile(work_dir, join(['PCA',groups,options.suffix],'_')+".fig");
volcano_path = fullfile(save_dir, "features_volcano_"+join(group_pairs(g,:),'_')+".fig");