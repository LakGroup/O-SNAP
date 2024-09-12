% function S = preprocess_for_classification(T)


T_norm = prepare_voronoi_table_data(T,cellstr(unique(T.group)),cellstr(unique(T.biological_replicate)),"normalize",true);

% save(fullfile(work_dir, analysis_name+".mat"),"T","T_norm");

% end
