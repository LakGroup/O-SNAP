function [T_train, T_test] = preprocess_for_classification(T, options)
arguments
    T table
    options.test_rep = "";
end

groups = unique(T.group);
bio_reps = unique(T.biological_replicate);
if options.test_rep == ""
    T_train = prepare_voronoi_table_data(T,"groups",cellstr(groups),"replicates",cellstr(bio_reps),"normalize",true);
    T_test = [];
else
    reps_train = bio_reps(~strcmp(bio_reps, options.test_rep));
    T_train = prepare_voronoi_table_data(T,"groups",cellstr(groups),"replicates",cellstr(reps_train),"normalize",true);
    T_test = prepare_voronoi_table_data(T,"groups",cellstr(groups),"replicates",cellstr(options.test_rep),"normalize",true);
end
