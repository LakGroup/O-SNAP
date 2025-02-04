function pca_result = run_SNAP_PCA_v2(T,options)
arguments
    T table
    options.vars_selected cell = {}
    options.num_components_explained = 0.60
    options.save_path string = ""
end
%% Apply variable selection and normalize
T =  preprocess_SNAP_table(T,....
    "normalize",true,"remove_NaN",true,...
    "keep_rep_sample_info",true);
if ~isempty(options.vars_selected)
    X = preprocess_SNAP_table(T(:,options.vars_selected),"numeric_only",true);
else
    X = preprocess_SNAP_table(T,"numeric_only",true);
end
pca_result.vars_selected = X.Properties.VariableNames;
X = table2array(X);

%% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
top_in_search_path('pca','stats\stats\pca.m');
[pca_result.pca_coefficients, pca_result.pca_scores, ~, ~, pca_result.explained, pca_result.pca_centers] = pca(X);

%% Keep enough components to explain the desired amount of variance.
if options.num_components_explained < 1
    explained_variance_to_keep_as_fraction = options.num_components_explained;
    pca_result.num_comps_to_keep = find(cumsum(pca_result.explained)/sum(pca_result.explained) >= explained_variance_to_keep_as_fraction, 1);
else
    pca_result.num_comps_to_keep = options.num_components_explained;
end

pca_result.pca_transformation_fcn = @(x) array2table((table2array(varfun(@double, x)) - pca_result.pca_centers) * pca_result.pca_coefficients(:,1:pca_result.num_comps_to_keep));

%% Plot
if options.save_path ~= ""
    try
        if exist(options.save_path,'file')
            delete(options.save_path)
        end
        plot_PCA(pca_result.pca_coefficients,pca_result.pca_scores,pca_result.explained,"groups_info",T.group,"bio_reps_info",T.biological_replicate,"save_path",options.save_path);
    catch ME
        disp(getReport(ME))
    end
end
end