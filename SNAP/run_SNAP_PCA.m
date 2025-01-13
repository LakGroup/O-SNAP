function pca_result = run_SNAP_PCA(T,options)
arguments
    T table
    options.vars_selected cell = {}
    options.num_components_explained = 0.95
    options.save_path = ""
end
%% Apply variable selection
X = prepare_voronoi_table_data(T,"numeric_only",true);
if ~isempty(options.vars_selected)
    X = T(:,options.vars_selected);
end
%% Normalize
X = normalize(table2array(X));

%% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
[pca_result.pca_coefficients, pca_result.pca_scores, ~, ~, pca_result.explained, pca_result.pca_centers] = pca(X);

%% Keep enough components to explain the desired amount of variance.
if options.num_components_explained < 1
    explained_variance_to_keep_as_fraction = options.num_components_explained;
    pca_result.num_components_to_keep = find(cumsum(pca_result.explained)/sum(pca_result.explained) >= explained_variance_to_keep_as_fraction, 1);
else
    pca_result.num_components_to_keep = options.num_components_explained;
end
pca_coefficients = pca_result.pca_coefficients(:,1:pca_result.num_components_to_keep);
pca_result.pca_transformation_fcn = @(x) array2table((table2array(varfun(@double, x)) - pca_centers) * pca_coefficients);

%% Plot
if options.save_path ~= ""
    plot_PCA(pca_result.pca_coefficients,pca_result.pca_scores,pca_result.explained,"groups_info",T.group,"bio_reps_info",T.biological_replicate,"save_path",options.save_path);
end
end