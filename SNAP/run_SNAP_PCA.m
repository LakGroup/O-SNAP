function pca_result=run_SNAP_PCA(T,options)
arguments
    T table
    options.vars_selected cell = {}
    options.num_components_explained = 0.75
    options.save_path string = ""
end
%% Apply variable selection and normalize
X =  preprocess_SNAP_table(T,...
    "normalize",true,"remove_NaN",true,...
    "keep_rep_sample_info",true);
if ~isempty(options.vars_selected)
    X = X{:,options.vars_selected};
    pca_result.vars_selected = options.vars_selected;
else
    X = X{:,vartype('numeric')};
    pca_result.vars_selected = T(:,vartype('numeric')).Properties.VariableNames;
end

%% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
top_in_search_path('pca','stats\stats\pca.m');
[pca_result.pca_coefficients, pca_result.pca_scores, ~, ~, pca_result.explained, pca_result.pca_centers] = pca(X);

%% Keep enough components to explain the desired amount of variance.
if options.num_components_explained < 1
    explained_variance_to_keep_as_fraction = options.num_components_explained;
    idx = find(cumsum(pca_result.explained)/sum(pca_result.explained) >= explained_variance_to_keep_as_fraction, 1);
    if idx <= numel(pca_result.explained)
        pca_result.num_comps_to_keep = idx;
    else
         pca_result.num_comps_to_keep = numel(pca_result.explained);
    end
else
    if options.num_components_explained <= numel(pca_result.explained)
        pca_result.num_comps_to_keep = options.num_components_explained;
    else
        pca_result.num_comps_to_keep = numel(pca_result.explained);
    end
end
pca_result.pca_coefficients = pca_result.pca_coefficients(:,1:pca_result.num_comps_to_keep);

pca_result.pca_transformation_fcn = @(x) pca_transformation_fcn(x,pca_result);
pca_result.pca_transformation_fcn(T);

%% Plot
if options.save_path ~= ""
    try
        if exist(options.save_path,'file')
            delete(options.save_path)
        end
        plot_SNAP_PCA(pca_result.pca_coefficients,pca_result.pca_scores,pca_result.explained,"groups_info",T.group,"bio_reps_info",T.biological_replicate,"save_path",options.save_path);
    catch ME
        disp(getReport(ME))
    end
end
end

function x_PCA = pca_transformation_fcn(x,pca_result)
    x_PCA = preprocess_SNAP_table(x,"normalize",true,"remove_NaN",true,"numeric_only",true);
    x_PCA = (x_PCA{:,pca_result.vars_selected} - pca_result.pca_centers) * pca_result.pca_coefficients(:,1:pca_result.num_comps_to_keep);
end