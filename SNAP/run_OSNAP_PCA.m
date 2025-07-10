% -------------------------------------------------------------------------
% run_OSNAP_PCA.m
% -------------------------------------------------------------------------
% Performs PCA transformation on O-SNAP feature data.
%
% Example on how to use it:
%   run_OSNAP_PCA(feature_data)
% -------------------------------------------------------------------------
% Input:
%   feature_data: The feature data table, where each row represents 
%                 a sample(nucleus). The first three columns represent (1) 
%                 Group/Phenotype, (2) replicate, and (3) Sample 
%                 Identifier. Each subsequent column is  an O-SNAP feature.
% Output:
%   pca_result: A struct array containing info on PCA transformation
%                 - vars_selected: Cell array of the features included in
%                                  the PCA
%                 - pca_coefficients: Coefficient values from PCA
%                 - pca_scores: Scores of the data used to generate PCA
%                 - explained: Variance explained for each component
%                 - pca_centers: Values for PCA centers
%                 - num_comps_to_keep: Number of components kept from PCA
%                 - pca_transformation_fcn: The PCA transformation function
% Options:
%   vars_selected: Cell array of selected features
%   num_components_explained: Minimum cutoff to establish the number of
%                             components to keep from PCA analysis. 
%   save_path: Name of file location to save to.
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function pca_result = run_OSNAP_PCA(feature_data,options)
arguments
    feature_data table
    options.vars_selected cell = {}
    options.num_components_explained = 0.75
    options.save_path string = ""
end
%% Apply variable selection and normalize
X =  preprocess_OSNAP_feature_data(feature_data,...
    "normalize",true,"remove_NaN",true,...
    "keep_rep_sample_info",true);
if ~isempty(options.vars_selected)
    X = X{:,options.vars_selected};
    pca_result.vars_selected = options.vars_selected;
else
    X = X{:,vartype('numeric')};
    pca_result.vars_selected = feature_data(:,vartype('numeric')).Properties.VariableNames;
end

%% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
top_in_search_path_OSNAP('pca','stats\stats\pca.m');
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
pca_result.pca_transformation_fcn(feature_data);

%% Plot
if options.save_path ~= ""
    try
        if exist(options.save_path,'file')
            delete(options.save_path)
        end
        plot_OSNAP_PCA(pca_result.pca_coefficients,pca_result.pca_scores,pca_result.explained,"groups_info",feature_data.group,"bio_reps_info",feature_data.biological_replicate,"save_path",options.save_path);
    catch ME
        disp(getReport(ME))
    end
end
end

function x_PCA = pca_transformation_fcn(x,pca_result)
    x_PCA = preprocess_OSNAP_feature_data(x,"normalize",true,"remove_NaN",true,"numeric_only",true);
    x_PCA = (x_PCA{:,pca_result.vars_selected} - pca_result.pca_centers) * pca_result.pca_coefficients(:,1:pca_result.num_comps_to_keep);
end
