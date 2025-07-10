% -------------------------------------------------------------------------
% run_OSNAP_classification_batch.m
% -------------------------------------------------------------------------
% Trains multiple sets of classifiers across folds and, if specified,
% multiple instances of models with the same architecture. Implements
% parallel processing for improved runtime.
%
% Example on how to use it:
%   run_OSNAP_classification_batch(feature_data)
% -------------------------------------------------------------------------
% Input:
%   train_data: The feature data of the training set, where each row  
%               represents a sample(nucleus). The first three columns 
%               represent (1) Group/Phenotype, (2) replicate, and (3)  
%               SampleIdentifier. Each subsequent column is  an O-SNAP 
%               feature.
% Output:
%   classifiers: A cell array where each cell stores a single instance of a
%                model trained on the O-SNAP feature data
% Options:
%   test_data: A set of O-SNAP feature data for test samples. If no test
%              data is provided, the classifier will undergo
%              cross-validation using feature_data
%   pca_result: A struct array containing info on PCA transformation
%   vars_selected: Cell array of selected features
%   n_models_per_type: Number of model instances to train within each fold
%                      for a given model architecture. Incorporated to
%                      evaluate the stability of a model architecture
%                      trained on the same data.
%   n_processes: Number of cores to run parallel processes on 
%   verbose: Flag for output to print
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%
function classifiers = run_OSNAP_classification_batch(train_data,options)
arguments
    train_data table
    options.test_data = []
    options.pca_result struct = []
    options.vars_selected cell = {}
    options.n_models_per_type double = 5
    options.n_processes double = 12
    options.verbose logical = false
end
groups = unique(train_data.group);
model_types = [...
    "tree_fine",...
    "tree_medium",...
    "tree_coarse",...
    "discriminant_linear",...
    "discriminant_quadratic",...
    "logistic_regression_efficient",...
    "svm_efficient_linear",...
    "naive_bayes_gaussian",...
    "naive_bayes_kernel",...
    "svm_linear",...
    "svm_quadratic",...
    "svm_gaussian_fine",...
    "svm_gaussian_medium",...
    "svm_gaussian_coarse",...
    "knn_fine",...
    "knn_medium",...
    "knn_coarse",...
    "knn_cosine",...
    "knn_cubic",...
    "knn_weighted",...
    "ensemble_boosted_trees",...
    "ensemble_bagged_trees",...
    "ensemble_subspace_discriminant",...
    "ensemble_rus_boosted_trees",...
    "NN_narrow",...
    "NN_medium",...
    "NN_wide",...
    "NN_bilayered",...
    "NN_trilayered",...
    "kernel_svm",...
    "kernel_logistic_regression"];
if numel(groups) == 2
    model_types = [model_types "logistic_regression_binary_glm"];
end

n_models_per_type = options.n_models_per_type;
pca_result = options.pca_result;
vars_selected = options.vars_selected;
test_data = options.test_data;
verbose = options.verbose;
[model_types_p,n_processes] = split_data_to_n_OSNAP(model_types,options.n_processes,"shuffle",false);
classifiers_p = cell(1,n_processes);
% generate classifiers of given model types
parfor p=1:n_processes
    warning('off','all')
    n_model_types = numel(model_types_p{p});
    classifiers_p{p} = cell(n_models_per_type,n_model_types);
    for j=1:n_model_types
        try
            for i=1:n_models_per_type
                classifiers_p{p}{i,j} = run_OSNAP_classification(...
                    train_data,...
                    'pca_result',pca_result,...
                    'vars_selected',vars_selected,...
                    'test_data',test_data,...
                    'model_type',model_types_p{p}(j),...
                    'verbose',false);
            end
        catch ME
            if verbose
                disp("Error with model type: " + model_types_p{p}(j))
                disp(getReport(ME))
            end
            for i=1:n_models_per_type
                classifiers_p{p}{i,j} = [];
            end
        end
    end
    warning('on','all')
end
classifiers = [classifiers_p{:}];

% filter models that threw errors
valid_idx = true(size(classifiers,2),1);
for j=1:numel(model_types)
    for i=1:options.n_models_per_type
        if isempty(classifiers{i,j})
            valid_idx(j) = false;
        end
    end
end
classifiers = classifiers(:,valid_idx);
end
