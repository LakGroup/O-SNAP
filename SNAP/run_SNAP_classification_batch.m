function classifiers = run_SNAP_classification_batch(train_data,options)
arguments
    train_data table
    options.test_data = []
    options.pca_result struct = []
    options.vars_selected cell = []
    options.verbose logical = false
end
groups = unique(train_data.group);
model_types = [...
    "tree_fine","tree_medium","tree_coarse",...
    "discriminant_linear","discriminant_quadratic",...
    "logistic_regression_efficient","svm_efficient_linear",...
    "naive_bayes_gaussian","naive_bayes_kernel",...
    "svm_linear","svm_quadratic",...
    "svm_gaussian_fine","svm_gaussian_medium","svm_gaussian_coarse",...
    "knn_fine","knn_medium","knn_coarse","knn_cosine","knn_cubic","knn_weighted",...
    "ensemble_boosted_trees","ensemble_bagged_trees","ensemble_subspace_discriminant","ensemble_rus_boosted_trees",...
    "NN_narrow","NN_medium","NN_wide","NN_bilayered","NN_trilayered",...
    "kernel_svm","kernel_logistic_regression"];
if numel(groups) == 2
    model_types = [model_types "logistic_regression_binary_glm"];
end
n_model_types = numel(model_types);

% generate classifiers of given model types
classifiers = cell(1,n_model_types);
for i=1:numel(model_types)
    warning('off','all')
    try
        classifiers{i} = run_SNAP_classification(...
            train_data,...
            'pca_result',options.pca_result,...
            'vars_selected',options.vars_selected,...
            'test_data',options.test_data,...
            'model_type',model_types(i),...
            'verbose',false);
    catch ME
        if options.verbose
            disp("Error with model type: " + model_types(i))
            disp(getReport(ME))
        end
        classifiers{i} = [];
    end
    warning('on','all')
end

% filter models that threw errors
valid_idx = true(n_model_types,1);
for i=1:numel(model_types)
    if isempty(classifiers(i))
        valid_idx(i) = false;
    end
end
classifiers = classifiers(~cellfun(@isempty, classifiers));
end