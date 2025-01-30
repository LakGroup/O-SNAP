function [classifiers,classification_summary] = run_SNAP_classification_batch(response,pca_result,vars_selected,options)
arguments
    response string
    pca_result struct
    vars_selected cell
    options.num_trials = 20;
    options.num_components double = 0.95
    options.compute_validation_accuracy logical = true;
    options.k_folds double = 5;
    options.verbose logical = false;
end
groups = unique(response);
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
classifiers = cell(n_model_types,options.num_trials);
num_components = options.num_components;
compute_validation_accuracy = options.compute_validation_accuracy;
k_folds = options.k_folds;
num_trials = options.num_trials;
verbose = options.verbose;
parfor i=1:numel(model_types)
    warning('off','all')
    for j=1:num_trials
        try
            classifiers{i,j} = run_SNAP_classification(...
                response,pca_result,vars_selected,...
                'num_components',num_components,...
                'compute_validation_accuracy',compute_validation_accuracy,...
                'model_type',model_types(i),...
                'k_folds',k_folds,...
                'verbose',false);
        catch ME
            if verbose
                disp("Error with model type: " + model_types(i))
                disp(getReport(ME))
            end
            classifiers{i,j} = [];
            break
        end
    end
    warning('on','all')
end

% filter out models that threw errors
valid_idx = true(n_model_types,1);
for i=1:numel(model_types)
    if any(cellfun('isempty',classifiers(i,:)))
        valid_idx(i) = false;
    end
end
classifiers= classifiers(valid_idx,:);
n_model_types = size(classifiers,1);

% sort by validation accuracies
validation_accuracies_mean = zeros(1,n_model_types);
validation_accuracies_std = zeros(1,n_model_types);
for i=1:n_model_types
    validation_accuracies_mean(i) = mean(cellfun(@(x) x.ValidationResults.ValidationAccuracy,classifiers(i,:),'uni',1));
    validation_accuracies_std(i) = std(cellfun(@(x) x.ValidationResults.ValidationAccuracy,classifiers(i,:),'uni',1));
end
[validation_accuracies_mean_sort,validation_accuracies_sort_idx]=sort(validation_accuracies_mean,'descend');
validation_accuracies_std_sort = validation_accuracies_std(validation_accuracies_sort_idx);
classifiers = classifiers(validation_accuracies_sort_idx,:);
if options.verbose
    for i=1:n_model_types
        fprintf('%34s: %6.3f +/- %3.3f\n',classifiers{i}.ModelType,validation_accuracies_mean_sort(i)*100,validation_accuracies_std_sort(i)*100);
    end
end
classifier_type = string(cellfun(@(x) x.ModelType, classifiers(:,1),'uni',0));
validation_accuracies_mean_sort = validation_accuracies_mean_sort';
validation_accuracies_std_sort = validation_accuracies_std_sort';

classification_summary = table(classifier_type,validation_accuracies_mean_sort,validation_accuracies_std_sort);
end