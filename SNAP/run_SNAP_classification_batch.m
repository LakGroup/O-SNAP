load("F:\ANALYSIS_KS_CM_H3K9me2\ANALYSIS_KS_CM_H3K9me2","T")
reps = {'0','1'};
groups = {'Ctrl','LMNAKD'};
T = T(all([ismember(T.biological_replicate,reps), ismember(T.group,groups)],2),:);

[~,summary] = run_SNAP_classification_batchX(T,"num_trials",1);

function [classifiers,summary] = run_SNAP_classification_batchX(T,options)
arguments
    T table
    options.method_feature_selection string = 'fscmrmr'
    options.num_trials = 10;
    options.num_features_selected double = 3;
    options.num_components_explained double = 0.95
    options.compute_validation_accuracy logical = true;
    options.k_folds double = 5;
    options.verbose logical = false;
end

groups = unique(T.group);

valid_model_types = [...
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
if length(groups) == 2
    valid_model_types = [valid_model_types "logistic_regression_binary_glm"];
end
n_model_types = length(valid_model_types);

% generate classifiers of given model types
classifiers = cell(n_model_types,options.num_trials);
method_feature_selection = options.method_feature_selection;
num_features_selected = options.num_features_selected;
num_components_explained = options.num_components_explained;
compute_validation_accuracy = options.compute_validation_accuracy;
k_folds = options.k_folds;
verbose = options.verbose;
num_trials = options.num_trials;
parfor i=1:length(valid_model_types)
    warning('off','all')
    for j=1:num_trials
        classifiers{i,j} = run_SNAP_classification(T,...
            'method_feature_selection',method_feature_selection,...
            'num_components_explained',num_components_explained,...
            'compute_validation_accuracy',compute_validation_accuracy,...
            'model_type',valid_model_types(i),...
            'k_folds',k_folds,...
            'num_features_selected',num_features_selected,...
            'verbose',verbose);
    end
    warning('on','all')
end

% sort by validation accuracies
validation_accuracies_mean = zeros(1,n_model_types);
validation_accuracies_std = zeros(1,n_model_types);
for i=1:n_model_types
    validation_accuracies_mean(i) = mean(cellfun(@(x) x.ValidationAccuracy,classifiers(i,:),'uni',1));
    validation_accuracies_std(i) = std(cellfun(@(x) x.ValidationAccuracy,classifiers(i,:),'uni',1));
end
[validation_accuracies_mean_sort,validation_accuracies_sort_idx]=sort(validation_accuracies_mean,'descend');
validation_accuracies_std_sort = validation_accuracies_std(validation_accuracies_sort_idx);
classifiers = classifiers(validation_accuracies_sort_idx);
for i=1:n_model_types
    fprintf('%34s: %6.3f +/- %3.3f\n',classifiers{i}.ModelType,validation_accuracies_mean_sort(i)*100,validation_accuracies_std_sort(i)*100);
end
classifier_type = string(cellfun(@(x) x.ModelType, classifiers,'uni',0));
validation_accuracies_mean_sort = validation_accuracies_mean_sort';
validation_accuracies_std_sort = validation_accuracies_std_sort';
summary = table(classifier_type,validation_accuracies_mean_sort,validation_accuracies_std_sort);
end