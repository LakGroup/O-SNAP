% -------------------------------------------------------------------------
% evaluate_OSNAP_classifiers.m
% -------------------------------------------------------------------------
% Evaluates and summarizes the performance of a set of classifiers trained
% on the same set of O-SNAP feature data. Produces ROC curves, confusion
% matrices, and a summary table of the overall classification performances
% of the models.
%
% Example on how to use it:
%   classification_summary = evaluate_OSNAP_classifiers(classifiers)
% -------------------------------------------------------------------------
% Input:
%   classifiers: A cell array where each cell stores a single instance of a
%                model trained on the O-SNAP feature data
% Output:
%   classification_summary: A table with the overall classification
%                           accuracies for all model types present in the
%                           classifier variable
%         Columns if test data used:
%                   - mean_batch_acc: Average accuracy across batches
%                   - std_batch_acc: Std dev across batches
%                   - std_model_acc: Std dev within a model type
%                   - n_batch: Number of batches (default: 5)
%         Columns if NO test data used:
%                   - mean_model_acc: Average accuracy within model type
%                   - std_model_acc: Std dev within a model type
% Options:
%   save_path: Name of file location to save to, excluding extension
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
%%
function classification_summary = evaluate_OSNAP_classifiers(classifiers,options)
arguments
    classifiers cell = {}
    options.save_path string = ""
end
n_batch = numel(classifiers);
n_model_replicates = size(classifiers{1},1);

%% Overall classification accuracy across all model types
model_type = cellfun(@(x) x.ModelType, [classifiers{:}], 'uni', 1);
model_type = model_type(1,:)';
% If the classifier was used against test data, report performance across
% the batches
if isfield(classifiers{1}{1},"TestResponse")
    groups = sort(unique(classifiers{1}{1}.TestResponse));
    mean_model_acc = mean(cellfun(@(x) x.TestAccuracy, [classifiers{:}], 'uni', 1),1)';
    std_model_acc = std(cellfun(@(x) x.TestAccuracy, [classifiers{:}], 'uni', 1),1)';
    classification_summary = table(model_type,mean_model_acc,std_model_acc);
    classification_summary = groupsummary(classification_summary,"model_type",["mean","std"]);
    classification_summary.std_std_model_acc = [];
    classification_summary = renamevars(classification_summary,...
        {'mean_mean_model_acc','std_mean_model_acc','mean_std_model_acc','GroupCount'},...
        {'mean_batch_acc','std_batch_acc','std_model_acc','n_batch'});
    classification_summary = sortrows(classification_summary(:,[1 3 4 5 2]),...
        {'mean_batch_acc','std_batch_acc'},{'descend','ascend'});
% If the classifier was only tested using validation, only show the
% accuracy statistics based on model type (not batches)
else
    groups = sort(unique(classifiers{1}{1}.ValidationResponse));
    mean_model_acc = mean(cellfun(@(x) x.ValidationAccuracy, [classifiers{:}], 'uni', 1))';
    std_model_acc = std(cellfun(@(x) x.ValidationAccuracy, [classifiers{:}], 'uni', 1),1)';
    classification_summary = sortrows(table(model_type,mean_model_acc,std_model_acc),...
        {'mean_model_acc','std_model_acc'},{'descend','ascend'});
end

% Find the best classifier model type
n_groups = numel(groups);
model_type_best = classification_summary{1,"model_type"};
classifiers_same_model = get_OSNAP_classifiers_of_model_type(model_type_best,classifiers);

%% Calculate confusion matrix
if isfield(classifiers{1}{1},"TestResponse")
    confusion_counts_cell = cellfun(@(x) ...
        confusionmat(string(x.TestResponse),...
                     string(x.TestPredictions),...
                     'order',groups), ...
        classifiers_same_model,'uni',0);
else
    confusion_counts_cell = cellfun(@(x) ...
        confusionmat(string(x.ValidationResponse),...
                     string(x.ValidationPredictions),...
                     'order',groups), ...
        classifiers_same_model,'uni',0);
end
% Stores counts of confusion matrices in a 3-D matrix (depth corresponds to
% batch)
confusion_counts = zeros(n_groups,n_groups,n_batch);
for j=1:n_batch
    confusion_counts(:,:,j) = mean(reshape(cell2mat(confusion_counts_cell(:,j)'),n_groups,n_groups,[]),3);
end

%% Calculate ROC data
ROC_data = cell(n_model_replicates,n_batch);
% x=classifiers_same_model{1,1};rocObj = rocmetrics(x.TestResponse,x.TestScores,x.ClassificationModel.ClassNames);
for b=1:n_batch
    for m=1:n_model_replicates
        ROC_data{m,b} = classifiers_same_model{m,b}.PerformanceCurve;
    end
end

%% Plot confusion matrix and ROC
if options.save_path ~= ""
    % Confusion matrix
    figure('visible','off');
    plot_OSNAP_confusion(gca,groups,confusion_counts,"model_type_name",model_type_best);
    savefig(gcf,options.save_path+"_confusion.fig");
    exportgraphics(gcf,options.save_path+"_confusion.png");
    % ROC
    figure('visible','off');
    plot_OSNAP_ROC(gca,groups,ROC_data);
    savefig(gcf,options.save_path+"_ROC.fig");
    exportgraphics(gcf,options.save_path+"_ROC.png");
else
    % Confusion matrix
    figure();
    plot_OSNAP_confusion(gca,groups,confusion_counts,"model_type_name",model_type_best);
    % ROC
    figure();
    plot_OSNAP_ROC(gca,groups,ROC_data);
end

end


