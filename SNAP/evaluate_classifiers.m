function classification_summary = evaluate_classifiers(classifiers,options)
arguments
    classifiers cell = {}
    options.save_path string = ""
end
n_batch = numel(classifiers);
n_model_replicates = size(classifiers{1},1);

%% overall classification accuracy across all model types
model_type = cellfun(@(x) x.ModelType, [classifiers{:}], 'uni', 1);
model_type = model_type(1,:)';
if isfield(classifiers{1}{1},"TestResponse")
    groups = sort(unique(classifiers{1}{1}.TestResponse));
    mean_model_acc = mean(cellfun(@(x) x.TestAccuracy, [classifiers{:}], 'uni', 1),1)';
    std_model_acc = std(cellfun(@(x) x.TestAccuracy, [classifiers{:}], 'uni', 1),1)';
    classification_summary = table(model_type,mean_model_acc,std_model_acc);
    classification_summary = groupsummary(classification_summary,"model_type",["mean","std"]);
    classification_summary.std_std_model_acc = [];
    classification_summary = renamevars(classification_summary,{'mean_mean_model_acc','std_mean_model_acc','mean_std_model_acc','GroupCount'},{'mean_batch_acc','std_batch_acc','std_model_acc','n_batch'});
    classification_summary = sortrows(classification_summary(:,[1 3 4 5 2]),{'mean_batch_acc','std_batch_acc'},{'descend','ascend'});
else
    groups = sort(unique(classifiers{1}{1}.ValidationResponse));
    mean_model_acc = mean(cellfun(@(x) x.ValidationAccuracy, [classifiers{:}], 'uni', 1))';
    std_model_acc = std(cellfun(@(x) x.ValidationAccuracy, [classifiers{:}], 'uni', 1),1)';
    classification_summary = sortrows(table(model_type,mean_model_acc,std_model_acc),{'mean_model_acc','std_model_acc'},{'descend','ascend'});
end

% find the best classifier model type
n_groups = numel(groups);
model_type_best = classification_summary{1,"model_type"};
classifiers_same_model = get_classifers_by_model_type(model_type_best,classifiers);

%% confusion matrix
if isfield(classifiers{1}{1},"TestResponse")
    confusion_counts_cell = cellfun(@(x) confusionmat(string(x.TestResponse),string(x.TestPredictions),'order',groups), ...
        classifiers_same_model,'uni',0);
else
    confusion_counts_cell = cellfun(@(x) confusionmat(string(x.ValidationResponse),string(x.ValidationPredictions),'order',groups), ...
            classifiers_same_model,'uni',0);
end
confusion_counts = zeros(n_groups,n_groups,n_batch);
for j=1:n_batch
    confusion_counts(:,:,j) = mean(reshape(cell2mat(confusion_counts_cell(:,j)'),n_groups,n_groups,[]),3);
end

%% ROC data
ROC_data = cell(n_model_replicates,n_batch);
% x=classifiers_same_model{1,1};rocObj = rocmetrics(x.TestResponse,x.TestScores,x.ClassificationModel.ClassNames);
for b=1:n_batch
    for m=1:n_model_replicates
        ROC_data{m,b} = classifiers_same_model{m,b}.PerformanceCurve;
    end
end

%% plot
if options.save_path ~= ""
    figure('visible','off');
    plot_SNAP_confusion(gca,groups,confusion_counts,"model_type_name",model_type_best);
    savefig(gcf,options.save_path+"_confusion.fig");
    exportgraphics(gcf,options.save_path+"_confusion.png");
    figure('visible','off');
    plot_SNAP_ROC(gca,groups,ROC_data);
    savefig(gcf,options.save_path+"_ROC.fig");
    exportgraphics(gcf,options.save_path+"_ROC.png");
else
    figure();
    plot_SNAP_confusion(gca,groups,confusion_counts,"model_type_name",model_type_best);
    figure();
    plot_SNAP_ROC(gca,groups,ROC_data);
end

end

