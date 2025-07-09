%% returns classifiers that match the model type
function classifiers_same_model = get_OSNAP_classifiers_of_model_type(model_type,classifiers)
arguments
    model_type string
    classifiers cell
end
n_batch = numel(classifiers);
classifiers_same_model = cell(1,n_batch);
for i=1:n_batch
    model_type_list = string(cellfun(@(x) x.ModelType, classifiers{i},'uni',0));
    classifiers_same_model{i} = classifiers{i}(strcmp(model_type,model_type_list));
end
classifiers_same_model = [classifiers_same_model{:}];
end