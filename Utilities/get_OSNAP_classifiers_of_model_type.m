% -------------------------------------------------------------------------
% get_OSNAP_classifiers_of_model_type.m
% -------------------------------------------------------------------------
% Returns classifiers that match the specified model type
%
% Example on how to use it:
%   classifiers_same_model = get_OSNAP_classifiers_of_model_type(model_type,...
%                                                                classifiers)
% -------------------------------------------------------------------------
% Input:
%   model_type: Desired model type to filter
%   classifiers: A cell array where each cell stores a single instance of a
%                model trained on the O-SNAP feature data
% Output:
%   classifiers_same_model: A cell array of models derived from classifiers
%                           that contain only those that match the
%                           specified model type
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
function classifiers_same_model = get_OSNAP_classifiers_of_model_type(model_type,classifiers)
arguments
    model_type string
    classifiers cell
end
n_batch = numel(classifiers);
classifiers_same_model = cell(1,n_batch);
%% For each batch, get all clasifiers that match the specified model type
for i=1:n_batch
    model_type_list = string(cellfun(@(x) x.ModelType, classifiers{i},'uni',0));
    classifiers_same_model{i} = classifiers{i}(strcmp(model_type,model_type_list));
end
classifiers_same_model = [classifiers_same_model{:}];
end
