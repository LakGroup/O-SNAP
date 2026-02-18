% -------------------------------------------------------------------------
% generate_OSNAP_feature_set.m
% -------------------------------------------------------------------------
% Exports information on a feature universe from a struct to an XLSX file
% that is subsequently loaded during FSEA.
%
% Example on how to use it:
%   generate_OSNAP_feature_set(feature_IDs,feature_universe)
% -------------------------------------------------------------------------
% Input:
%   feature_IDs: String array of the feature identifiers (O-SNAP feature 
%                names)
%   feature_universes: A struct array with information on each feature 
%                      universe
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
function generate_OSNAP_feature_set(feature_IDs,feature_universe)
arguments
    feature_IDs cell
    feature_universe struct
end
    filepath = feature_universe.name+".xlsx";
    %% Clean feature universe if features are missing (i.e. from NaN filtering) from data
    feature_set_names = fieldnames(feature_universe.feature_set);
    n_feature_sets = length(feature_set_names);
    for s=1:n_feature_sets
        idx = cellfun(@(x) find(strcmpi(feature_IDs,x)),feature_universe.feature_set.(feature_set_names{s}),'uni',0);
        idx_valid = ~cellfun('isempty',idx);
        if all(~idx_valid)
            feature_universe.feature_set.(feature_set_names{s}) = rmfield(feature_universe.feature_set.(feature_set_names{s}));
        elseif any(~idx_valid)
            feature_universe.feature_set.(feature_set_names{s}) = feature_universe.feature_set.(feature_set_names{s})(idx_valid);
        end
    end
    %% Organize feature set information based on the feature set universe
    % Includes features belonging to feature set and direction of how a
    % feature increase aligns to a feature set description
    feature_set_names = fieldnames(feature_universe.feature_set);
    n_feature_sets = length(feature_set_names);
    feature_idxs = cell(n_feature_sets,1);
    for s=1:n_feature_sets
        idx = cellfun(@(x) find(strcmpi(feature_IDs,x)),feature_universe.feature_set.(feature_set_names{s}),'uni',0);
        idx_valid = ~cellfun('isempty',idx);
        idx = cell2mat(idx(idx_valid));
        direction = feature_universe.direction.(feature_set_names{s});
        direction = direction(idx_valid);
        feature_idxs{s} = (idx.*direction)';
    end
    max_n_features = max(cellfun(@(x) length(x), feature_idxs, 'uni',1));
    data_feature_set = table('Size',[n_feature_sets,2+max_n_features],...
        'VariableTypes',[{'string','string'}  repmat({'double'},1,max_n_features)]);
    for s=1:n_feature_sets
        data_feature_set{s,1} = feature_set_names(s);
        data_feature_set{s,2} = {NaN};
        data_feature_set{s,3:max_n_features+2} = [feature_idxs{s} NaN(1,max_n_features - length(feature_idxs{s}))];
    end
    %% Save results
    writetable(data_feature_set,filepath,"WriteVariableNames",false,"WriteMode","overwritesheet");
    update_FS(char(filepath))
end

