% -------------------------------------------------------------------------
% filter_OSNAP_feature_data.m
% -------------------------------------------------------------------------
% Filters an O-SNAP feature data table based on the specified phenotypes
% and replicates. Useful if certain groups or replicates need to be
% excluded for an analysis.
%
% Example on how to use it:
%   feature_data_filtered = filter_OSNAP_feature_data(feature_data,...
%                                      {'Control','KO'},...
%                                      {'20250101'})
% -------------------------------------------------------------------------
% Input:
%   feature_data: The feature data output table, where each row represents 
%      a sample (nucleus). The first three columns represent (1) 
%      Group/Phenotype, (2) replicate, and (3) Sample Identifier. Each 
%      subsequent column is an O-SNAP feature. 
%   groups: Cell array containing char array of the identifiers of the 
%           phenotypes/cell states
%           *ENSURE THAT EACH SAMPLE FILENAME CONTAINS EXACTLY ONE 
%            IDENTIFIER FROM groups*
%   replicates: Cell array containing char array of the identifiers of the 
%               replicates.
% Output:
%   feature_data_filtered: The filtered feature data based on requested
%                          phenotypes and replicates
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
function feature_data_filtered = filter_OSNAP_feature_data(feature_data,groups,replicates,options)
arguments
    feature_data table
    groups cell
    replicates cell
    options.remove_NaN logical = false
end
% filter for groups and replicates
fprintf("  Filtering feature data table for desired groups and replicates...\n")
feature_data_filtered = preprocess_OSNAP_feature_data(feature_data,...
    "groups",groups,"replicates",replicates,...
    "normalize",false,"remove_NaN",options.remove_NaN,...
    "keep_rep_sample_info",true);
% end execution if not enough groups after filtering
if numel(unique(feature_data_filtered.group)) < 2
    ME = MException('SNAP:insufficient_group_number', ...
                    sprintf('Only one group found: %s',groups{1}));
    throw(ME)
end
end
