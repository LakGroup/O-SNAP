% -------------------------------------------------------------------------
% find_contains_variable_OSNAP.m
% -------------------------------------------------------------------------
% Returns the indices for samples whose individual MAT analysis files
% (stored in the directories "OSNAP_nucleus_info") contain all of the
% variables specified
%
% Example on how to use it:
%   valid_idx = find_contains_variable_OSNAP(SNAP_nucleus_file_list,...
%                                            {'major_axis',...
%                                             'voronoi_cluster_radius'})
% -------------------------------------------------------------------------
% Input:
%   SNAP_nucleus_file_list: Table with details on the samples of interest
%                           to extract feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
%   vars_to_load: Cell array listing the variable names that the MAT file 
%                 for the sample must contain in order for it to be 
%                 considered valid to load from
% Output:
%   valid_idx: Indices for the samples/rows in SNAP_nucleus_file_list that
%              correspond to sample MAT analysis files that do contain all
%              the variables specified in vars_to_load
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
function valid_idx = find_contains_variable_OSNAP(SNAP_nucleus_file_list,vars_to_load)
    filepaths = SNAP_nucleus_file_list{:,'filepath'};
    valid_idx = ones(size(SNAP_nucleus_file_list,1),1);
    for i=1:length(filepaths)
        if ~has_variables_OSNAP(SNAP_nucleus_file_list{i,'filepath'},vars_to_load)
            valid_idx(i) = 0;
        end
    end
    if ~all(valid_idx)
        disp('Problem files:')
        f = SNAP_nucleus_file_list{~valid_idx,'filepath'};
        for i=1:length(f) 
            has_variables_OSNAP(f(i),vars_to_load,"verbose",1);
        end
    end
    valid_idx = find(valid_idx);
end


