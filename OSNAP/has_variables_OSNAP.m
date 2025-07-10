% -------------------------------------------------------------------------
% has_variables_OSNAP.m
% -------------------------------------------------------------------------
% Checks whether a file contains given variables.
%
% Example on how to use it:
%   [result, vars_missing] = has_variables_OSNAP(...
%                                   "D:\Analysis\ExperimentA\Rep1\sample_KO_1.mat",...
%                                   {'major_axis',...
%                                    'voronoi_cluster_radius'})
% -------------------------------------------------------------------------
% Input:
%   file_path: File path to specific sample (nucleus) to analyze.
%   vars_to_check: Cell array listing the variable names that the MAT file 
%                  for the sample is being evaluated for
% Output:
%   result: Logical flag that is true when file contains all variables in
%           vars_to_check
%   vars_missing: A cell array listing missing variables, if there are any
% Options:
%   verbose: Flag for output to print
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
function [result, vars_missing] = has_variables_OSNAP(file_path, vars_to_check, options)
    arguments
        file_path string
        vars_to_check cell
        options.verbose logical = false
    end
    vars_in_file = who("-file",file_path);
    vars_missing = setdiff(vars_to_check,vars_in_file);
    if ~isempty(vars_missing)
        result = 0;
        if options.verbose
            fprintf('''%s'' not in: %s\n',string(join(vars_missing, ''', ''')),file_path)
        end
    else
        result = 1;
    end
end


