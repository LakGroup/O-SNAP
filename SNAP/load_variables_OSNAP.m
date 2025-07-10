% -------------------------------------------------------------------------
% load_variables_OSNAP.m
% -------------------------------------------------------------------------
% Load specified variables from file
%
% Example on how to use it:
%   data = load_variables_OSNAP("D:\Analysis\ExperimentA\Rep1\sample_KO_1.mat",...
%                               {'major_axis',...
%                               'voronoi_cluster_radius'})
% -------------------------------------------------------------------------
% Input:
%   file_path: File path to specific sample (nucleus) to analyze.
%   vars_to_check: Cell array listing the variable names to load from MAT
%                  file
% Output:
%   data: Struct array containing loaded variables 
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
function data = load_variables_OSNAP(file_path, vars_to_load, options)
    arguments
        file_path string
        vars_to_load cell
        options.verbose logical = true
    end
    try
        [file_has_variables_OSNAP, vars_missing] = has_variables_OSNAP(file_path, vars_to_load, "verbose", options.verbose);
        if file_has_variables_OSNAP
            for i=1:length(vars_to_load)
                try
                    vars_value = load(file_path,vars_to_load{i});
                    data.(vars_to_load{i}) = vars_value.(vars_to_load{i});
                catch
                    ME = MException('OSNAP:failed_to_load_existing_variable', ...
                        'Failed to load ''%s'' from %s',vars_to_load{i},file_path);
                    throw(ME)
                end
            end
        else
            ME = MException('OSNAP:variable_not_in_file', ...
                        'File is missing ''%s''',string(join(vars_missing,''',''')));
            throw(ME)
        end
    catch ME
        throw(ME);
    end
end

