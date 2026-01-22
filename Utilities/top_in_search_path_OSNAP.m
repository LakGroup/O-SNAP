% -------------------------------------------------------------------------
% top_in_search_path_OSNAP.m
% -------------------------------------------------------------------------
% In the event that a function is overshadowed, this function sets the
% priority of the function definition in the file path containing 
% "path_sub_string" to the top by removing and re-adding it to the path;
% this is meant to effectively select the function definition if multiple
% versions exist in the path.
%
% The version actively being used can be checked using the which() function
%
% Example on how to use it:
%   top_in_search_path_OSNAP('pca','stats\stats\pca.m');
% -------------------------------------------------------------------------
% Input:
%   function_name: Function call that needs to have its priority changed
%   path_substring: A substring indicates the desired path of interest. The
%                   path to the file with the desired function definition
%                   should contain this string.
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
function top_in_search_path_OSNAP(function_name,path_substring)
arguments
    function_name string
    path_substring string = []
end
    filepaths = which(function_name,'-all');
    if contains(filepaths{1},path_substring)
        return
    end
    for i=2:numel(filepaths)
        if contains(filepaths{i},path_substring)
            dir_top = fileparts(filepaths{i});
            rmpath(dir_top);
            addpath(dir_top);
            break
        end
    end
end
