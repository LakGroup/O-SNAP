% -------------------------------------------------------------------------
% save_OSNAP_sample.m
% -------------------------------------------------------------------------
% Function to facilitate saving sample data to a corresponding analysis
% file with error handling to facilitate batch processing
%
% Example on how to use it:
%   save_OSNAP_sample("D:\Analysis\ExperimentA\Rep1\sample_KO_1.mat",...
%                     sample_data);
% -------------------------------------------------------------------------
% Input:
%   save_path: Name of file location to save to
%   data: Struct arry containing information to save to MAT file for a
%         given sample (nucleus)
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
function save_OSNAP_sample(save_path, data)
    try
        save(save_path,...
            "-struct",...
            "data")
    catch ME
        disp(getReport(ME));
        disp("Error saving: " + save_path)
    end
end


