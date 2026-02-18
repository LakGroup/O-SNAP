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
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
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


