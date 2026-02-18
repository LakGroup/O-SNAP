% -------------------------------------------------------------------------
% generate_OSNAP_samples.m
% -------------------------------------------------------------------------
% Reads a data table that stores an Nx2 array with each row corresponding 
% to each localization and the columns are the x, y localizations in PIXELS
%
% Example on how to use it:
%   generate_OSNAP_samples("D:\Analysis\ExperimentA\sample_001.csv",117)
% -------------------------------------------------------------------------
% Input:
%   file_path: Full path to the sample file location
%   pixel_size: Conversion of pixel to nm
% Output:
%   data: A structure with the converted x,y localizations in nm, ready for
%         further processing in the OSNAP pipeline
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

function data = generate_OSNAP_samples(file_path,pixel_size)
arguments
    file_path string
    pixel_size double
end
    %% Load data
    locs_loaded = readtable(file_path,'ReadRowNames',0,'VariableNamingRule','preserve');
    %% Flip to match conventional images
    % y = max(y) - y; 
    % if pixel_size == 117 % pixel size in nm
    %     y = 684 - y;
    % end
    %% Clean data
    [~,data.name,~] = fileparts(file_path);
    data_to_unique(:,1) = locs_loaded{:,1};
    data_to_unique(:,2) = locs_loaded{:,2};
    data_to_unique = unique(data_to_unique,'rows');
    x = data_to_unique(:,1)*pixel_size;
    y = data_to_unique(:,2)*pixel_size;
    %% Assign localizations
    data.x = x;
    data.y = y;
    data.pixel_size_nm = pixel_size;
end

