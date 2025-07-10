% -------------------------------------------------------------------------
% get_OSNAP_n_samples_per_split.m
% -------------------------------------------------------------------------
% Outputs the number of samples present in each split
%
% Example on how to use it:
%   n_samp_per_split = get_OSNAP_n_samples_per_split(n_samples,n_splits)
% -------------------------------------------------------------------------
% Input:
%   n_samples: Total number of samples
%   n_split: Number of splits to divide samples into
% Output:
%   n_samp_per_split: [1xn_split] array of the number of samples per
%                     split
% -------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function n_samp_per_split = get_OSNAP_n_samples_per_split(n_samples,n_splits)
    n_samp_per_process_n_min1 = floor(n_samples/n_splits);
    n_samp_per_process_n = n_samples - n_samp_per_process_n_min1*(n_splits);
    n_samp_per_split = [repmat(n_samp_per_process_n_min1, 1, n_splits)];
    n_samp_per_split(1:n_samp_per_process_n) = n_samp_per_split(1:n_samp_per_process_n) + 1;
end
