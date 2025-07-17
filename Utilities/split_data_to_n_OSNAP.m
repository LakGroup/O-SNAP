% -------------------------------------------------------------------------
% split_data_to_n_OSNAP.m
% -------------------------------------------------------------------------
% Splits the contents of data, which can take multiple forms, and divides
% it into a cell array with n_out cells. Designed to partition data for
% parallel computing.
%
% Example on how to use it:
%   [split_data, n] = split_data_to_n_OSNAP(data,12)
% -------------------------------------------------------------------------
% Input:
%   data: The data structure whose contents are to be partitioned into
%         n_out divisions
%   n_in: The number of splits to divide the data into
% Output:
%   split_data: The [1xn_out] cell array where each cell contains a 
%               subdivision of the original data. 
%   n_out: The final number of splits after partitioning the data. If the
%          number of samples in data is less than n_in, then n_out is
%          simply the number of samples present in data.
% Options:
%   shuffle: Flag for whether the contents of data should be shuffled prior
%            to partitioning
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
function [split_data, n_out] = split_data_to_n_OSNAP(data,n_in,options)
    arguments
        data
        n_in double
        options.shuffle logical = true
    end
    n_out = n_in;
    class_data = class(data);
    if isempty(data)
        split_data = cell(1,n_in);
        return
    end
    switch class_data
        case 'double'
            n_samples = numel(data);
            if options.shuffle
                data = data(randperm(n_samples));
            end
            if n_samples < n_in
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i);
                end
                n_out = n_samples;
            else
                % divide into groups for parallel processing
                n_samp_per_process = get_OSNAP_n_samples_per_split(n_samples,n_in);
                split_data = cell(1,n_in);
                idx = 1;
                for i=1:n_in
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1);
                    idx = idx + n_samp_per_process(i);
                end
            end
        case 'cell'
            n_samples = size(data,2);
            if options.shuffle
                data = data(randperm(n_samples));
            end
            if n_samples < n_in
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i);
                end
                n_out = n_samples;
            else
                % divide into groups for parallel processing
                n_samp_per_process = get_OSNAP_n_samples_per_split(n_samples,n_in);
                split_data = cell(1,n_in);
                idx = 1;
                for i=1:n_in
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1);
                    idx = idx + n_samp_per_process(i);
                end
            end
        case 'table'
            n_samples = size(data,1);
            if options.shuffle
                data = data(randperm(n_samples),:);
            end
            if n_samples < n_in
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i,:);
                end
                n_out = n_samples;
            else
                % else, divide into groups for parallel processing
                n_samp_per_process = get_OSNAP_n_samples_per_split(n_samples,n_in);
                split_data = cell(1,n_in);
                idx = 1;
                for i=1:n_in
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1,:);
                    idx = idx + n_samp_per_process(i);
                end
            end
        case 'string'
            n_samples = numel(data);
            if options.shuffle
                data = data(randperm(n_samples),:);
            end
            if n_samples < n_in
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i,:);
                end
                n_out = n_samples;
            else
                % else, divide into groups for parallel processing
                n_samp_per_process = get_OSNAP_n_samples_per_split(n_samples,n_in);
                split_data = cell(1,n_in);
                idx = 1;
                for i=1:n_in
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1);
                    idx = idx + n_samp_per_process(i);
                end
            end
        otherwise
            ME = MException('SNAP:invalid_data_class_for_split', ...
                sprintf("Parameter 'data' has invalid class: %s",class_data));
            throw(ME)
    end
end
