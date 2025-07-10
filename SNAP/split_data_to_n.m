function [split_data, n] = split_data_to_n(data,n,options)
    arguments
        data
        n double
        options.shuffle logical = true
    end
    class_data = class(data);
    if isempty(data)
        split_data = cell(1,n);
        return
    end
    switch class_data
        case 'double'
            n_samples = numel(data);
            if options.shuffle
                data = data(randperm(n_samples));
            end
            if n_samples < n
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i);
                end
            else
                % divide into groups for parallel processing
                n_samp_per_process = get_n_samples_per_split(n_samples,n);
                split_data = cell(1,n);
                idx = 1;
                for i=1:n
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1);
                    idx = idx + n_samp_per_process(i);
                end
            end
        case 'cell'
            n_samples = size(data,2);
            if options.shuffle
                data = data(randperm(n_samples));
            end
            if n_samples < n
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i);
                end
            else
                % divide into groups for parallel processing
                n_samp_per_process = get_n_samples_per_split(n_samples,n);
                split_data = cell(1,n);
                idx = 1;
                for i=1:n
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1);
                    idx = idx + n_samp_per_process(i);
                end
            end
        case 'table'
            n_samples = size(data,1);
            if options.shuffle
                data = data(randperm(n_samples),:);
            end
            if n_samples < n
                % if the number of samples is less than parallel processing
                split_data = cell(1,n_samples);
                for i=1:n_samples
                    split_data{i} = data(i,:);
                end
            else
                % else, divide into groups for parallel processing
                n_samp_per_process = get_n_samples_per_split(n_samples,n);
                split_data = cell(1,n);
                idx = 1;
                for i=1:n
                    split_data{i} = data(idx:idx+n_samp_per_process(i)-1,:);
                    idx = idx + n_samp_per_process(i);
                end
            end
        otherwise
            ME = MException('SNAP:invalid_data_class_for_split', ...
                sprintf("Parameter 'data' has invalid class: %s",class_data));
            throw(ME)
    end
end

function n_samp_per_split = get_n_samples_per_split(n_samples,n_splits)
    n_samp_per_process_n_min1 = floor(n_samples/n_splits);
    n_samp_per_process_n = n_samples - n_samp_per_process_n_min1*(n_splits);
    n_samp_per_split = [repmat(n_samp_per_process_n_min1, 1, n_splits)];
    n_samp_per_split(1:n_samp_per_process_n) = n_samp_per_split(1:n_samp_per_process_n) + 1;
end