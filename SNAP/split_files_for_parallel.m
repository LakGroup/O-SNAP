function [split_data, n] = split_table_to_n(data,n,options)
    arguments
        data table
        n double
        options.shuffle logical = true
    end
    n_samples = size(data,1);
    if options.shuffle
        data = data(randperm(n_samples),:);
    end
    % if the number of samples is less than parallel processing
    if n_samples < n
        disp('Number of samples is less than number of processes; not very efficient for parallel...')
        n = n_samples;
        split_data = cell(1,n_samples);
        for i=1:n_samples
            split_data{i} = data(i,:);
        end
    end
    % else, divide into groups for parallel processing
    n_samp_per_process = get_n_samples_per_process(n_samples,n);
    split_data = cell(1,n);
    idx = 1;
    for i=1:n
        split_data{i} = data(idx:idx+n_samp_per_process(i)-1,:);
        idx = idx + n_samp_per_process(i);
    end
end

