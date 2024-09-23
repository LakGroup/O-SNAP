function [split_data, n_workers] = split_files_for_parallel(data_info_table, n_workers, shuffle)
    arguments
        data_info_table table
        n_workers double
        shuffle logical = true
    end
    n_samples = size(data_info_table,1);
    data_info_table_shuffled = data_info_table;
    if shuffle
        data_info_table_shuffled = data_info_table_shuffled(randperm(n_samples),:);
    end
    % if the number of samples is less than parallel processing
    if n_samples < n_workers
        disp('Number of samples is less than number of processes; not very efficient for parallel...')
        n_workers = n_samples;
        split_data = cell(1,n_samples);
        for i=1:n_samples
            split_data{i} = data_info_table_shuffled(i,:);
        end
    end
    % else, divide into groups for parallel processing
    n_samp_per_process = get_n_samples_per_process(n_samples,n_workers);
    split_data = cell(1,n_workers);
    idx = 1;
    for i=1:n_workers
        split_data{i} = data_info_table_shuffled(idx:idx+n_samp_per_process(i)-1,:);
        idx = idx + n_samp_per_process(i);
    end
end

