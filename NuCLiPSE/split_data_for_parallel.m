function [split_data, n_processes] = split_data_for_parallel(data, n_processes, shuffle)
    arguments
        data cell
        n_processes double = 12
        shuffle logical = true
    end
    n_samples = size(data,2);
    data_shuffled = data;
    if shuffle
        data_shuffled = data_shuffled(randperm(n_samples));
    end
    % if the number of samples is less than parallel processing
    if n_samples < n_processes
        disp('Number of samples is less than number of processes; not very efficient for parallel...')
        n_processes = n_samples;
        split_data = cell(1,n_samples);
        for i=1:n_samples
            split_data{i} = data(i);
        end
    end
    % else, divide into groups for parallel processing
    n_samp_per_process = get_n_samples_per_process(n_samples,n_processes);
    split_data = cell(1,n_processes);
    idx = 1;
    for i=1:n_processes
        split_data{i} = data_shuffled(idx:idx+n_samp_per_process(i)-1);
        idx = idx + n_samp_per_process(i);
    end
end
