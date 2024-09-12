function n_samp_per_process = get_n_samples_per_process(n_samples,n_processes)
    n_samp_per_process_n_min1 = floor(n_samples/n_processes);
    n_samp_per_process_n = n_samples - n_samp_per_process_n_min1*(n_processes);
    n_samp_per_process = [repmat(n_samp_per_process_n_min1, 1, n_processes)];
    n_samp_per_process(1:n_samp_per_process_n) = n_samp_per_process(1:n_samp_per_process_n) + 1;
end

