function n_samp_per_split = get_OSNAP_n_samples_per_split(n_samples,n_splits)
    n_samp_per_process_n_min1 = floor(n_samples/n_splits);
    n_samp_per_process_n = n_samples - n_samp_per_process_n_min1*(n_splits);
    n_samp_per_split = [repmat(n_samp_per_process_n_min1, 1, n_splits)];
    n_samp_per_split(1:n_samp_per_process_n) = n_samp_per_split(1:n_samp_per_process_n) + 1;
end