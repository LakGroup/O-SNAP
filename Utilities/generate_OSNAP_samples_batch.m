function runtime = generate_OSNAP_samples_batch(sample_info, work_dir, pixel_size)
n_processes = maxNumCompThreads;
% split data
[split_data, n_processes] = split_data_to_n_OSNAP(sample_info, n_processes);
% save scaled localizations
disp('Saving localizations...');
tic
parfor p=1:n_processes
    data_p = split_data{p};
    for s=1:size(data_p,1)
        % get sample info
        sample_path = data_p{s,"Path"};
        sample_file = data_p{s,"Sample"};
        sample_replicate = data_p{s,"Replicate"};
        [~,sample_name,~] = fileparts(sample_file);
        % directory management
        rep_dir = fullfile(work_dir,"Rep"+sample_replicate);
        save_dir = fullfile(rep_dir,"OSNAP_nucleus_data");
        if ~exist(rep_dir,'dir')
            mkdir(rep_dir)
        end
        if ~exist(save_dir,"dir")
            mkdir(save_dir)
        end
        save_path = fullfile(save_dir,(sample_name + ".mat"));
        % create data object
        sample_data = generate_OSNAP_samples(fullfile(sample_path,sample_file),pixel_size);
        % save
        save_OSNAP_sample(save_path, sample_data);
        disp("   " + sample_name + ": " + string(datetime))
    end
end
runtime = toc;
disp('Done!')
end 

