% -------------------------------------------------------------------------
% generate_OSNAP_samples_batch.m
% -------------------------------------------------------------------------
% Organizes running generate_OSNAP_samples.m on multiple files in parallel
%
% Example on how to use it:
%   generate_OSNAP_samples_batch(OSNAP_sample_file_list,"D:\Analysis\ExperimentA",117)
% -------------------------------------------------------------------------
% Input:
%    OSNAP_sample_file_list: Table with details on the samples of interest
%                           to extract feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
%   work_dir: O-SNAP analysis directory. Should contain one or more folders
%             with the MAT O_SNAP feature files for each sample (nucleus),
%             divided by replicates with each folder prefixed by "Rep"
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
function runtime = generate_OSNAP_samples_batch(sample_info, work_dir, pixel_size)
n_processes = maxNumCompThreads;
% Split data
[split_data, n_processes] = split_data_to_n_OSNAP(sample_info, n_processes);
% Save scaled localizations
disp('Saving localizations...');
tic
parfor p=1:n_processes
    data_p = split_data{p};
    for s=1:size(data_p,1)
        % Get sample info
        sample_path = data_p{s,"Path"};
        sample_file = data_p{s,"Sample"};
        sample_replicate = data_p{s,"Replicate"};
        [~,sample_name,~] = fileparts(sample_file);
        % Directory management
        rep_dir = fullfile(work_dir,"Rep"+sample_replicate);
        save_dir = fullfile(rep_dir,"OSNAP_nucleus_data");
        if ~exist(rep_dir,'dir')
            mkdir(rep_dir)
        end
        if ~exist(save_dir,"dir")
            mkdir(save_dir)
        end
        save_path = fullfile(save_dir,(sample_name + ".mat"));
        % Create data object
        if ~exist(save_path,'file')
            sample_data = generate_OSNAP_samples(fullfile(sample_path,sample_file),pixel_size);
            % Save information
            save_OSNAP_sample(save_path, sample_data);
        end
        disp("   " + sample_name + ": " + string(datetime))
    end
end
runtime = toc;
disp('Done!')
end 

