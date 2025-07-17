% -------------------------------------------------------------------------
% get_all_valid_OSNAP_samples.m
% -------------------------------------------------------------------------
% Creates a table with file path information for all samples that fall
% under the specified phenotype and replicate combinations.
%
% Example on how to use it:
%   OSNAP_sample_file_list = get_all_valid_OSNAP_samples("D:\Analysis\ExperimentA",...
%                                 {'major_axis','voronoi_cluster_radius'})
% -------------------------------------------------------------------------
% Input:
%   work_dir: O-SNAP analysis directory. Should contain one or more folders
%             with the MAT O-SNAP feature files for each sample (nucleus),
%             divided by replicates with each folder prefixed by "Rep"
%   vars_to_load: Cell array listing the variable names that the MAT file 
%                 for the sample must contain in order for it to be 
%                 considered valid to load from
% Output:
%   OSNAP_sample_file_list: Table with details on the samples of interest
%                           to extract feature information from, including 
%                           the file path, replicate, and sample
%                           identifier. Notably, does not include
%                           phenotype, unlike get_valid_OSNAP_samples
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
function OSNAP_sample_file_list = get_all_valid_OSNAP_samples(work_dir,vars_to_load)
    if(exist(work_dir,"dir"))
        folder_info = dir(work_dir);
        folder_info = folder_info([folder_info(:).isdir]);
        dir_names = {folder_info(3:end).name};
        valid_dirs_idx = cellfun(@(x) startsWith(x, "Rep"),dir_names,'uni',1);
        valid_dirs = dir_names(valid_dirs_idx);
        save_name = "OSNAP_sample_file_list.csv";
        %% get preliminary information on valid files to load based on input
        n_valid_files = zeros(1,numel(valid_dirs));
        % get number of valid files to instantiate variables
        for i=1:numel(valid_dirs)
            valid_dir = string(fullfile(work_dir,valid_dirs(i),'OSNAP_nucleus_data'));
            folder_info = dir(valid_dir);
            folder_info = folder_info(~[folder_info(:).isdir]);
            file_names_i = {folder_info(:).name};
            valid_idx = endsWith(file_names_i,".mat");
            n_valid_files(i) = sum(valid_idx);
        end
        % loop back through to add values
        file_names = strings(sum(n_valid_files),1);
        rep_names = strings(sum(n_valid_files),1);
        n=1;
        for i=1:numel(valid_dirs)
            valid_dir = string(fullfile(work_dir,valid_dirs(i),'OSNAP_nucleus_data'));
            folder_info = dir(valid_dir);
            folder_info = folder_info(~[folder_info(:).isdir]);
            file_names_i = {folder_info(:).name};
            valid_idx = endsWith(file_names_i,".mat");
            file_names(n:n+n_valid_files(i)-1) = fullfile(valid_dir,string(file_names_i(valid_idx)));
            rep_names(n:n+n_valid_files(i)-1) = repmat(extractAfter(valid_dirs(i),"Rep"),n_valid_files(i),1);
            n = n+n_valid_files(i);
        end
        [~,sample_names,~] = fileparts(file_names);
        OSNAP_sample_file_list = table(sample_names,rep_names,file_names);
        OSNAP_sample_file_list.Properties.VariableNames={'name','replicate','filepath'};
        OSNAP_sample_file_list = OSNAP_sample_file_list(find_contains_variable_OSNAP(OSNAP_sample_file_list,vars_to_load),:);
        if ~size(OSNAP_sample_file_list,1)
            ME = MException('OSNAP:no_valid_files_found', ...
                'No files were found matching analysis parameters');
            throw(ME)
        end
        writetable(OSNAP_sample_file_list,fullfile(work_dir,save_name),'WriteVariableNames',true);
    else
        disp("No valid directory: " + work_dir)
        ME = MException('OSNAP:no_valid_directory_found', ...
                        sprintf('Directory not found: %s',work_dir));
        throw(ME)
    end
end