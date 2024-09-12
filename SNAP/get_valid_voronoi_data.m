function data_info_table = get_valid_voronoi_data(work_dir,groups,reps,vars_to_load)
    if(exist(work_dir,"dir"))
        folder_info = dir(work_dir);
        folder_info = folder_info([folder_info(:).isdir]);
        dir_names = {folder_info(3:end).name};
        reps_name = cellfun(@(x) ['Rep' x], reps,'uni',0);
        valid_dirs_idx = cellfun(@(x) contains(x, reps_name),dir_names,'uni',1);
        valid_dirs = dir_names(valid_dirs_idx);
        %% get preliminary information on valid files to load based on input
        n_samples = 0;
        dir_names = [];
        sample_names = [];
        rep_names = [];
        group_names = [];
        for d=1:length(valid_dirs)
            valid_dir = string(fullfile(work_dir,valid_dirs(d),'voronoi_data'));
            folder_info = dir(valid_dir);
            folder_info = folder_info(~[folder_info(:).isdir]);
            file_names = {folder_info(:).name};
            valid_idx = all([endsWith(file_names,".mat"); contains(file_names, 'storm','IgnoreCase',true); contains(file_names, groups,'IgnoreCase',true)]);
            n_valid_files = sum(valid_idx);
            sample_names_d = string(extractBefore(file_names(valid_idx),".mat"))';
            dir_names = [dir_names; repmat(valid_dir,n_valid_files,1)];
            sample_names = [sample_names; sample_names_d];
            rep_names = [rep_names; repmat(string(reps{d}),n_valid_files,1)];
            group_names = [group_names; arrayfun(@(x) find_group(x,groups), sample_names_d,'uni',1)];
            n_samples = n_samples + n_valid_files;
        end
        filepaths = arrayfun(@(x,y) fullfile(x,y+".mat"),dir_names,sample_names,'uni',1);
        data_info_table = table(dir_names,sample_names,rep_names,group_names,filepaths);
        data_info_table.Properties.VariableNames={'dir','name','rep','group','filepath'};
        data_info_table = data_info_table(idx_contains_variable(data_info_table,vars_to_load),:);
        if ~size(data_info_table,1)
            ME = MException('SNAP:no_valid_files_found', ...
                'No files were found matching analysis parameters');
            throw(ME)
        end
    else
        disp("No valid directory: " + work_dir)
        ME = MException('SNAP:no_valid_directory_found', ...
                        sprintf('Directory not found: %s',work_dir));
        throw(ME)
    end
end

