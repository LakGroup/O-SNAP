function SNAP_nucleus_file_list = get_valid_OSNAP_samples(work_dir,groups,replicates,vars_to_load)
    save_name = "nucleus_SNAP_file_list.csv";
    if exist(fullfile(work_dir,save_name),'file')
        SNAP_nucleus_file_list = readtable(fullfile(work_dir,save_name),'NumHeaderLines',0,'TextType','string', 'VariableNamingRule', 'preserve');
        SNAP_nucleus_file_list.Properties.VariableTypes = repmat("string",1,5);
        reps_loaded = unique(SNAP_nucleus_file_list.rep);
        groups_loaded = unique(SNAP_nucleus_file_list.group);
        if (numel(replicates) > numel(reps_loaded)) || (numel(groups) > numel(groups_loaded))
            SNAP_nucleus_file_list = generate_valid_SNAP_nucleus_files(work_dir,groups,replicates,vars_to_load);
        end
        for i=1:numel(replicates)
            if ~ismember(replicates{i}, reps_loaded)
                SNAP_nucleus_file_list = generate_valid_SNAP_nucleus_files(work_dir,groups,replicates,vars_to_load);
            end
        end
        for i=1:numel(groups)
            if ~ismember(groups{i}, groups_loaded)
                SNAP_nucleus_file_list = generate_valid_SNAP_nucleus_files(work_dir,groups,replicates,vars_to_load);
            end
        end
    else
        SNAP_nucleus_file_list = generate_valid_SNAP_nucleus_files(work_dir,groups,replicates,vars_to_load);
    end
    % filter for groups and replicates
    SNAP_nucleus_file_list = SNAP_nucleus_file_list(ismember(SNAP_nucleus_file_list.rep,replicates),:);
    SNAP_nucleus_file_list = SNAP_nucleus_file_list(ismember(SNAP_nucleus_file_list.group,groups),:);
end

function SNAP_nucleus_file_list = generate_valid_SNAP_nucleus_files(work_dir,groups,reps,vars_to_load)
    if(exist(work_dir,"dir"))
        folder_info = dir(work_dir);
        folder_info = folder_info([folder_info(:).isdir]);
        dir_names = {folder_info(3:end).name};
        reps_name = cellfun(@(x) ['Rep' x], reps,'uni',0);
        valid_dirs_idx = cellfun(@(x) contains(x, reps_name),dir_names,'uni',1);
        valid_dirs = dir_names(valid_dirs_idx);
        save_name = "nucleus_SNAP_file_list.csv";
        %% get preliminary information on valid files to load based on input
        n_samples = 0;
        dir_names = [];
        sample_names = [];
        rep_names = [];
        group_names = [];
        for d=1:length(valid_dirs)
            valid_dir = string(fullfile(work_dir,valid_dirs(d),'SNAP_nucleus_data'));
            folder_info = dir(valid_dir);
            folder_info = folder_info(~[folder_info(:).isdir]);
            file_names = {folder_info(:).name};
            valid_idx = all([endsWith(file_names,".mat"); contains(file_names, groups,'IgnoreCase',true)]);
            n_valid_files = sum(valid_idx);
            sample_names_d = string(extractBefore(file_names(valid_idx),".mat"))';
            dir_names = [dir_names; repmat(valid_dir,n_valid_files,1)];
            sample_names = [sample_names; sample_names_d];
            rep_names = [rep_names; repmat(string(reps{d}),n_valid_files,1)];
            group_names = [group_names; arrayfun(@(x) find_group(x,groups), sample_names_d,'uni',1)];
            n_samples = n_samples + n_valid_files;
        end
        filepaths = arrayfun(@(x,y) fullfile(x,y+".mat"),dir_names,sample_names,'uni',1);
        SNAP_nucleus_file_list = table(dir_names,sample_names,rep_names,group_names,filepaths);
        SNAP_nucleus_file_list.Properties.VariableNames={'dir','name','rep','group','filepath'};
        SNAP_nucleus_file_list = SNAP_nucleus_file_list(find_contains_variable_OSNAP(SNAP_nucleus_file_list,vars_to_load),:);
        if ~size(SNAP_nucleus_file_list,1)
            ME = MException('SNAP:no_valid_files_found', ...
                'No files were found matching analysis parameters');
            throw(ME)
        end
        writetable(SNAP_nucleus_file_list,fullfile(work_dir,save_name),'WriteVariableNames',true);
    else
        disp("No valid directory: " + work_dir)
        ME = MException('SNAP:no_valid_directory_found', ...
                        sprintf('Directory not found: %s',work_dir));
        throw(ME)
    end
end

function group = find_group(sample_name,groups)
    group_idx = cell2mat(cellfun(@(x) contains(sample_name, x), groups,'uni',0));
    if sum(group_idx) == 1
        group = string(groups{group_idx});
    else
        group_idx = cell2mat(cellfun(@(x) contains(sample_name, ['_' x '-']), groups,'uni',0));
        if sum(group_idx) == 1
            group = string(groups{group_idx});
        else
            group_idx = cell2mat(cellfun(@(x) contains(sample_name, ['-' x '-']), groups,'uni',0));
            if sum(group_idx) == 1
                group = string(groups{group_idx});
            else
                group_idx = cell2mat(cellfun(@(x) contains(sample_name, ['-' x '_']), groups,'uni',0));
                if sum(group_idx) == 1
                    group = string(groups{group_idx});
                else
                    group_idx = cell2mat(cellfun(@(x) contains(sample_name, x), groups,'uni',0));
                    if any(group_idx)
                        group = string(groups{group_idx});
                    else
                        error('Error: At least one group specified could not be found in the data set')
                    end
                end
            end
        end
    end
end
