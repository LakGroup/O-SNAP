% -------------------------------------------------------------------------
% get_valid_OSNAP_samples.m
% -------------------------------------------------------------------------
% Creates a table with file path information for all samples that fall
% under the specified phenotype and replicate combinations.
%
% Example on how to use it:
%   OSNAP_nucleus_file_list = get_valid_OSNAP_samples("D:\Analysis\ExperimentA",...
%                                 {'Control','KO'},...
%                                 {'20250101','20250108','20250201'},...
%                                 {'major_axis','voronoi_cluster_radius'})
% -------------------------------------------------------------------------
% Input:
%   work_dir: O-SNAP analysis directory. Should contain one or more folders
%             with the MAT O-SNAP feature files for each sample (nucleus),
%             divided by replicates with each folder prefixed by "Rep"
%   groups: Cell array containing char array of the identifiers of the 
%           phenotypes/cell states
%           *ENSURE THAT EACH SAMPLE FILENAME CONTAINS EXACTLY ONE 
%            IDENTIFIER FROM groups*
%   replicates: Cell array containing char array of the identifiers of the 
%               replicates.
%               *ENSURE THAT EACH REPLCIATE IS STORED IN A SEPARATE FOLDER 
%                PREFIXED BY "Rep"*
%   vars_to_load: Cell array listing the variable names that the MAT file 
%                 for the sample must contain in order for it to be 
%                 considered valid to load from
% Output:
%   OSNAP_sample_file_list: Table with details on the samples of interest
%                           to extract feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
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
function OSNAP_sample_file_list = get_valid_OSNAP_samples(work_dir,groups,replicates,vars_to_load)
    % save_name = "OSNAP_sample_file_list.csv";
    % if exist(fullfile(work_dir,save_name),'file')
    %     OSNAP_nucleus_file_list = readtable(fullfile(work_dir,save_name),'NumHeaderLines',0,'TextType','string', 'VariableNamingRule', 'preserve');
    %     OSNAP_nucleus_file_list.Properties.VariableTypes = repmat("string",1,5);
    %     reps_loaded = unique(OSNAP_nucleus_file_list.replicates);
    %     groups_loaded = unique(OSNAP_nucleus_file_list.group);
    %     if (numel(replicates) > numel(reps_loaded)) || (numel(groups) > numel(groups_loaded))
    %         OSNAP_nucleus_file_list = search_and_find_valid_OSNAP_samples(work_dir,groups,replicates,vars_to_load);
    %     end
    %     for i=1:numel(replicates)
    %         if ~ismember(replicates{i}, reps_loaded)
    %             OSNAP_nucleus_file_list = search_and_find_valid_OSNAP_samples(work_dir,groups,replicates,vars_to_load);
    %         end
    %     end
    %     for i=1:numel(groups)
    %         if ~ismember(groups{i}, groups_loaded)
    %             OSNAP_nucleus_file_list = search_and_find_valid_OSNAP_samples(work_dir,groups,replicates,vars_to_load);
    %         end
    %     end
    % else
    OSNAP_sample_file_list = get_all_valid_OSNAP_samples(work_dir,vars_to_load);
    OSNAP_sample_file_list.group = arrayfun(@(x) find_group(x,groups), OSNAP_sample_file_list{:,"name"},'uni',1);
    OSNAP_sample_file_list = OSNAP_sample_file_list(:,["name","group","replicate","filepath"]);
    % end
    % filter for groups and replicates
    OSNAP_sample_file_list = OSNAP_sample_file_list(ismember(OSNAP_sample_file_list.replicate,replicates),:);
    OSNAP_sample_file_list = OSNAP_sample_file_list(ismember(OSNAP_sample_file_list.group,groups),:);
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

