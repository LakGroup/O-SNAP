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
