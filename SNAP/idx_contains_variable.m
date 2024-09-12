function valid_idx = idx_contains_variable(data_info_table,vars_to_load)
    filepaths = data_info_table{:,'filepath'};
    valid_idx = ones(size(data_info_table,1),1);
    for i=1:length(filepaths)
        if ~has_variables(data_info_table{i,'filepath'},vars_to_load)
            valid_idx(i) = 0;
        end
    end
    if ~all(valid_idx)
        disp('Problem files:')
        f = data_info_table{~valid_idx,'filepath'};
        for i=1:length(f) 
            has_variables(f(i),vars_to_load,"verbose",1);
        end
    end
    valid_idx = find(valid_idx);
end
