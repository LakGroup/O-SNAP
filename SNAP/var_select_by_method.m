function t = var_select_by_method(var_select_result)
    n_batch = size(var_select_result,2);
    var_select_result_nonempty=vertcat(var_select_result{:})';
    var_select_result_nonempty=var_select_result_nonempty(~cellfun('isempty',var_select_result_nonempty(:,1)),:);
    if isempty(var_select_result_nonempty)
        t = [];
    else
        var_sel_method = string(cellfun(@(x) x.Method, var_select_result_nonempty(:,1), 'uni',0));
        for i=1:size(var_select_result_nonempty,1)
            var_selected_method_by_batch = cellfun(@(x) x.VarSel, var_select_result_nonempty(i,:), 'uni',0);
            select_freq_by_batch = sum(vertcat(var_selected_method_by_batch{:}),1)'/n_batch;
            if ~isempty(select_freq_by_batch)
                t_data.(var_sel_method(i)) = select_freq_by_batch;
            end
        end
        t = struct2table(t_data);
    end
end