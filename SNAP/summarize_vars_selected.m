function vars_selected_summary = summarize_vars_selected(feature_data,vars_selected, var_select_result)
arguments
    feature_data table
    vars_selected cell
    var_select_result cell
end
[freq_per_batch,v]=groupcounts(convertCharsToStrings([vars_selected{:}])');
freq_per_batch = freq_per_batch/size(vars_selected,2);
t = var_select_by_method(var_select_result);
t.Properties.RowNames = preprocess_SNAP_table(feature_data,"numeric_only",true,"normalize",true).Properties.VariableNames;
t = outerjoin(t,table(freq_per_batch,'RowNames',v),'Keys','Row');
vars_selected_summary = sortrows(t(:,[end 1:end-1]),"freq_per_batch","descend","MissingPlacement","last");
end