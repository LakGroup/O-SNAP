function [T_norm,group_values] = preprocess_SNAP_table(T,options)
arguments
    T table
    options.groups cell = {};
    options.replicates cell = {};
    options.normalize logical = true;
    options.remove_NaN logical = true;
    options.keep_rep_sample_info logical = false;
    options.numeric_only logical = false; % only retains numeric columns
end

T_norm = T;

% remove NaNs
if options.remove_NaN
    is_nan_col = all(ismissing(T_norm),1);
    T_norm = T_norm(:,~is_nan_col);
    is_nan_row = any(ismissing(T_norm),2);
    T_norm = T_norm(~is_nan_row,:);
end

% select only desired groups and replicates for analysis
if ~isempty(options.replicates)
    T_norm = T_norm(ismember(T_norm.biological_replicate,options.replicates),:);
end
if ~isempty(options.groups)
    T_norm = T_norm(ismember(T_norm.group,options.groups),:);
end
if ismember('group',T_norm.Properties.VariableNames)
    group_values = T_norm.group;
else
    group_values = [];
end

% normalize T
if options.normalize
    % T_norm{:,vartype('numeric')} = normalize(abs(T_norm{:,vartype('numeric')}),1);
    T_norm{:,vartype('numeric')} = normalize(T_norm{:,vartype('numeric')},1);
end

% remove unnecessary columns (group, bio replicate)
if options.numeric_only
    T_norm = T_norm(:,vartype('numeric'));
elseif ~options.keep_rep_sample_info
    if ismember('biological_replicate',T_norm.Properties.VariableNames)
        T_norm(:,'biological_replicate') = [];
    end
    if ismember('name',T_norm.Properties.VariableNames)
        T_norm(:,'name') = [];
    end
end
end
