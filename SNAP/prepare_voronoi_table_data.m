function [T_norm,group_values] = prepare_voronoi_table_data(T,options)
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

% select only desired groups and replicates for analysis
if ~isempty(options.replicates)
    T_norm = T_norm(ismember(T_norm.biological_replicate,options.replicates),:);
end
if ~isempty(options.groups)
    T_norm = T_norm(ismember(T_norm.group,options.groups),:);
end
group_values = T_norm.group;

% normalize T; note that all information stored in signs is erased
if options.normalize
    T_norm{:,vartype('numeric')} = normalize(abs(T_norm{:,vartype('numeric')}),1);
end

% remove NaNs
if options.remove_NaN
    is_nan_col = all(ismissing(T_norm),1);
    T_norm = T_norm(:,~is_nan_col);
    is_nan_row = any(ismissing(T_norm),2);
    T_norm = T_norm(~is_nan_row,:);
    group_values = group_values(~is_nan_row);
end

% remove unnecessary columns (group, bio replicate)
if options.numeric_only
    T_norm = T_norm(:,4:end);
elseif ~options.keep_rep_sample_info
    T_norm = T_norm(:,[1 4:end]);
end
end
