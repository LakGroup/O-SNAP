function T_norm = prepare_voronoi_table_data(T,options)
arguments
    T table
    options.groups cell = {};
    options.replicates cell = {};
    options.normalize logical = true;
    options.remove_NaN logical = true;
    options.keep_rep_sample_info logical = false;
    options.numeric_only logical = false;
end

% find cells from groups of interest
if ~isempty(options.groups)
    in_groups_idx = cellfun(@(x) strcmpi(T.group,x), options.groups,'uni',0);
    in_groups_idx = any([in_groups_idx{:}],2);
else
    in_groups_idx = ones(size(T,1),1);
end
if ~isempty(options.replicates)
    in_reps_idx = cellfun(@(x) strcmpi(T.biological_replicate,x), options.replicates,'uni',0);
    in_reps_idx = any([in_reps_idx{:}],2);
else
    in_reps_idx = ones(size(T,1),1);
end
idx = all([in_groups_idx in_reps_idx],2);
T_norm = T(idx,:);

% normalize
if options.normalize
    T_norm(:,4:end) = normalize(T_norm(:,4:end));
end

% remove NaNs
if options.remove_NaN
    is_nan_col = all(ismissing(T_norm),1);
    T_norm = T_norm(:,~is_nan_col);
    is_nan_row = any(ismissing(T_norm),2);
    T_norm = T_norm(~is_nan_row,:);
end

% remove unnecessary columns (group, bio replicate)
if options.numeric_only
    T_norm = T_norm(:,4:end);
elseif ~options.keep_rep_sample_info
    T_norm = T_norm(:,[1 4:end]);
end
end