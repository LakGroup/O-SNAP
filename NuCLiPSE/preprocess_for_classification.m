% function S = preprocess_for_classification(T)


% find cells from groups of interest
if any(ismember({'Group','BioRep'},T.Properties.VariableNames))
    T = renamevars(T,["BioRep","Group"],["biological_replicate","group"]);
end
groups = unique(T.group);
row_idx = cellfun(@(x) strcmpi(T.group,x), groups,'uni',0);
row_idx = any([row_idx{:}],2);
T_norm = T(row_idx,:);

% normalize and remove unnecessary columns
T_norm(:,4:end) = normalize(T_norm(:,4:end));
T_norm = T_norm(:,[1 4:end]);

% save(fullfile(work_dir, analysis_name+".mat"),"T","T_norm");

% end
