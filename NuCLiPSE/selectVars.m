function varSelected = selectVars(T,nVarSel,method)
%% set up
groups = T.Group;
T_num = T(:,vartype('numeric')); % get only numeric columns
labels = replace(T_num.Properties.VariableNames,"_"," ");
X = normalize(table2array(T_num)); % normalize T
% remove nans from normalization (ex. if all values are the same)
validCol = all(~isnan(X));
X = X(:,validCol);
labels = labels(validCol);

%% variable selection
if strcmpi(method,"relieff")
    [varSel_idx,varSel_scores] = relieff(X,cellstr(groups),10,'method','classification');
else
    [varSel_idx,varSel_scores] = feval(lower(method),X,groups);
end
varSelected = labels(varSel_idx(1:nVarSel));
importanceScore = varSel_scores(varSel_idx(1:12));

%% plot
nexttile
b = bar(importanceScore/sum(importanceScore),'FaceColor','flat');
b.CData(1:nVarSel,:) = repmat([0.8500 0.3250 0.0980],nVarSel,1);
ylabel({'Predictor importance score','(probability)'})
ylim([0 0.45])
set(gca,'XTick',1:12, 'XTickLabel',labels(varSel_idx(1:12)))
title('Variable Score')

end