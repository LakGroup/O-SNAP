function [vars_selected, var_select_result] = select_features_SNAP(T,options)
arguments
    T table
    options.methods cell = {...
        'biPLS',...
        'rPLS',...
        'ChiSquare',...
        'MRMR',...
        'ReliefF',...
        'Boruta',...
        'biPLS_PLSToolbox',... %(Interval PLS from the PLS Toolbox; Eigenvector Inc.)
        'rPLS_PLSToolbox',... %(Recursive PLS from the PLS Toolbox; Eigenvector Inc.)
        'iPLS_PLSToolbox',... %(Interval PLS from the PLS Toolbox; Eigenvector Inc.)
        'GA_PLSToolbox'... %(Genetic Algorithm from the PLS Toolbox; Eigenvector Inc.)
        }
    options.threshold double = 0.5
    options.verbose logical = true
end

%% set up
[T_num,groups] =  prepare_voronoi_table_data(T,"numeric_only",true,"normalize",true); % get only numeric columns
X = table2array(T_num);
[~,~,groups_idx] = unique(groups);
var_select_result = cell(1,numel(options.methods));
var_select_summary = zeros(1,size(X,2));
n_selection_methods = 0;
methods = options.methods;
%% variable selection
parfor i=1:length(options.methods)
    try
        var_select_result{i} = VariableSelection(X,groups_idx,methods{i});
        if ~isempty(var_select_result{i}.VarSel)
            var_select_summary = var_select_summary + var_select_result{i}.VarSel;
            n_selection_methods = n_selection_methods + 1;
        end
    catch ME
        disp(getReport(ME))
    end    
end
%% find features that satisfy frequency threshold across different selection method   
ratio_found = var_select_summary./(n_selection_methods*ones(1,size(X,2)));
idx_selected = find(ratio_found>options.threshold);
vars_selected = T_num.Properties.VariableNames(idx_selected);
ratio_found_selected = ratio_found(idx_selected);
[ratio_found_selected_sorted,idx_selected_sorted] = sort(ratio_found_selected,'descend');
vars_selected = vars_selected(idx_selected_sorted);
if options.verbose
    if numel(idx_selected) > 0
        fprintf("    Features selected: %2.0f\n",numel(idx_selected))
        for i=1:numel(idx_selected)
            fprintf("      - %s (%2.1f%%)\n",vars_selected{i},ratio_found_selected_sorted(i)*100);
        end
    else
        disp("    No features found")
    end
end
end