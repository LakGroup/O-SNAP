% -------------------------------------------------------------------------
% select_OSNAP_features.m
% -------------------------------------------------------------------------
% Given information on MRMR scores, combines the information to return a
% summary on the selected feature and, when relevant, the result of the
% aggregated feature selection.
%
% Example on how to use it:
%   select_OSNAP_features(scores,vars_names)
% -------------------------------------------------------------------------
% Input:
%   scores: [NxM] array where N is the number of O-SNAP features in which
%           every row corresponds to the score for the same feature. The
%           columns M represent each fold in which the MRMR was performed,
%           if multiple training sets were used to generate multiple sets
%           of feature selections that need to be aggregated.
%   vars_names: Cell array containing information on the feature names that
%               correspond to the rows in scores
% Output:
%   result: A table with information summarizing the scores across the
%           features and folds as well as the final aggreagated score, if
%           calculated
%   vars_selected: Table of boolean values indicating whether a feature is
%                  selected within a batch and from the aggregated
%                  selection. If only a single set of feature selection is
%                  performed (i.e. no folds), then vars_selected is simply
%                  a cell array of the selected feature identifiers 
% Options:
%   max_idx: Upper limit on number of features selected
%   save_path: Name of file location to save to.
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   ....
% -------------------------------------------------------------------------
%%
function [result,vars_selected] = select_OSNAP_features(scores,vars_names,options)
arguments
    scores double
    vars_names cell
    options.max_idx = 12;
    options.save_path string = "";
end
n_batch = size(scores,2);
n_feature = numel(vars_names);

%% batch-wise
if n_batch > 1
    result = table('Size',[n_feature n_batch+1], ...
        'VariableNames',[compose("batch_%02.0f",1:n_batch) "batch_sum"],...
        'VariableTypes',repmat({'logical'},1,n_batch+1));
    result.Properties.RowNames = vars_names';
    vars_selected = table('Size',[n_feature n_batch+1], ...
        'VariableNames',[compose("batch_%02.0f",1:n_batch) "batch_sum"],...
        'VariableTypes',repmat({'logical'},1,n_batch+1));
    vars_selected.Properties.RowNames = vars_names';
    
    % variable selection based on each batch
    scores_sort = zeros(n_feature,n_batch);
    sort_idx = zeros(n_feature,n_batch);
    scores_sort_diff = zeros(options.max_idx-1,n_batch);
    n_vars_selected = zeros(1,n_batch);
    pks_idxs = cell(1,n_batch);
    for i=1:n_batch
        [scores_sort(:,i),sort_idx(:,i)] = sort(scores(:,i),'descend');
        scores_sort_diff(:,i) = diff(scores_sort(1:options.max_idx,i));

        [~,knee_idx] = knee_pt(scores_sort(1:options.max_idx,i));
        [pks,pks_idxs{i}] = findpeaks(abs(scores_sort_diff(:,i)));
        [~,m] = max(pks);
        % defines number of ranked vars to select 
        if pks_idxs{i}(m) < options.max_idx
            pks_idx = pks_idxs{i}(m);
        else
            [~,pks_idx] = max(abs(scores_sort_diff(:,i)));
        end
        n_vars_selected(i) = min([knee_idx, pks_idx]);
        vars_selected_idx = zeros(size(scores,1),1);
        vars_selected_idx(sort_idx(1:n_vars_selected(i),i)) = 1;
        result.(sprintf("batch_%02.0f",i)) = scores(:,i);
        vars_selected.(sprintf("batch_%02.0f",i)) = vars_selected_idx;
    end
    
    % variable selection based on summation over all batches
    [scores_sum_sort,sort_sum_idx] = sort(sum(scores,2),'descend');
    scores_sorted_by_sum = scores(sort_sum_idx,:);
    scores_sum_sort_diff = diff(scores_sum_sort(1:options.max_idx));

    [~,knee_idx] = knee_pt(abs(scores_sum_sort(1:options.max_idx)));
    [pks,pks_sum_idx] = findpeaks(abs(scores_sum_sort_diff));
    [~,m] = max(pks);
    if pks_sum_idx(m) < options.max_idx
        pks_idx = pks_sum_idx(m);
    else
        [~,pks_idx] = max(abs(scores_sum_sort_diff));
    end
    n_vars_selected_sum = min([knee_idx pks_idx]);
    vars_selected_idx = zeros(size(scores,1),1);
    vars_selected_idx(sort_sum_idx(1:n_vars_selected_sum)) = 1;
    result.("batch_sum") = sum(scores,2);
    vars_selected.("batch_sum") = vars_selected_idx;
    
    % plot
    if options.save_path ~= ""
        try
            if exist(options.save_path,'file')
                delete(options.save_path)
            end
            f_batch=figure('visible','off','position',[1,50,1600,800]);
            tiledlayout('flow');
            % each batch
            c_data = lines(n_batch);
            for i=1:n_batch
                nexttile
                bar(scores_sort(1:options.max_idx,i),"FaceColor",c_data(i,:));
                hold on
                plot(1:options.max_idx-1,abs(scores_sort_diff(:,i)),"LineWidth",1.5,"Color","y","LineStyle",":")
                hold on
                scatter(1:n_vars_selected(i),scores_sort(1:n_vars_selected(i),i),'filled','pentagram','SizeData',300,'MarkerFaceColor','yellow','MarkerEdgeColor','k');
                xticklabels(vars_names(sort_idx(1:options.max_idx,i)));
                set(gca,"TickLabelInterpreter",'none')
                set(gca,"FontSize",8)
                ylabel(gca,"MRMR Score")
                % xlabel(gca,"Features","FontSize",14)
                title(sprintf("Batch %02.0f",i))
                legend({'MRMR score','\Delta score','Selected features'})
            end
            axs=findobj(f_batch,'Type','axes');
            linkaxes(axs,'xy');
            % sum
            f_sum=figure('visible','off','position',[1,50,800,600]);
            bar(scores_sorted_by_sum(1:options.max_idx,:),'stacked');
            hold on
            plot(1:options.max_idx-1,abs(scores_sum_sort_diff),"LineWidth",1.5,"Color","y","LineStyle",":")
            hold on
            scatter(1:n_vars_selected_sum,scores_sum_sort(1:n_vars_selected_sum),'filled','pentagram','SizeData',300,'MarkerFaceColor','yellow','MarkerEdgeColor','k');
            xticklabels(vars_names(sort_sum_idx(1:options.max_idx)));
            set(gca,"TickLabelInterpreter",'none')
            set(gca,"FontSize",11)
            ylabel(gca,"Aggregate MRMR Score","FontSize",14)
            xlabel(gca,"Features","FontSize",14)
            title("Batch Sum")
            legend([compose('Batch %.0f score',1:n_batch) {'\Delta score','Selected features'}],...
                'Direction','normal','IconColumnWidth',15)
            % sgtitle('MRMR scores');
            savefig(f_batch,options.save_path+"_batch.fig");
            saveas(f_batch,options.save_path+"_batch.png");
            savefig(f_sum,options.save_path+".fig");
            saveas(f_sum,options.save_path+".png");
        catch ME
            disp(getReport(ME))
        end
    end

%% ALL DATA
else
    result = [];
    % variable selection based on summation over all batches
    [scores_sort,sort_idx] = sort(scores,'descend');
    scores_sorted = scores(sort_idx,:);
    scores_sort_diff = diff(scores_sort(1:options.max_idx));
    [pks,pks_idxs] = findpeaks(abs(scores_sort_diff));
    [~,m] = max(pks);
    if pks_idxs(m) < options.max_idx
        n_vars_selected = pks_idxs(m);
    else
        [~,n_vars_selected] = max(abs(scores_sort_diff));
    end
    vars_selected = vars_names(sort_idx(1:n_vars_selected));

    % plot
    if options.save_path ~= ""
        try
            if exist(options.save_path,'file')
                delete(options.save_path)
            end
            figure('visible','off','position',[1,50,600,400]);
            % sum
            bar(scores_sorted(1:options.max_idx,:),'stacked');
            hold on
            plot(1:options.max_idx-1,abs(scores_sort_diff),"LineWidth",1.5,"Color","y","LineStyle",":")
            hold on
            scatter(pks_idxs,scores_sort(pks_idxs),'filled','MarkerFaceColor','k');
            hold on
            scatter(1:n_vars_selected,scores_sort(1:n_vars_selected),'filled','pentagram','SizeData',300,'MarkerFaceColor','yellow','MarkerEdgeColor','k');
            xticklabels(vars_names(sort_idx(1:options.max_idx)));
            set(gca,"TickLabelInterpreter",'none')
            set(gca,"FontSize",8)
            legend({'Score','\Delta score','\Delta score peaks','Selected features'},...
                'Direction','normal','IconColumnWidth',15)
            % sgtitle('MRMR scores');
            savefig(options.save_path+".fig");
            saveas(gcf,options.save_path+".png");
        catch ME
            disp(getReport(ME))
        end
    end
end
