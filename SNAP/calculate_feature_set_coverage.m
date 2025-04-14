function data = calculate_feature_set_coverage(feature_comparisons,save_dir,options)
arguments
    feature_comparisons cell
    save_dir string
    options.plot = true;
    options.alpha double = 0.05;
    options.fold_change_threshold double = 2;
    options.feature_universe_names = ["universe_1"];
end
%% define feature_universes
feature_universes = get_feature_universes(options.feature_universe_names);

for i=1:numel(feature_comparisons)
    group_ctrl = feature_comparisons{i}.groups(1);
    group_case = feature_comparisons{i}.groups(2);
    group_pair_string = replace(group_ctrl + "__vs__" + group_case,"-","_");
    data.(group_pair_string) = sortrows(feature_comparisons{i}.feature_table,"adj_p_value");
    %% run loop on all universes
    for f = 1:size(feature_universes,2)
        data.(group_pair_string) = analyze_feature_family(data.(group_pair_string),[group_ctrl group_case],feature_universes{f},save_dir,options);
    end
end
end

%% perform ontology analysis for a given feature universe
function S = analyze_feature_family(S,group_pair,feature_universe_info,save_dir,options)
    feature_universe_name = feature_universe_info.name;
    feature_universe = feature_universe_info.feature_set;
    %% get feature families (defined in functions below)
    S = sort_features(S,feature_universe,feature_universe_name);
    S_universe = groupcounts(S,feature_universe_name);
    idx_sig = all([abs(S.log2_fold_change) > log2(options.fold_change_threshold) S.adj_p_value < options.alpha],2);
    S_diff = S(idx_sig,:);
    if size(S_diff,1) == 0
        fprintf("   %s: No features satisfy adj p value < %0.2f \n", join(group_pair(:), " vs "), options.alpha);
        return
    end
    S_diff_grouped = sortrows(groupcounts(S_diff,feature_universe_name),"GroupCount","descend");
    n_rows = 0;
    S_down = S_diff(S_diff.fold_change < 0, :);
    if size(S_down,1) > 0
        S_down_grouped = sortrows(groupcounts(S_down,feature_universe_name),"GroupCount","descend");
        n_rows = n_rows + 1;
    end
    S_up = S_diff(S_diff.fold_change > 0,:);
    if size(S_up,1) > 0
        S_up_grouped = sortrows(groupcounts(S_up,feature_universe_name),"GroupCount","descend");
        n_rows = n_rows + 1;
    end
    if n_rows == 2
        n_rows = n_rows + 1;
    end
    % BgRatio = (size of the feature set)/(total number of features)
    for i=1:size(S_universe,1)
        idx = strcmpi(S{:,feature_universe_name},S_universe{i,feature_universe_name});
        S{idx,"bg_ratio_"+feature_universe_name} = repmat(S_universe{i,"Percent"}/100,sum(idx),1);
    end
    % FeatureRatio = (number of hits in data for a specific feature set)/(number of hits in data across all features)
    for i=1:size(S_diff_grouped,1)
        S_diff{strcmpi(S_diff{:,feature_universe_name},S_diff_grouped{i,feature_universe_name}),"feature_ratio_"+feature_universe_name} = sprintf("%d/%d",S_diff_grouped{i,"GroupCount"},sum(S_diff_grouped.GroupCount));
    end
    for i=1:size(S_diff,1)
        idx = strcmpi(S{:,"feature"},S_diff{i,"feature"});
        S{idx,"feature_ratio_"+feature_universe_name} = repmat(S_diff{i,"feature_ratio_"+feature_universe_name},sum(idx),1);
    end
    idx = strcmpi(S{:,"feature_ratio_"+feature_universe_name},"");
    if sum(idx) > 0
        S{idx,"feature_ratio_"+feature_universe_name} = NaN(sum(idx),1);
    end
    %% plot
    if options.plot
        figure('Position',[10 10 910 10+300*n_rows]);
        sgtitle(group_pair(1) + " vs " + group_pair(2))
        tiledlayout(n_rows,2);
        if n_rows == 3
            % counts
            nexttile(5);
            x = replace(S_diff_grouped{:,feature_universe_name},"_"," ");
            x = reordercats(categorical(x),x);
            bar(x,S_diff_grouped.GroupCount);
            ylabel("Counts")
            y_lim_counts = ylim();
            title({"Differential features", "Counts"},"Interpreter","none")
            % coverage (p only)
            nexttile(6);
            S_join = join(S_diff_grouped,S_universe,'Keys',feature_universe_name);
            S_join.coverage = S_join.GroupCount_S_diff_grouped./S_join.GroupCount_S_universe;
            S_join = sortrows(S_join,"coverage","descend");
            x = replace(S_join{:,feature_universe_name},"_"," ");
            x = reordercats(categorical(x),x);
            bar(x,S_join.coverage);
            ylim([0 1])
            xlabel("Feature Families")
            ylabel("Coverage of family (%)")
            title({"Differential features","Coverage"},"Interpreter","none")
        end
        % decreased features in group 2
        if exist("S_down_grouped","var")
            % counts
            nexttile();
            x = replace(S_down_grouped{:,feature_universe_name},"_"," ");
            x = reordercats(categorical(x),x);
            bar(x,S_down_grouped.GroupCount);
            xlabel("Feature Families")
            ylabel("Counts")
            if n_rows == 3
                ylim(y_lim_counts)
            end
            title({"Features increased for " + group_pair(1),"Counts"},"Interpreter","none")
            % coverage (p only)
            nexttile();
            S_join = join(S_down_grouped,S_universe,'Keys',feature_universe_name);
            S_join.coverage = S_join.GroupCount_S_down_grouped./S_join.GroupCount_S_universe;
            S_join = sortrows(S_join,"coverage","descend");
            x = replace(S_join{:,feature_universe_name},"_"," ");
            x = reordercats(categorical(x),x);
            bar(x,S_join.coverage);
            ylim([0 1])
            xlabel("Feature Families")
            ylabel("Coverage of family (%)")
            title({"Features increased for " + group_pair(1),"Coverage"},"Interpreter","none")
        end
        % increased features in group 1
        if exist("S_up_grouped","var")
            % counts
            nexttile();
            x = replace(S_up_grouped{:,feature_universe_name},"_"," ");
            x = reordercats(categorical(x),x);
            bar(x,S_up_grouped.GroupCount);
            xlabel("Feature Families")
            ylabel("Counts")
            if n_rows == 3
                ylim(y_lim_counts)
            end
            title({"Features increased for " + group_pair(2),"Counts"},"Interpreter","none")
            % coverage (p only)
            nexttile();
            S_join = join(S_up_grouped,S_universe,'Keys',feature_universe_name);
            S_join.coverage = S_join.GroupCount_S_up_grouped./S_join.GroupCount_S_universe;
            S_join = sortrows(S_join,"coverage","descend");
            x = replace(S_join{:,feature_universe_name},"_"," ");
            x = reordercats(categorical(x),x);
            bar(x,S_join.coverage);
            ylim([0 1])
            xlabel("Feature Families")
            ylabel("Coverage of family (%)")
            title({"Features increased for " + group_pair(2),"Coverage"},"Interpreter","none")
        end
        sgtitle({replace(feature_universe_name,"_"," "),join(group_pair(:),' vs ')});
        savefig(fullfile(save_dir, "feature_set_coverage_"+join(group_pair(:),'_')+"_"+ feature_universe_name + ".fig"));
        saveas(gcf,fullfile(save_dir, "feature_set_coverage_"+join(group_pair(:),'_')+"_"+ feature_universe_name +".png"));
    end
end

%% sort features into families from the feature universe
function S = sort_features(S,feature_universe,feature_universe_name)
    S{:,feature_universe_name} = repmat("",size(S,1),1);
    feature_families = fieldnames(feature_universe);
    for i=1:size(S,1)
        feature_families_idx = cellfun(@(x) ismember(S{i,"feature"},feature_universe.(x)), feature_families,'uni',1);
        if ~any(feature_families_idx)
            S{i,feature_universe_name} = NaN;
        else
            S{i,feature_universe_name} = feature_families(feature_families_idx);
        end
    end
    if ismissing(S(:,feature_universe_name))
        ME = MException('SNAP:feature_universe_is_missing_features', ...
        'Universe %s does not fully describe features in data',feature_universe_name);
        throw(ME)
    end
end


