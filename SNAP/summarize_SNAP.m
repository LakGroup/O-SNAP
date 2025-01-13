close all force
rng("default") % reproducibility

% 
root_dir = "F:\";
analysis_name = "ANALYSIS_HK_BJ_TSA_H2B";
% analysis_name = "ANALYSIS_AM_Heterokaryon_H2B";
work_dir = fullfile(root_dir,analysis_name);
image_dir = fullfile(root_dir,"ANALYSIS_HK_BJ_TSA_H2B");
% image_dir = fullfile(root_dir,"ANALYSIS_AM_Heterokaryon_H2B");
save_analysis_path = fullfile(work_dir,analysis_name+".mat");
save_classification_path = fullfile(work_dir,analysis_name+"_classifiers.mat");

%% load data
loaded_data = load_variables(save_analysis_path,{...
    'groups',...
    'replicates',...
    'feature_comparisons'...
    'feature_data',...
    'var_select_result',...
    'vars_selected',...
    'pca_result'});
groups = loaded_data.groups;
replicates = loaded_data.replicates;
feature_comparisons = loaded_data.feature_comparisons;
feature_data = loaded_data.feature_data;
var_select_result = loaded_data.var_select_result;
vars_selected = loaded_data.vars_selected;
pca_result = loaded_data.pca_result;
classifier_data = load_variables(save_classification_path,{...
    'classifiers',...
    'classification_summary'});
classifiers = classifier_data.classifiers;
classification_summary = classifier_data.classification_summary;
% clearvars loaded_data classifier_data
n_groups = numel(groups);
n_compare = numel(feature_comparisons);
n_each_group = groupcounts(feature_data.group);
SNAP_nucleus_file_list = get_valid_SNAP_nucleus_files(image_dir,groups,replicates,{'voronoi_areas'});
%% generate feature selection table
vars_selected_methods = cellfun(@(x) x.Method,var_select_result,'uni',0);
n_vars_selected_methods = numel(vars_selected_methods);
n_vars_selected = numel(vars_selected);
vars_all_table = array2table(reshape(cell2mat(cellfun(@(x) x.VarSel,var_select_result,'uni',0)),[],n_vars_selected_methods),...
    'VariableNames',vars_selected_methods,...
    'RowNames',feature_data.Properties.VariableNames(4:end));
vars_all_table.Properties.VariableTypes=repmat({'logical'},1,n_vars_selected_methods);
vars_all_table.Ratio = table2array(sum(vars_all_table,2))/n_vars_selected_methods;
vars_all_table = vars_all_table(:,[end end-1 1:end-2]);
%% SET UP FIGURE
fig_width = 1920;
fig_height = 920;
fig = uifigure('Name',analysis_name +" - SNAP summary",'Position',[100 100 fig_width fig_height]);
g = uigridlayout(fig);
g.ColumnWidth = {'1x','1x'};
g.RowHeight = {2*fig_height/5,3*fig_height/5,'1x'};

%% REPRESENTATIVE IMAGES
p_rep_images = uipanel(g);
g_rep_images = uigridlayout(p_rep_images,[4,2]);
g_rep_images.RowHeight = {30,20,'1x',60};
g_rep_images.ColumnWidth =  repmat({'1x'},1,n_groups);
title = uilabel(g_rep_images,'Text','REPRESENTATIVE IMAGES');
title.Layout.Column = [1 n_groups];
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
for i=1:n_groups
    title = uilabel(g_rep_images,'Text',sprintf('%s (N=%.0f nuclei)',groups{i},n_each_group(i)));
    title.HorizontalAlignment = 'center';
    title.FontSize = 18;
end
filepaths_by_groups = cellfun(@(x) replace(replace(SNAP_nucleus_file_list{SNAP_nucleus_file_list.group == x,'filepath'},"SNAP_nucleus_data","SNAP_nucleus_images"),".mat",".png"),groups,'uni',0);
for i=1:n_groups
    im = uiimage(g_rep_images);
    im.ImageSource = filepaths_by_groups{i}(1);
end
for i=1:n_groups
    im = g_rep_images.Children(1+n_groups+i);
    g_path_label = uigridlayout(g_rep_images,[1,3]);
    g_path_label.ColumnWidth =  {50,'1x',50};
    button_prev = uibutton(g_path_label,...
        "Text",sprintf("Prev<<\n(%3.0f)",0),...
        "Enable","off");            
    path_label = uilabel(g_path_label);
    c = uicontrol(fig,'Style','text');
    path_label.Text = textwrap(c,{im.ImageSource},strlength(extractBefore(im.ImageSource,'SNAP_nucleus_images'))+14);
    path_label.HorizontalAlignment = 'center';
    path_label.FontSize = 11;
    button_next = uibutton(g_path_label,...
        "Text",sprintf("Next>>\n(%3.0f)",numel(filepaths_by_groups{i})-1));
    button_prev.ButtonPushedFcn = @(src,event) rep_image_prev_callback(fig,button_prev,button_next,im,path_label,filepaths_by_groups{i},im.ImageSource);
    button_next.ButtonPushedFcn = @(src,event) rep_image_next_callback(fig,button_prev,button_next,im,path_label,filepaths_by_groups{i},im.ImageSource);
end
%% VOLCANO PLOT
p_volcano = uipanel(g);
g_volcano = uigridlayout(p_volcano,[3,1]);
g_volcano.RowHeight = {30,20,'1x'};
title = uilabel(g_volcano,'Text','VOLCANO PLOT');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
feature_comparisons_list = cellfun(@(x) join(x.groups," vs "),feature_comparisons,'uni',0);
feature_comparisons_list = [feature_comparisons_list{:}];
gt_volcano = uigridlayout(g_volcano,[1,2]);
gt_volcano.Layout.Row = 3;
dd_volcano = uidropdown(g_volcano,"Items",feature_comparisons_list,...
    "ValueChangedFcn",@(src,event) plot_volcano_table(src,fig,gt_volcano,feature_comparisons),...
    "FontWeight","bold",...
    "FontSize",14);
dd_volcano.Layout.Row = 2;
gt_volcano.ColumnWidth = {'2x','3x'};
plot_volcano_table(dd_volcano,fig,gt_volcano,feature_comparisons,vars_all_table);
%% FEATURE SELECTION
vars_selected_table = vars_all_table(ismember(vars_all_table.Properties.RowNames,vars_selected),:);
vars_selected_table = sortrows(vars_selected_table,"Ratio","descend");
p_feature_select = uipanel(g);
g_feature_select = uigridlayout(p_feature_select,[1 2]);
g_feature_select.RowHeight = {30,'1x'};
g_feature_select.ColumnWidth =  {'1x','1.5x'};
title = uilabel(g_feature_select,'Text','FEATURE SELECTION');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
title.Layout.Row = 1;
title.Layout.Column = [1 2];
vars_selected_uitable = uitable(g_feature_select,...
    'Data',vars_selected_table,...
    'SelectionType','row',...
    'Multiselect','on');
p_violin=uipanel(g_feature_select,'AutoResizeChildren','off');
g_violin = uigridlayout(p_violin,[2,1]);
create_violin_panel(g_violin,vars_selected_table,feature_data,groups)
vars_selected_uitable.SelectionChangedFcn = @(src,event) update_feature_selections(src,event,g_violin,feature_data);
vars_selected_uitable.Selection = 1:4;
% format cells
[ratio,~,ratio_idx] = unique(vars_selected_table{:,1});
for i=1:numel(ratio)
    idx = ratio_idx==i;
    cell_color = [0 1 0]+([1 0 1]*(1-ratio(i))*2);
    addStyle(vars_selected_uitable,uistyle("BackgroundColor",cell_color),"cell",[find(idx),ones(sum(idx),1)]);
end
% adjust first column
vars_selected_uitable.Data{:,1}=vars_selected_uitable.Data{:,1}*100;
vars_selected_uitable.ColumnName{1} = "Selected Freq.";
%% PCA/CLASSIFICATION
g_PCA_classification = uigridlayout(g,[1,2]);
g_PCA_classification.ColumnWidth = {'1x','1x'};
%% PCA
p_PCA = uipanel(g_PCA_classification);
g_PCA = uigridlayout(p_PCA,[2,1]);
g_PCA.RowHeight = {30,'1x'};
title = uilabel(g_PCA,'Text','PRINCIPAL COMPONENT ANALYSIS');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
g_PCA_ax = axes(g_PCA);
plot_PCA(g_PCA_ax,pca_result.pca_coefficients,pca_result.pca_scores,pca_result.explained,"groups_info",feature_data.group,"bio_reps_info",feature_data.biological_replicate);
%% CLASSIFICATION
p_classification = uipanel(g_PCA_classification);
g_classification = uigridlayout(p_classification,[4, 1]);
g_classification.RowHeight =  {30,16,'1x','3x'};
title = uilabel(g_classification,'Text','CLASSIFICATION');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
title = uilabel(g_classification,'Text',sprintf('Validation accuracies (N=%.0f models)',size(classifiers,2)));
title.HorizontalAlignment = 'center';
title.FontSize = 12;
g_classification_table = uitable(g_classification,...
    'Data',classification_summary{:,2:end},...
    'RowName',classification_summary{:,1},...
    'ColumnName',{'Mean','Std'},...
    'SelectionType','row');
ax_confusion = axes(g_classification);
plot_confusion(ax_confusion,classifiers(1,1:end)) % only top classifiers
g_classification_table.SelectionChangedFcn = @(~,event) update_classification_selections(event,ax_confusion,classifiers);
%% SAVE
% exportgraphics(fig,fullfile(work_dir,analysis_name+"_summary.pdf"));
% disp("Results saved to: "+analysis_name+"_summary.pdf")
%% Previous image in file list
function rep_image_prev_callback(fig, button_prev,button_next,ax_image,ax_label,filepaths,current_path)
    current_idx = find(strcmp(filepaths,current_path));
    if current_idx > 1
        try
            button_prev.Enable = 'on';
            button_prev.Enable = 'on';
            new_idx = current_idx-1;
            c = uicontrol(fig,'Style','text');
            ax_image.ImageSource = filepaths(new_idx);
            ax_label.Text = textwrap(c,filepaths(new_idx),strlength(extractBefore(filepaths(new_idx),'SNAP_nucleus_images'))+14);           
            button_prev.Text = sprintf("<<Prev\n(%3.0f)",new_idx);
            button_next.Text = sprintf("Next>>\n(%3.0f)",numel(filepaths)-new_idx);
        catch
        end
    else
        button_prev.Enable = 'off';
    end
end
function rep_image_next_callback(fig, button_prev,button_next,ax_image,ax_label,filepaths,current_path)
    current_idx = find(strcmp(filepaths,current_path));
    if current_idx < numel(filepaths)
        try
            button_prev.Enable = 'on';
            button_next.Enable = 'on';
            new_idx = current_idx+1;
            c = uicontrol(fig,'Style','text');
            ax_image.ImageSource = filepaths(new_idx);
            ax_label.Text = textwrap(c,filepaths(new_idx),strlength(extractBefore(filepaths(new_idx),'SNAP_nucleus_images'))+14);           
            button_prev.Text = sprintf("<<Prev\n(%3.0f)",new_idx);
            button_next.Text = sprintf("Next>>\n(%3.0f)",numel(filepaths)-new_idx);
        catch ME
            disp(getReport(ME))
        end
    else
        button_next.Enable = 'off';
    end
end
function plot_volcano_table(src,fig,gridlayout,feature_comparisons,vars_selection_table)
    idx = src.ValueIndex;
    comparison_table = feature_comparisons{idx}.feature_table;
    comparison_table.("Selected Freq.") = vars_selection_table{:,1};
    feature_comparisons{idx}.feature_table;
    comparison_table = sortrows(comparison_table,"log2_fold_change","descend");
    comparison_table = comparison_table(:,...
        ["feature","Selected Freq.","log2_fold_change","adj_p_value",comparison_table.Properties.VariableNames(2:end-4),"fold_change","p_value"]);
    if numel(gridlayout.Children) < 2
        g_volcano_axes = axes(gridlayout);
        g_volcano_table = uitable(gridlayout);
    else
        g_volcano_axes = gridlayout.Children(1);
        g_volcano_table = gridlayout.Children(2);
    end
    plot_volcano(fig,g_volcano_axes,feature_comparisons{idx});
    g_volcano_table.Data=comparison_table{:,2:end};
    g_volcano_table.ColumnName=comparison_table.Properties.VariableNames(2:end);
    g_volcano_table.RowName=comparison_table.feature;
    g_volcano_table.ColumnSortable = true;
    % highlight red
    highlight_idx = all([(g_volcano_table.Data(:,2)>1 ) (g_volcano_table.Data(:,3) < 0.05)],2);
    addStyle(g_volcano_table,uistyle("BackgroundColor",'#FF6262'),"row",find(highlight_idx));
    % highlight blue
    highlight_idx = all([(g_volcano_table.Data(:,2)<-1 ) (g_volcano_table.Data(:,3) < 0.05)],2);
    addStyle(g_volcano_table,uistyle("BackgroundColor",'#62C8FF'),"row",find(highlight_idx));
    % highlight green
    [ratio,~,ratio_idx] = unique(g_volcano_table.Data(:,1));
    for i=1:numel(ratio)
        idx = ratio_idx==i;
        if ratio(i) >= 0.5
            cell_color = [0 1 0]+([1 0 1]*(1-ratio(i))*2);
        else
            cell_color = [1 1 1];
        end
        addStyle(g_volcano_table,uistyle("BackgroundColor",cell_color),"cell",[find(idx),ones(sum(idx),1)]);
    end
end
%% generate volcano plot
function plot_volcano(fig,ax,feature_comparisons,options)
    arguments
        fig
        ax
        feature_comparisons struct;
        options.alpha double = 0.05;
        options.fold_change_threshold double = 2;
        options.plot_labels logical = true;
    end
    feature_comparison_table = feature_comparisons.feature_table;
    group_1 = feature_comparisons.groups(1);
    group_2 = feature_comparisons.groups(2);
    log2_fold_change_threshold = log2(options.fold_change_threshold);
    idx_sig = all([abs(feature_comparison_table.log2_fold_change)>log2_fold_change_threshold feature_comparison_table.adj_p_value < options.alpha],2);
    idx_sig_up = all([feature_comparison_table.log2_fold_change>0 idx_sig],2);
    idx_sig_down = all([feature_comparison_table.log2_fold_change<0 idx_sig],2);
    ax_size = 20*ones(size(feature_comparison_table,1),1);
    ax_size(idx_sig) = 35*ones(sum(idx_sig),1);
    ax_color = repmat([0.5 0.5 .5],size(feature_comparison_table,1),1);
    ax_color(idx_sig_up,:) = repmat([1 0 0],sum(idx_sig_up),1);
    ax_color(idx_sig_down,:) = repmat([0 0 1],sum(idx_sig_down),1);
    scatter(ax,feature_comparison_table.log2_fold_change,-log10(feature_comparison_table.adj_p_value),ax_size,ax_color,'filled','MarkerEdgeColor','k');
    hold(ax,'on')
    plot(ax,[log2_fold_change_threshold log2_fold_change_threshold],ylim(ax),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    hold(ax,'on')
    plot(ax,-[log2_fold_change_threshold log2_fold_change_threshold],ylim(ax),'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    hold(ax,'on')
    plot(ax,xlim(ax), [-log10(options.alpha) -log10(options.alpha)],'LineStyle','-.','Color',[1 0.3 0],'LineWidth',1.2)
    ylabel(ax,"-log_{10}(adj. p-value)")
    xlabel(ax,"log_{2}(Fold Change)")
    title(ax,replace("\color{blue}" + group_1 + " \color{black}vs " + "\color{red}" + group_2,"_","-"));
    if options.plot_labels
        hold(ax,'on')
        x = feature_comparison_table{idx_sig,"log2_fold_change"}+0.025;
        y = -log10(feature_comparison_table{idx_sig,"adj_p_value"});
        feature_name = replace(feature_comparison_table{idx_sig,"feature"},"_"," ");
        for f=1:size(feature_name,1)
            if x(f) < 0
                t_color = 'b';
            elseif x(f) > 0
                t_color = 'r';
            else
                t_color = 'k';
            end
            ht = text(ax,x(f),y(f),'    ','HorizontalAlignment','center',...
                'tag',num2str(f),'Color',t_color,'FontSize',10,'Interpreter','none');
            pointerBehavior.enterFcn = ...
                @(ax, cpp)set(findobj(fig,'tag',num2str(f)), ...
                'string', {feature_name(f),""});
            pointerBehavior.traverseFcn = [];
            pointerBehavior.exitFcn = ...
                @(ax, cpp)set(findobj(fig, 'tag',num2str(f)), ...
                'string', '    ');
            iptSetPointerBehavior(ht, pointerBehavior);
            iptPointerManager(fig, 'enable');
            hold(ax,'on')
        end
    end
end
function create_violin_panel(gridlayout,vars_selected_table,feature_data,groups)
    ax_violin_legend = axes(gridlayout);
    n_legend_rows = create_violin_legend(ax_violin_legend,groups);
    gridlayout.RowHeight = {30*n_legend_rows,'1x'};
    p_violin_axes = uipanel(gridlayout,"BorderWidth",0);
    p_violin_axes.AutoResizeChildren = 'off';
    ax_violin_plot = tiledlayout(p_violin_axes,2,2,'TileSpacing','Compact','Padding','Compact');
    % plot new axes
    for i=1:4
        ax_violin_plot_i = nexttile(ax_violin_plot);
        plot_violin(ax_violin_plot_i,feature_data,vars_selected_table.Properties.RowNames{i});
    end
end
%% create one legend for violin plots
function n_legend_rows = create_violin_legend(ax,groups)
    set(ax,'ydir','reverse');
    n_groups = numel(groups);
    c_groups = lines(numel(groups));
    x0 = 32;
    text_ratio = 1.4;
    x_center = x0/2;
    y_center = 1/2;
    icon_length = 0.6;
    offset_x = x0/100;
    x_min = offset_x;
    text_graphics = cell(1,n_groups+2);
    for i=1:n_groups
        text_graphics{i}=text(ax,x_min,y_center,groups{i},'FontSize',10);
    end
    text_graphics{n_groups+1}=text(ax,x_min,y_center,'Mean','FontSize',10);
    text_graphics{n_groups+2}=text(ax,x_min,y_center,'Median','FontSize',10);
    n_legend_rows = ceil(sum([cellfun(@(t) t.Extent(3)/t.Extent(4)*text_ratio,text_graphics,'uni',1),...
        icon_length*(n_groups+1),offset_x*(n_groups+3+2)])/x0);
    daspect(ax,[1 1 1]);
    xlim(ax,[0 x0]);
    ylim(ax,[0 n_legend_rows]);
    x = zeros(1,n_legend_rows);
    row_idx = 1;
    graphics_idx = 1;
    graphics_by_row = cell(1,n_legend_rows);
    for i=1:n_groups+2
        % calculate length of legend entry
        t_length = text_graphics{i}.Extent(3)/text_graphics{i}.Extent(4)*text_ratio;
        if i <= n_groups
            width_i = icon_length + t_length + offset_x;
        elseif i==n_groups+1
            width_i = t_length + 2*offset_x;
        else
            width_i = icon_length + t_length + 2*offset_x;
        end
        % create new row if necessary
        if (x(row_idx) + width_i) > x0
            row_idx = row_idx + 1;
            graphics_idx = 1;
        end
        % add legend entry
        if i <= n_groups
            icon = rectangle(ax,'Position',[x(row_idx) 0.15+(row_idx-1) icon_length 0.7],...
                'FaceColor',c_groups(i,:),'LineWidth',0.1);
            hold(ax,'on');
            text_graphics{i}.Position(1) = x(row_idx) + icon_length + offset_x;
            text_graphics{i}.Position(2) = text_graphics{i}.Position(2) + row_idx - 1;
        elseif i == n_groups+1
            icon = scatter(ax,x(row_idx)+0.5*offset_x,y_center+(row_idx-1),50,'pentagramk','filled');
            hold(ax,'on');
            text_graphics{i}.Position(1) = x(row_idx) + 2*offset_x;
            text_graphics{i}.Position(2) = text_graphics{i}.Position(2) + row_idx - 1;
        else
            icon = line(ax,[x(row_idx) (x(row_idx) + icon_length)],[y_center y_center]+(row_idx-1),"Color","k","LineWidth",2);
            hold(ax,'on');
            text_graphics{i}.Position(1) = x(row_idx) + icon_length + offset_x;
            text_graphics{i}.Position(2) = text_graphics{i}.Position(2) + row_idx - 1;
        end
        % update final row x position and graphics by row for next legend entry
        x(row_idx) = x(row_idx) + width_i;
        graphics_by_row{row_idx}(graphics_idx) = icon;
        graphics_by_row{row_idx}(graphics_idx+1) = text_graphics{i};
        graphics_idx = graphics_idx + 2;
    end
    drawnow()
    for row_idx=1:n_legend_rows
        row_center = x(row_idx)/2;
        x_adjust = x_center - row_center;
        for graphics_idx=1:numel(graphics_by_row{row_idx})
            graphics = graphics_by_row{row_idx};
            if any(strcmp(class(graphics(graphics_idx)),["matlab.graphics.primitive.Line","matlab.graphics.primitive.Data","matlab.graphics.chart.primitive.Scatter"]))
                graphics(graphics_idx).XData = graphics(graphics_idx).XData + x_adjust;
            else
                graphics(graphics_idx).Position(1) = graphics(graphics_idx).Position(1) + x_adjust;
            end
        end
    end
    set(ax,'Visible','off');
end
%% update violin plot panel
function update_feature_selections(src,event,panel,feature_data)
    idx_select = event.Selection;
    if numel(idx_select) > 4
        src.Selection = idx_select(1:4);
        n_rows = 4;
    else
        n_rows = numel(idx_select);
    end
    p_violin_axes = panel.Children(2);
    % reset axes
    for i=1:numel(p_violin_axes.Children)
        delete(p_violin_axes.Children(1));
    end
    % plot new axes
    for i=1:n_rows
        ax_violin_plot_i = subplot(2,2,i,'Parent',p_violin_axes);
        plot_violin(ax_violin_plot_i,feature_data,event.Source.RowName{idx_select(i)});
    end
end
%% generate violin plot
function plot_violin(ax,T,feature)
arguments
    ax
    T table
    feature string
end
% get info on groups
[groups,~,group_idx] = unique(T.group);
group_cat = categorical(T.group);
n_groups = numel(groups);
% create violinplot
c_map = lines(n_groups);
feature_data = T{:,feature};
for g=1:n_groups
    idx = group_idx==g;
    group_idx_g = group_cat(idx);
    feature_data_g = feature_data(idx);
    [feature_data_g_no_outliers, outliers_idx] = rmoutliers(feature_data_g);
    group_idx_g = group_idx_g(~outliers_idx);
    % violin
    v = violinplot(ax,group_idx_g,feature_data_g_no_outliers);
    v.FaceColor = c_map(g,:);
    hold(ax,"on")
    % included points
    swarmchart(ax,group_idx_g,feature_data_g_no_outliers,10,c_map(g,:),'filled');
    hold(ax,"on")
    % median
    median_val = median(feature_data_g,"omitmissing");
    plot(ax,[g-0.4 g+0.4], [median_val median_val],"-k","LineWidth",2);
    hold(ax,"on")
    % mean - NOTE: OUTLIERS ARE REMOVED WHEN CALCULATING MEAN, BUT
    % NOT MEDIAN
    scatter(ax,g,mean(feature_data_g_no_outliers,"omitmissing"),50,'pentagramk','filled')
    hold(ax,"on")
    % outliers
    if sum(outliers_idx > 0)
        scatter(ax,repmat(g,1,sum(outliers_idx)),feature_data_g(outliers_idx),"xk");
    end
    ax.XTick = [];
end
scatter(ax,nan,nan,50,'pentagramk','filled');
plot(ax,nan, nan,"-k","LineWidth",2);
title(ax,replace(feature,"_"," "))
end
%% update confusion matrix
function update_classification_selections(event,ax,classifiers)
    idx_select = event.Selection;
    plot_confusion(ax,classifiers(idx_select,:))
end
%% plot confusion matrix
function plot_confusion(ax,classifiers)
    confusion_data = cellfun(@(x) confusionmat(x.ValidationResults.Response,x.ValidationResults.Predictions), classifiers,'uni',0);
    confusion_data = reshape(cell2mat(confusion_data),2,2,[]);
    confusion_mean = mean(confusion_data,3);
    confusion_std = std(confusion_data,0,3);
    groups = sort(unique(classifiers{1}.ValidationResults.Response));
    n_groups = numel(groups);
    set(ax,'ydir','reverse')
    xlim(ax,[0 n_groups])
    ylim(ax,[0 n_groups])
    xticks(ax,(1:n_groups)-0.5)
    yticks(ax,(1:n_groups)-0.5)
    ax.FontSize = 14;
    xticklabels(ax, groups);
    yticklabels(ax, groups);
    ytickangle(ax,90);
    ax.XLabel.FontSize = 12;
    ax.YLabel.FontSize = 12;
    xlabel(ax,'Predicted Class','FontSize',14)
    ylabel(ax,'True Class','FontSize',14,'Rotation',90)
    color_diagonal = [0.00,0.45,0.74];
    color_offdiagonal = [0.85,0.33,0.10];
    title(ax,replace(classifiers{1}.ModelType,"_"," "),'FontSize',12)
    for i_true=1:n_groups
        confusion_rowsum = sum(confusion_data(i_true,:,1));
        for i_pred=1:n_groups
            confusion_ratio = confusion_mean(i_pred,i_true)/confusion_rowsum;
            if i_true == i_pred
                c = color_diagonal+(1-color_diagonal)*(1-confusion_ratio);
            else
                c = color_offdiagonal+(1-color_offdiagonal)*(1-confusion_ratio);
            end
            if confusion_ratio > 0.5
                t = 'w';
            else
                t = 'k';
            end
            rectangle(ax,...
                'Position',[i_pred-1 i_true-1 1 1],...
                'FaceColor',c)
            text(ax,i_pred-0.5,i_true-0.5-0.08,...
                sprintf('%.0f \x00B1 %.0f',...
                confusion_mean(i_pred,i_true),...
                confusion_std(i_pred,i_true)),...
                'HorizontalAlignment','center',...
                'FontSize',12,...
                'Color',t);
            text(ax,i_pred-0.5,i_true-0.5+0.08,...
                sprintf('(%.1f%% \x00B1 %.1f%%)',...
                confusion_ratio*100,...
                confusion_std(i_pred,i_true)/confusion_rowsum*100),...
                'HorizontalAlignment','center',...
                'FontSize',10,...
                'FontAngle','italic',...
                'Color',t);
        end
    end
    disableDefaultInteractivity(ax)
end
%% Plot PCA
function plot_PCA(ax,coeff,score,explained,options)
arguments
    ax
    coeff double
    score double
    explained double
    options.bio_reps_info string = []
    options.groups_info string = []
    options.color_by_group logical = true
    options.show_bio_reps logical = true
    options.marker_size = 20
end
    if size(score,2)<3
        return
    end
    coeff = coeff(:,1:3);
    score = score(:,1:3);
    %% plot
    marker_list = {'o','^','square','v','diamond','pentagram','hexagram'};
    if ~isempty(options.groups_info)
        [groups,~,groups_idx] = unique(options.groups_info);
    end
    if options.show_bio_reps && ~isempty(options.bio_reps_info)
        [bio_reps,~,bio_reps_idx] = unique(options.bio_reps_info);
    end
    if ~options.color_by_group && exist(bio_reps,"var")
        s_face_colors = lines(length(bio_reps));
    elseif ~isempty(groups)
        s_face_colors = lines(length(groups));
    else
        s_face_colors = [0 0 0];
    end
    % Force each column of the coefficients to have a positive largest element.
    % This tends to put the large var vectors in the top and right halves of
    % the plot.
    [p,d] = size(coeff);
    [~,maxind] = max(abs(coeff),[],1);
    colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
    coeff = coeff .* colsign;
    % Scale the scores so they fit on the plot, and change the sign of
    % their coordinates according to the sign convention for the coefs.
    max_coef_len = sqrt(max(sum(coeff.^2,2)));
    score = (max_coef_len.*(score ./ max(abs(score(:))))).*colsign;
    % only group information is availble
    if ~exist("bio_reps","var") && ~isempty(groups)
        for g=1:length(groups)
            idx = groups_idx == g;
            s = scatter3(ax,score(idx,1),score(idx,2),score(idx,3),options.marker_size,'filled');
            s.MarkerFaceColor = s_face_colors(g,:);
            hold(ax,'on')
        end
        legend(ax,groups,'Location','southoutside','Interpreter','none')
    % color = groups; marker type = replicate
    elseif options.color_by_group && ~isempty(groups)
        for b=1:length(bio_reps)
            for g=1:length(groups)
                idx = all([(groups_idx == g) (bio_reps_idx == b)],2);
                s = scatter3(ax,score(idx,1),score(idx,2),score(idx,3),'filled');
                s.MarkerFaceColor = s_face_colors(g,:);
                s.Marker = marker_list(b);
                if b == 1
                    s.SizeData = options.marker_size;
                else
                    s.SizeData = 2*options.marker_size;
                end
                hold(ax,'on')
            end
        end
        legend(ax,groups,'Location','southoutside','Interpreter','none')
    % color = replicate; marker type = color
    elseif ~isempty(groups)
        for b=1:length(bio_reps)
            for g=1:length(groups)
                idx = all([(groups_idx == g) (bio_reps_idx == b)],2);
                s = scatter3(ax,score(idx,1),score(idx,2),score(idx,3),'filled');
                s.MarkerFaceColor = s_face_colors(b,:);
                s.Marker = marker_list(g);
                if g == 1
                    s.SizeData = options.marker_size;
                else
                    s.SizeData = 2*options.marker_size;
                end
                hold(ax,'on')
            end
        end
        legend(ax,bio_reps,'Location','southoutside','Interpreter','none')
    else 
        scatter3(ax,score(:,1),score(:,2),score(:,3),'filled');
    end
    % biplot(ax,coeff,'VarLabels',T_num.Properties.VariableNames);
    grid(ax,'on');
    set(ax, "XTickLabel",[]);
    set(ax, "YTickLabel",[]);
    set(ax, "ZTickLabel",[]);
    xlabel(ax,"PC1 (" +  sprintf("%4.2f",explained(1)) + "%)")
    ylabel(ax,"PC2 (" +  sprintf("%4.2f",explained(2)) + "%)")
    zlabel(ax,"PC3 (" +  sprintf("%4.2f",explained(3)) + "%)")
end
