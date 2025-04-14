close all force
rng("default") % reproducibility

% 
root_dir = "F:\";
analysis_name = "ANALYSIS_AM_H3Variant_H3-1_H2B";
% analysis_name = "ANALYSIS_AM_Heterokaryon_H2B_CtrlhFb_CtrlmESC_HK48h";
% analysis_name = "ANALYSIS_YZ_hChon_H2B";
work_dir = fullfile(root_dir,analysis_name);
image_dir = fullfile(root_dir,"ANALYSIS_AM_H3Variant_H3-1_H2B");
analysis_path = fullfile(work_dir,analysis_name+".mat");

%% load data
if ~exist("SNAP_data","var")
    SNAP_data = load_variables(analysis_path,{...
        'groups',...
        'replicates',...
        'train_idxs',...
        'test_idxs',...
        'feature_comparisons'...
        'feature_data',...
        'vars_selected_batch',...
        'vars_select_result_batch',...
        'pca_result_batch_all',...
        'pca_result_batch_each',...
        'classifiers_batch_all',...
        'classification_summary_batch_all', ...
        'classifiers_batch_each',...
        'classification_summary_batch_each'});
end
n_groups = numel(SNAP_data.groups);
n_compare = numel(SNAP_data.feature_comparisons);
n_each_group = groupcounts(SNAP_data.feature_data.group);
n_batch = numel(SNAP_data.train_idxs);

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
g_rep_images.RowHeight = {30,16,'1x',100};
g_rep_images.ColumnWidth =  repmat({'1x'},1,n_groups);
title = uilabel(g_rep_images,'Text','REPRESENTATIVE IMAGES');
title.Layout.Column = [1 n_groups];
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
try
    SNAP_nucleus_file_list = get_valid_SNAP_nucleus_files(image_dir,SNAP_data.groups,SNAP_data.replicates,{'voronoi_areas'});
    for i=1:n_groups
        title = uilabel(g_rep_images,'Text',sprintf('%s (N=%.0f nuclei)',SNAP_data.groups{i},n_each_group(i)));
        title.HorizontalAlignment = 'center';
        title.FontSize = 14;
    end
    filepaths_by_groups = cellfun(@(x) replace(replace(SNAP_nucleus_file_list{SNAP_nucleus_file_list.group == x,'filepath'},"SNAP_nucleus_data","SNAP_nucleus_images"),".mat",".png"),SNAP_data.groups,'uni',0);
    for i=1:n_groups
        im = uiimage(g_rep_images);
        im.ImageSource = filepaths_by_groups{i}(floor(numel(filepaths_by_groups{i})/2));
    end
    for i=1:n_groups
        im = g_rep_images.Children(1+n_groups+i);
        g_path_label = uigridlayout(g_rep_images,[2,1]);
        g_path_label.RowHeight =  {'1x','1.3x'}; 
        path_label = uilabel(g_path_label);
        idx = floor(numel(filepaths_by_groups{i})/2);
        c = uicontrol(fig,'Style','text');
        path_label.Text = textwrap(c,{im.ImageSource},strlength(extractBefore(im.ImageSource,'SNAP_nucleus_images'))+20);
        path_label.HorizontalAlignment = 'center';
        path_label.FontSize = 12;
        im_slider =  uislider(g_path_label, ...
            "Limits",[1 numel(filepaths_by_groups{i})],...
            "Value",idx);
        im_slider.ValueChangedFcn = @(src,event) update_rep_image_callback(fig,im_slider,im,path_label,filepaths_by_groups{i});
    end
catch
end
%% VOLCANO PLOT
p_volcano = uipanel(g);
g_volcano = uigridlayout(p_volcano,[3,1]);
g_volcano.RowHeight = {30,20,'1x'};
title = uilabel(g_volcano,'Text','VOLCANO PLOT');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
feature_comparisons_list = cellfun(@(x) join(x.groups," vs "),SNAP_data.feature_comparisons,'uni',0);
feature_comparisons_list = [feature_comparisons_list{:}];
gt_volcano = uigridlayout(g_volcano,[1,2]);
gt_volcano.Layout.Row = 3;
dd_volcano = uidropdown(g_volcano,...
    "Items",feature_comparisons_list,...
    "ValueChangedFcn",@(src,event) plot_volcano_table(src,fig,gt_volcano,SNAP_data.feature_comparisons,SNAP_data.vars_selected_batch),...
    "FontWeight","bold",...
    "FontSize",14);
dd_volcano.Layout.Row = 2;
gt_volcano.ColumnWidth = {'2x','3x'};
plot_volcano_table(dd_volcano,fig,gt_volcano,SNAP_data.feature_comparisons,SNAP_data.vars_selected_batch);
%% FEATURE SELECTION
p_feature_select = uipanel(g);
g_feature_select = uigridlayout(p_feature_select,[1 2]);
g_feature_select.RowHeight = {30,'1x'};
g_feature_select.ColumnWidth =  {'1.2x','1x'};
title = uilabel(g_feature_select,'Text','FEATURE SELECTION');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
title.Layout.Row = 1;
title.Layout.Column = [1 2];
v = renamevars(SNAP_data.vars_select_result_batch, ...
    SNAP_data.vars_select_result_batch.Properties.VariableNames, ...
    SNAP_data.vars_select_result_batch.Properties.VariableNames+"_MRMR");
v = sortrows([SNAP_data.vars_selected_batch(:,[end, 1:end-1]) v(:,[end, 1:end-1])],"batch_sum_MRMR","descend");
writetable(v,fullfile(work_dir,'feature_select_summary.csv'),'WriteRowNames',true);
vars_selected_uitable = uitable(g_feature_select,...
    'Data',v{:,:},...
    'SelectionType','row',...
    'ColumnFormat',[repmat({'logical'},1,n_batch+1) repmat({'numeric'},1,n_batch+1)],...
    'ColumnName',v.Properties.VariableNames,...
    'ColumnSortable',true,...
    'RowName',v.Properties.RowNames,...
    'Multiselect','on');
p_violin=uipanel(g_feature_select,'AutoResizeChildren','off');
g_violin = uigridlayout(p_violin,[2,1]);
create_violin_panel(fig,g_violin,vars_selected_uitable,SNAP_data.feature_data,SNAP_data.groups)
vars_selected_uitable.SelectionChangedFcn = @(src,event) update_feature_selections(fig,g_violin,src,SNAP_data.feature_data);
vars_selected_uitable.Selection = 1:4;
% format cells green
idx = find(vars_selected_uitable.Data(:,1));
addStyle(vars_selected_uitable,...
    uistyle("BackgroundColor",[0 1 0]),...
    "cell",[idx ones(numel(idx),1)]);
idx = find(~vars_selected_uitable.Data(:,1));
addStyle(vars_selected_uitable,...
    uistyle("BackgroundColor",[0.5 0.5 0.5]),...
    "cell",[idx ones(numel(idx),1)]);
% format cells yellow
for b=1:n_batch+1
    values_MRMR = normalize(v{:,(n_batch+1)+b},'range',[0 1]);
    for i=1:size(v,1)
        addStyle(vars_selected_uitable,...
            uistyle("BackgroundColor",[1 1 0]+([0 0 1]*(1-values_MRMR(i)))),...
            "cell",[i (n_batch+1)+b]);
    end
end
%% PCA/CLASSIFICATION
g_PCA_classification = uigridlayout(g,[1,2]);
g_PCA_classification.ColumnWidth = {'1x','1x'};
%% PCA
p_PCA = uipanel(g_PCA_classification);
g_PCA = uigridlayout(p_PCA,[3,1]);
g_PCA.RowHeight = {30,40,'1x'};
title = uilabel(g_PCA,'Text','PRINCIPAL COMPONENT ANALYSIS');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
g_PCA_ax = axes(g_PCA);
g_PCA_ax.Layout.Row = 3;
g_PCA_ui_comps = uigridlayout(g_PCA,[1,3]);
g_PCA_ui_comps.ColumnWidth = {'3x','1x','2x'};
g_PCA_ui_comps.Layout.Row = 2;
dd_PCA = uidropdown(g_PCA_ui_comps,...
    "Items","BATCH "+(1:n_batch),...
    "FontWeight","bold",...
    "FontSize",14);
sb_PCA = uibutton(g_PCA_ui_comps,"state",...
    "Text","Hide Test Data",...
    "VerticalAlignment","Center",...
    "Value",1,...
    "FontSize",14);
bg_PCA = uibuttongroup(g_PCA_ui_comps);
b_PCA_all = uiradiobutton(bg_PCA,"Text","All","Position",[10 0 50 22]);
b_PCA_all.Value = true;
b_PCA_each = uiradiobutton(bg_PCA,"Text","Each","Position",[70 0 50 22]);
dd_PCA.ValueChangedFcn = @(src,event) update_PCA(g_PCA_ax,SNAP_data.feature_data,SNAP_data.train_idxs,SNAP_data.test_idxs,SNAP_data.pca_result_batch_each,SNAP_data.pca_result_batch_all,sb_PCA,src,bg_PCA);
sb_PCA.ValueChangedFcn = @(src,event) update_PCA_sb(g_PCA_ax,SNAP_data.feature_data,SNAP_data.train_idxs,SNAP_data.test_idxs,SNAP_data.pca_result_batch_each,SNAP_data.pca_result_batch_all,src,dd_PCA,bg_PCA);
bg_PCA.SelectionChangedFcn = @(src,event) update_PCA(g_PCA_ax,SNAP_data.feature_data,SNAP_data.train_idxs,SNAP_data.test_idxs,SNAP_data.pca_result_batch_each,SNAP_data.pca_result_batch_all,sb_PCA,dd_PCA,src);
batch_idx = 1;
train_data = preprocess_SNAP_table(SNAP_data.feature_data(SNAP_data.train_idxs{batch_idx},:),"remove_NaN",true,'keep_rep_sample_info',true);
test_data = preprocess_SNAP_table(SNAP_data.feature_data(SNAP_data.test_idxs{batch_idx},:),"remove_NaN",true,'keep_rep_sample_info',true);
plot_PCA(g_PCA_ax,...
    SNAP_data.pca_result_batch_all{batch_idx},...
    train_data,...
    "test_data",test_data);
%% CLASSIFICATION
p_classification = uipanel(g_PCA_classification);
g_classification = uigridlayout(p_classification,[4, 1]);
g_classification.RowHeight =  {30,42,'1x','2x'};
title = uilabel(g_classification,'Text','CLASSIFICATION');
title.HorizontalAlignment = 'center';
title.FontSize = 24;
title.FontWeight = 'bold';
g_classification_ui_comps = uigridlayout(g_classification,[1, 2]);
l_classification = uilabel(g_classification_ui_comps,'Text',sprintf('Validation accuracies (N=%.0f batches)',numel(SNAP_data.train_idxs)));
l_classification.HorizontalAlignment = 'center';
l_classification.FontSize = 12;
bg_classification = uibuttongroup(g_classification_ui_comps);
b_classification_all = uiradiobutton(bg_classification,"Text","All","Position",[30 0 50 22]);
b_classification_each = uiradiobutton(bg_classification,"Text","Each","Position",[110 0 50 22]);
b_classification_all.Value = true;
t_classification_all = SNAP_data.classification_summary_batch_all(:,["mean_batch_acc" "std_batch_acc"]);
t_classification_all.Properties.RowNames = SNAP_data.classification_summary_batch_all.model_type;
t_classification_all = renamevars(t_classification_all,...
    t_classification_all.Properties.VariableNames,...
    t_classification_all.Properties.VariableNames+"_all");
t_classification_each = SNAP_data.classification_summary_batch_each(:,["mean_batch_acc" "std_batch_acc"]);
t_classification_each.Properties.RowNames = SNAP_data.classification_summary_batch_each.model_type;
t_classification_each = renamevars(t_classification_each,...
    t_classification_each.Properties.VariableNames,...
    t_classification_each.Properties.VariableNames+"_each");
t_classification = innerjoin(t_classification_all,t_classification_each,"Key","Row");
t_classification.group_count = SNAP_data.classification_summary_batch_all.n_batch;
t_classification = sortrows(t_classification,"mean_batch_acc_all","descend");
t_classification{:,1:4} = t_classification{:,1:4}*100;
writetable(t_classification,fullfile(work_dir,'classification_summary.csv'),'WriteRowNames',true);
classification_table = uitable(g_classification,...
    'Data',t_classification{:,:},...
    'RowName',t_classification.Properties.RowNames,...
    'ColumnSortable',true,...
    'ColumnName',t_classification.Properties.VariableNames,...
    'SelectionType','row');
classification_table.Selection = 1;
ax_confusion = axes(g_classification);
plot_SNAP_confusion(ax_confusion,SNAP_data.classification_summary_batch_all.model_type(1),SNAP_data.classifiers_batch_all) % plot top classifiers
classification_table.SelectionChangedFcn = @(src,event) update_classification_selection(ax_confusion,src,SNAP_data.classifiers_batch_all,SNAP_data.classifiers_batch_each,bg_classification);
bg_classification.SelectionChangedFcn = @(src,event)  update_classification_selection(ax_confusion,classification_table,SNAP_data.classifiers_batch_all,SNAP_data.classifiers_batch_each,src);
%% SAVE
% exportgraphics(fig,fullfile(work_dir,analysis_name+"_summary.pdf"));
% disp("Results saved to: "+analysis_name+"_summary.pdf")
%% Previous image in file list
function update_rep_image_callback(fig,im_slider,ax_image,ax_label,filepaths)
    new_idx = floor(im_slider.Value);
    c = uicontrol(fig,'Style','text');
    ax_image.ImageSource = filepaths(new_idx);
    ax_label.Text = textwrap(c,filepaths(new_idx),strlength(extractBefore(filepaths(new_idx),'SNAP_nucleus_images'))+20);           
end
function plot_volcano_table(src,fig,gridlayout,feature_comparisons,vars_selected)
    idx = src.ValueIndex;
    comparison_table = feature_comparisons{idx}.feature_table(:,2:end);
    comparison_table.Properties.RowNames = feature_comparisons{idx}.feature_table{:,1};
    comparison_table = comparison_table(:,...
        ["log2_fold_change","adj_p_value","fold_change","p_value" comparison_table.Properties.VariableNames(1:end-4)]);
    comparison_table = innerjoin(comparison_table,vars_selected,"Key","Row");
    comparison_table = sortrows(comparison_table,"log2_fold_change","descend");
    comparison_table = comparison_table(:,[end 1:end-1]);
    comparison_table = renamevars(comparison_table, "batch_sum","selected_all");
    if numel(gridlayout.Children) < 2
        g_volcano_axes = axes(gridlayout);
        g_volcano_table = uitable(gridlayout);
    else
        g_volcano_axes = gridlayout.Children(1);
        g_volcano_table = gridlayout.Children(2);
    end
    plot_volcano(fig,g_volcano_axes,feature_comparisons{idx});
    g_volcano_table.Data=comparison_table{:,:};
    g_volcano_table.ColumnFormat = ['logical' repmat({'numeric'},1,size(g_volcano_table.Data,2)-1)];
    g_volcano_table.ColumnName=comparison_table.Properties.VariableNames;
    g_volcano_table.RowName=comparison_table.Properties.RowNames;
    g_volcano_table.ColumnSortable = true;
    % highlight red
    highlight_idx = all([(g_volcano_table.Data(:,2)>1) (g_volcano_table.Data(:,3) < 0.05)],2);
    addStyle(g_volcano_table, ...
        uistyle("BackgroundColor",'#FF6262'), ...
        "row",find(highlight_idx));
    % highlight blue
    highlight_idx = all([(g_volcano_table.Data(:,2)<-1) (g_volcano_table.Data(:,3) < 0.05)],2);
    addStyle(g_volcano_table,...
        uistyle("BackgroundColor",'#62C8FF'),...
        "row",find(highlight_idx));
    %highlight green
    idx = find(g_volcano_table.Data(:,1));
    addStyle(g_volcano_table,...
        uistyle("BackgroundColor",[0 1 0]),...
        "cell",[idx ones(numel(idx),1)]);
    idx = find(~g_volcano_table.Data(:,1));
    addStyle(g_volcano_table,...
        uistyle("BackgroundColor",[0.5 0.5 0.5]),...
        "cell",[idx ones(numel(idx),1)]);
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
                @(ax, cpp)set(findobj(fig, 'tag',num2str(f)), ...
                'string', {feature_name(f), ""});
            pointerBehavior.traverseFcn = [];
            pointerBehavior.exitFcn = ...
                @(ax, cpp)set(findobj(fig,'tag',num2str(f)), ...
                'string', '    ');
            iptSetPointerBehavior(ht, pointerBehavior);
            iptPointerManager(fig, 'enable');
            hold(ax,'on')
        end
    end
    hold(ax,'off')
end
function create_violin_panel(fig,gridlayout,ui_table,feature_data,groups)
    ax_violin_legend = axes(gridlayout);
    n_legend_rows = create_violin_legend(ax_violin_legend,groups);
    gridlayout.RowHeight = {30*n_legend_rows,'1x'};
    p_violin_axes = uipanel(gridlayout,"BorderWidth",0);
    p_violin_axes.AutoResizeChildren = 'off';
    ax_violin_plot = tiledlayout(p_violin_axes,2,2,'TileSpacing','Compact','Padding','Compact');
    % plot new axes
    for i=1:4
        ax_violin_plot_i = nexttile(ax_violin_plot);
        plot_violin(fig,ax_violin_plot_i,feature_data,ui_table.RowName{i});
    end
end
%% create one legend for violin plots
function n_legend_rows = create_violin_legend(ax,groups)
    set(ax,'ydir','reverse');
    n_groups = numel(groups);
    c_groups = lines(numel(groups));
    x0 = 28;
    text_ratio = 1.8;
    x_center = x0/2;
    y_center = 1/2;
    icon_length = 0.6;
    offset_x = x0/50;
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
    % drawnow()
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
function update_feature_selections(fig,panel,ui_table,feature_data)
    idx_select = ui_table.Selection;
    if numel(idx_select) > 4
        ui_table.Selection = idx_select(1:4);
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
        plot_violin(fig,ax_violin_plot_i,feature_data,ui_table.RowName{idx_select(i)});
    end
end
%% generate violin plot
function plot_violin(fig,ax,T,feature)
arguments
    fig
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
    % mean -OUTLIERS ARE INCLUDED WHEN CALCULATING MEAN
    scatter(ax,g,mean(feature_data_g,"omitmissing"),50,'pentagramk','filled')
    hold(ax,"on")
    % outliers
    if sum(outliers_idx > 0)
        scatter(ax,repmat(g,1,sum(outliers_idx)),feature_data_g(outliers_idx),"xk");
    end
    ax.XTick = [];
end
scatter(ax,nan,nan,50,'pentagramk','filled');
plot(ax,nan, nan,"-k","LineWidth",2);
c = uicontrol(fig,'Style','text');
title(ax,textwrap(c,{replace(feature,"_"," ")},30))
end
%% update confusion matrix
function update_classification_selection(ax,classification_table,classifiers_batch_all,classifiers_batch_each,bg_classification)
    model_type = classification_table.RowName{classification_table.Selection};
    if strcmp(bg_classification.SelectedObject.Text,"Each")
        classifiers = classifiers_batch_each;
    else
        classifiers = classifiers_batch_all;
    end
    plot_SNAP_confusion(ax,model_type,classifiers)
end
%% update PCA state button
function update_PCA_sb(ax,feature_data,train_idxs,test_idxs,pca_result_batch_each,pca_result_batch_all,sb,dd,bg)
    update_PCA(ax,feature_data,train_idxs,test_idxs,pca_result_batch_each,pca_result_batch_all,sb,dd,bg)
    disp(sb.Value)
    if sb.Value
        sb.Text = "Hide Test";
        sb.Value = 1;
    else
        sb.Text = "Show Test";
        sb.Value = 0;
    end
end

%% update PCA
function update_PCA(ax,feature_data,train_idxs,test_idxs,pca_result_batch_each,pca_result_batch_all,sb,dd,bg)
    batch_idx = dd.ValueIndex;
    if sb.Value
        test_data = preprocess_SNAP_table(feature_data(test_idxs{batch_idx},:),"remove_NaN",true,'keep_rep_sample_info',true);
    else
        test_data = [];
    end
    if strcmp(bg.SelectedObject.Text,"Each")
        pca_result = pca_result_batch_each;
    else
        pca_result = pca_result_batch_all;
    end
    train_data = preprocess_SNAP_table(feature_data(train_idxs{batch_idx},:),"remove_NaN",true,'keep_rep_sample_info',true);
    plot_PCA(ax,...
        pca_result{batch_idx},...
        train_data,...
        "test_data",test_data);
end
%% Plot PCA
function plot_PCA(ax,pca_result,train_data,options)
arguments
    ax
    pca_result struct
    train_data table
    options.test_data table = []
end
    marker_list = {'o','^','square','v','diamond','pentagram','hexagram'};
    marker_size = 20;
    train_PCA = pca_result.pca_transformation_fcn(train_data(:,vartype('numeric')));
    %% plot
    [groups,~,groups_idx] = unique(train_data.group);
    [bio_reps,~,bio_reps_idx] = unique(train_data.biological_replicate);
    s_face_colors = lines(length(groups));
    for b=1:length(bio_reps)
        for g=1:length(groups)
            idx = all([(groups_idx == g) (bio_reps_idx == b)],2);
            if size(train_PCA,2) > 2
                s = scatter3(ax,train_PCA(idx,1),train_PCA(idx,2),train_PCA(idx,3),'filled',...
                    "MarkerFaceColor",s_face_colors(g,:),...
                    "Marker",marker_list(b));
            elseif size(train_PCA,2) == 2
                s = scatter(ax,train_PCA(idx,1),train_PCA(idx,2),'filled',...
                    "MarkerFaceColor",s_face_colors(g,:),...
                    "Marker",marker_list(b));
            elseif size(train_PCA,2) == 1
                hold off
                s = swarmchart(ax,ones(sum(idx),1),train_PCA(idx,1),'filled',...
                    'YJitter','density',...
                    "MarkerFaceColor",s_face_colors(g,:),...
                    "Marker",marker_list(b));
            else
                return
            end
            if b == 1
                s.SizeData = marker_size;
            else
                s.SizeData = 2*marker_size;
            end
            hold(ax,'on')
        end
    end
    if ~isempty(options.test_data)
        [~,~,groups_idx] = unique(options.test_data.group);
        test_PCA = pca_result.pca_transformation_fcn(options.test_data(:,vartype('numeric')));
        for g=1:length(groups)
            idx = (groups_idx == g);
            if ~isempty(idx)
                if size(test_PCA,2) > 2
                    scatter3(ax,test_PCA(idx,1),test_PCA(idx,2),test_PCA(idx,3),...
                        "MarkerEdgeColor", s_face_colors(g,:),...
                        "Marker","o","SizeData",marker_size,"LineWidth",1);
                elseif size(test_PCA,2) == 2
                    scatter(ax,test_PCA(idx,1),test_PCA(idx,2),...
                        "MarkerEdgeColor", s_face_colors(g,:),...
                        "Marker","o","SizeData",marker_size,"LineWidth",1);
                elseif size(test_PCA,2) == 1
                    swarmchart(ax,ones(sum(idx),1),test_PCA(idx,1),...
                        'YJitter','density',...
                        "MarkerEdgeColor",s_face_colors(g,:),...
                        "Marker","o","SizeData",marker_size,"LineWidth",1);
                else
                    return
                end
                hold(ax,'on')
            end
        end
    end
    legend(ax,groups,'Location','southoutside','Interpreter','none')
    % biplot(ax,coeff,'VarLabels',T_num.Properties.VariableNames);
    grid(ax,'on');
    set(ax, "XTickLabel",[]);
    xlabel(ax,"PC1 (" +  sprintf("%4.2f",pca_result.explained(1)) + "%)")
    set(ax,'XLimitMethod','tight')
    if size(train_PCA,2) > 1
        set(ax, "YTickLabel",[]);
        ylabel(ax,"PC2 (" +  sprintf("%4.2f",pca_result.explained(2)) + "%)")
        set(ax,'YLimitMethod','tight')
    end
    if size(train_PCA,2) > 2
        set(ax, "ZTickLabel",[]);
        zlabel(ax,"PC3 (" +  sprintf("%4.2f",pca_result.explained(3)) + "%)")
        set(ax,'ZLimitMethod','tight')
    end
    % hold(ax,"off");
    set(ax,'nextplot','replacechildren');
end
