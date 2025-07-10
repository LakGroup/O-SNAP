% -------------------------------------------------------------------------
% plot_OSNAP_FSEA.m
% -------------------------------------------------------------------------
% Produces a summary plot for the FSEA of a feature universe, displaying
% NES, adjusted q-values, and feature set sizes
%
% Example on how to use it:
%   plot_OSNAP_FSEA(T_FSEA,"Control","KO",save_path,options)
% -------------------------------------------------------------------------
% Input:
%   T_FSEA: The output from FSEA with information on feature family ES,
%           NES, and statistics
%   group_ctrl: String identifier for control phenotype
%   group_case: String identifier for case phenotype
%   save_path: Name of file location to save to, excluding extension
% Options:
%   alpha: Significance testing threshold
%   p_show: Upper bound for feature significance for a feature set to be
%           displayed
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
function plot_OSNAP_FSEA(T_FSEA,group_ctrl,group_case,save_path,options)
arguments
    T_FSEA table
    group_ctrl string
    group_case string
    save_path string
    options.alpha double = 0.05;
    options.q_show double =  0.5
end
up_flag = any(all([T_FSEA.ES > 0, T_FSEA.("NES_q-val") < options.q_show],2));
down_flag = any(all([T_FSEA.ES < 0, T_FSEA.("NES_q-val") < options.q_show],2));
if up_flag && down_flag
    figure('Position',[10 10 1210 300],'visible','off');
    %% feature value increase
    subplot(1,2,1)
    plot_OSNAP_FSEA_axis(T_FSEA,group_case,"case_greater_than_ctrl",true,"alpha",options.alpha);
    %% feature value decrease
    subplot(1,2,2)
    plot_OSNAP_FSEA_axis(T_FSEA,group_ctrl,"case_greater_than_ctrl",false,"alpha",options.alpha);
elseif up_flag
    figure('Position',[10 10 610 300],'visible','off');
    plot_OSNAP_FSEA_axis(T_FSEA,group_case,"case_greater_than_ctrl",true,"alpha",options.alpha);
elseif down_flag
    figure('Position',[10 10 610 300],'visible','off');
    plot_OSNAP_FSEA_axis(T_FSEA,group_ctrl,"case_greater_than_ctrl",false,"alpha",options.alpha);
else
    return
end
% sgtitle(replace(feature_universe_name,"_"," "))
%% save
savefig(save_path +".fig");
saveas(gcf,save_path +".png");
end

function plot_OSNAP_FSEA_axis(T_FSEA,group,options)
arguments
    T_FSEA table
    group string = "";
    options.case_greater_than_ctrl logical = true;
    options.alpha double = 0.05;
    options.p_show double =  0.5;
end
if options.case_greater_than_ctrl 
    T_FSEA_filt = T_FSEA(all([T_FSEA.ES > 0, T_FSEA.("NES_q-val") < options.p_show],2),:);
    T_FSEA_filt = sortrows(T_FSEA_filt,"NES","descend");
else
    T_FSEA_filt = T_FSEA(all([T_FSEA.ES < 0, T_FSEA.("NES_q-val") < options.p_show],2),:);
    T_FSEA_filt = sortrows(T_FSEA_filt,"NES","ascend");
end
n_feature_families = size(T_FSEA_filt,1);
x = flip(1:n_feature_families);
y = abs(T_FSEA_filt.NES)';
s = T_FSEA_filt{:,"FS size"}'*10+10;
c = T_FSEA_filt.("NES_q-val");
colormap(flipud(jet))
for i=1:n_feature_families
    l = plot([0,y(i)],[x(i) x(i)],"--k");
    if c(i) < options.alpha
        set(l,"LineWidth", 2);
    end
    hold on;
end
scatter(y,x,s,c,'filled');
hold on
cbar = colorbar();
yticks(1:n_feature_families);
yticklabels(flip(replace(T_FSEA_filt.("FS.NAME"),"_"," ")));
ytickangle(30)
xlabel("Normalized Enrichment Score")
ylim([-0.5,n_feature_families+0.25])
title(group,"Interpreter","none")
cbar.Label.String = "FDR q-value";
clim([0 options.p_show]);
cbar.Direction = 'reverse';
% cbar.TickLabels = compose("%0.2f",cbar.Ticks);
cbar.Ruler.TickLabelFormat = '%.2f';
marker_refs = 5*unique(round(T_FSEA.("FS size")/5),'sorted');
if marker_refs(1) == 0
    marker_refs(1) = 5;
end
if length(marker_refs) > 2
    marker_labels = [marker_refs(1) round(mean(marker_refs)) marker_refs(end)];
else
    marker_labels = marker_refs;
end
n_markers = length(marker_labels);
marker_sizes = marker_labels*10+10;
legentry=repelem({''},1,n_feature_families+1+n_markers);
for i = 1:n_markers
   scatter(nan,nan,marker_sizes,'k','filled');
   legentry{n_feature_families+1+i} = num2str(marker_labels(i));
end
[leg,att] = legend(legentry);
for i = 1:n_markers
   set(att(n_markers+i).Children,'MarkerSize',sqrt(marker_sizes(i)));
end
leg.Title.Visible = 'on';
leg.Title.String = "# features";
leg.Title.NodeChildren.Position = [0.5 1.1 0];
leg.Location = "southeast";
leg.Position(4) = leg.Position(4) * 1.4;
leg.Position(2) = leg.Position(2) + n_feature_families*0.005;

end

