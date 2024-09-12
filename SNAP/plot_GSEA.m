function plot_GSEA(work_dir,T_GSEA,feature_universe_name,options)
arguments
    work_dir string
    T_GSEA table
    feature_universe_name string
    options.alpha double = 0.05;
end
up_or_down = any([T_GSEA.ES > 0 T_GSEA.ES < 0]);
if all(up_or_down)
    figure('Position',[10 10 1210 910]);
    %% upregulation
    subplot(1,2,1)
    plot_GSEA_axis(T_GSEA,"dir","up","alpha",options.alpha);
    %% downregulation
    subplot(1,2,2)
    plot_GSEA_axis(T_GSEA,"dir","down","alpha",options.alpha);
elseif up_or_down(1)
    figure('Position',[10 10 610 910]);
    plot_GSEA_axis(T_GSEA,"dir","up","alpha",options.alpha);
else
    figure('Position',[10 10 610 910]);
    plot_GSEA_axis(T_GSEA,"dir","down","alpha",options.alpha);
end
sgtitle(replace(feature_universe_name,"_"," "))
%% save
savefig(fullfile(work_dir,"GSEA_plots","GSEA_plots_"+feature_universe_name +".fig"));
saveas(gcf,fullfile(work_dir,"GSEA_plots","GSEA_plots_"+feature_universe_name +".png"));
end

function plot_GSEA_axis(T_GSEA,options)
arguments
    T_GSEA table
    options.dir string = "up";
    options.alpha double = 0.05;
end
switch options.dir 
    case "up"
        T_GSEA_filt = T_GSEA(T_GSEA.ES > 0,:);
        T_GSEA_filt = sortrows(T_GSEA_filt,"NES","descend");
    case "down"
        T_GSEA_filt = T_GSEA(T_GSEA.ES < 0,:);
        T_GSEA_filt = sortrows(T_GSEA_filt,"NES","ascend");
    otherwise
        return
end
n_feature_families = size(T_GSEA_filt,1);
x = flip(1:n_feature_families);
y = abs(T_GSEA_filt.NES)';
[feature_family_size,center,scale] = normalize(T_GSEA_filt.("GS size"),'range');
s = feature_family_size'*200+50;
c = T_GSEA_filt.("NES_q-val");
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
yticklabels(flip(replace(T_GSEA_filt.("GS.NAME"),"_"," ")));
ytickangle(30)
xlabel("Normalized Enrichment Score")
ylim([0.5,n_feature_families+0.5])
title(upper(options.dir)+"-REGULATION")
cbar.Label.String = "adjusted p-value";
clim([0, max(c)]);
cbar.Direction = 'reverse';
cbar.TickLabels = compose("%0.1e",cbar.Ticks);
marker_refs = 5*unique(round(T_GSEA_filt.("GS size")/5),'sorted');
if marker_refs(1) == 0
    marker_refs(1) = 5;
end
if length(marker_refs) > 2
    marker_labels = [marker_refs(1) round(mean(marker_refs)) marker_refs(end)];
else
    marker_labels = marker_refs;
end
n_markers = length(marker_labels);
marker_sizes = (marker_labels - 5*round(center/5))./scale*200+50;
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
leg.Title.NodeChildren.Position = [0.5 1.2 0];
leg.Location = "southeast";
leg_height = leg.Position(4) * 1.3;
leg.Position(2) = leg.Position(2) - (leg_height - leg.Position(4)) + n_feature_families*0.005;
leg.Position(4) = leg_height;
<<<<<<< HEAD:NuCLiPSE/plot_GSEA.m
end
=======
end
>>>>>>> experimental:SNAP/plot_GSEA.m
