function generateVoronoiData(filepath, min_logVorDensity, max_logVorDensity)
    data = [];
    if ~hasVars(filepath,{'voronoi_areas_all','voronoi_neighbors'})
        varsToLoad = {'x','y'};
        data = loadVars(filepath, varsToLoad);
        % calculate voronoi segmentation
        x = data.x;
        y = data.y;
        % calculate voronoi cells and areas
        dt = delaunayTriangulation(x,y);
        [vertices,connections] = voronoiDiagram(dt);
        voronoi_cells = cellfun(@(x) vertices(x,:),connections,'uni',0);
        voronoi_areas_all = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'uni',0);
        voronoi_areas_all = vertcat(voronoi_areas_all{:});
        idx = ~isnan(voronoi_areas_all);
        voronoi_areas = voronoi_areas_all(idx);
        connections = connections(idx);
        reduced_voronoi_areas = voronoi_areas / mean(voronoi_areas);
        reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
        % find neighbors
        connectivity_list = dt.ConnectivityList;
        attached_triangles = vertexAttachments(dt);
        neighbors = cellfun(@(x) connectivity_list(x,:),attached_triangles,'uni',0);
        neighbors = cellfun(@(x) unique(x),neighbors,'uni',0);
        % return relevant info
        data.voronoi_areas_all = voronoi_areas_all;
        data.voronoi_areas = voronoi_areas;
        data.voronoi_neighbors = neighbors;
        data.reduced_log_voronoi_density = reduced_log_voronoi_density;
        % colorbar setup
        scale_bar_size_x = 2;
        scale_bar_size_y = 0.5;
        CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap');
        map = (colormap(CT));
        %% WITH REDUCED DENSITIES
        % clean data
        idx = reduced_log_voronoi_density>=min_logVorDensity;    
        connections = connections(idx);    
        reduced_log_voronoi_density_plot = reduced_log_voronoi_density(idx);  
        reduced_log_voronoi_density_plot(reduced_log_voronoi_density_plot>max_logVorDensity) = max_logVorDensity;
        max_connections = max(cellfun(@(x) length(x), connections));
        faces = NaN(length(connections),max_connections);
        for j = 1:length(connections)
            faces(j,1:length(connections{j})) = connections{j};
        end
        % plot voronoi segmentation for this sample
        f = figure('visible','off');
        set(f,'name','Log Voronoi Density Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.6 0.7]);
        patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',reduced_log_voronoi_density_plot,'FaceColor','flat','edgecolor','none')    
        h = colorbar();    
        pbaspect([1 1 1])
        h.Position = [0.9 0.1 0.02 0.8];
        h.TickLabelInterpreter = 'latex';
        h.FontSize = 14;
        h.Color = 'w';
        logVorDensityTicks = unique(sort([0 0.1 0.25, 0.41, 0.7, 1, 1.7, 2, min_logVorDensity, max_logVorDensity]));
        logVorDensityTicks = logVorDensityTicks(all([(logVorDensityTicks >= min_logVorDensity); (logVorDensityTicks <= max_logVorDensity)]));
        h.Ticks = logVorDensityTicks;
        h.TickLabels = arrayfun(@(elem) sprintf('%.2f',elem), logVorDensityTicks,  'uni',  0);
        ylabel(h, 'Reduced Log Voronoi polygon density $(\textrm{nm}^{-2})$','FontSize',14,'interpreter','latex');
        caxis([min_logVorDensity max_logVorDensity]);
        pixels_um_x = scale_bar_size_x*1000;
        pixels_um_y = scale_bar_size_y*1000;
        rectangle('Position',[min(x) min(y) pixels_um_x pixels_um_y],'facecolor','w')
        set(gca,'colormap',map,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
        axis equal
        xlim([min(x) max(x)])
        ylim([min(y) max(y)])    
        set(f, 'InvertHardCopy', 'off'); 
        [dirName,name,~] = fileparts(filepath);
        saveas(f,fullfile(replace(dirName,"VoronoiAreaData","ReducedVoronoiDensityMaps"),name+"_reduced.png"));
        close();
        clearvars -except data filepath
    elseif ~hasVars(filepath,{'voronoi_areas','reduced_log_voronoi_density'})
        varsToLoad = {'x','y','voronoi_areas_all','voronoi_neighbors'};
        data = loadVars(filepath, varsToLoad);
        idx = ~isnan(data.voronoi_areas_all);
        data.voronoi_areas = data.voronoi_areas_all(idx);
        reduced_voronoi_areas = data.voronoi_areas / mean(data.voronoi_areas);
        data.reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
        clearvars -except data filepath
    end
    saveVoronoiAnalysisData(filepath, data);
end