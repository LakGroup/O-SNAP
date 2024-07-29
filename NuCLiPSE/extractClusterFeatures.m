function extractClusterFeatures(filepath)
    if ~hasVars(filepath,{'clusterNLocs','clusterArea','clusterDensity','clusterGyrationR'})
        varsToLoad = {'x','y','voronoi_areas_all',...
            'area_thresholds','min_number_of_localizations','clusters'};
        data = loadVars(filepath, varsToLoad);
        nThresholds = length(data.area_thresholds);
        data.clusterNLocs = cell(nThresholds,1);
        data.clusterArea = cell(nThresholds,1);
        data.clusterDensity = cell(nThresholds,1);
        data.clusterGyrationR = cell(nThresholds,1);
        for i=1:nThresholds
            clusters = data.clusters(:,i);
            nCluster = max(clusters);
            filtered_points = [data.x data.y data.voronoi_areas_all clusters];
            filtered_points = filtered_points(clusters > 0,:);
            data.clusterNLocs{i} = arrayfun(@(x) sum(filtered_points(:,4) == x), 1:nCluster,'uni',1)';
            data.clusterArea{i} = arrayfun(@(x) sum(filtered_points(filtered_points(:,4) == x,3),'omitmissing'), 1:nCluster,'uni',1)';
            data.clusterDensity{i} = arrayfun(@(x) data.clusterNLocs{i}(x)/data.clusterArea{i}(x), 1:nCluster,'uni',1)';
            data.clusterGyrationR{i} = arrayfun(@(x) sqrt(...
                mean((filtered_points(filtered_points(:,4) == x,1)...
                -mean(filtered_points(filtered_points(:,4) == x,1))).^2 ...
                + (filtered_points(filtered_points(:,4) == x,2) ...
                -mean(filtered_points(filtered_points(:,4) == x,2))).^2)), 1:nCluster,'uni',1);
        end
        clearvars -except filepath data varsToLoad
        data = rmfield(data,varsToLoad);
        saveVoronoiAnalysisData(filepath, data);
    end
end
