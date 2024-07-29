function batchPlotRadialDensity(filePath, maxRadialDensity, maxRadialHeteroDensity)
    varsToLoad = {'locs_norm','x_length','y_length',...
        'components','boundary'...
        'radialDensity','radialHeteroDensity',...
        'hetero_center','voronoi_areas_all'};
    data = loadVars(filePath,varsToLoad);
    [saveDir,name,~] = fileparts(filePath);
    saveDir = replace(saveDir,"VoronoiAreaData","ReducedVoronoiDensityMaps");
    %% plot
    plotFile = fullfile(saveDir,name+"_radialDensity.png");
    % if ~exist(plotFile,'file')
        plotRadialDensity(data.radialDensity, data.locs_norm, data.voronoi_areas_all, data.x_length, data.y_length, maxRadialDensity, 1, plotFile);
    % end
    
    plotFile = fullfile(saveDir,name+"_radialHeteroDensity.png");
    % if ~exist(plotFile,'file')
        plotRadialDensity(data.radialHeteroDensity,...
            normalizePoints(data.hetero_center,data.components,data.boundary),...
            zeros(size(data.hetero_center,1),1), data.x_length, data.y_length, maxRadialHeteroDensity, 10, plotFile);
    % end
end

function points_norm = normalizePoints(points,components,Boundary)
    rotatedBoundary = Boundary*components;
    points_norm = points*components - min(rotatedBoundary) - range(rotatedBoundary)/2;
end