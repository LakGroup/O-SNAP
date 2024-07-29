function saveVoronoi(filePath, data)
    % name  = data.name;
    % x = data.x;
    % y = data.y;
    % voronoi_areas = data.voronoi_areas;
    % voronoi_areas_all = data.voronoi_areas_all;
    % voronoi_neighbors = data.voronoi_neighbors;
    % reduced_log_voronoi_density = data.reduced_log_voronoi_density;
    if exist(filePath,'file')
        try
            save(filePath,...
                "-struct",...
                "data",...
                "-append")
        catch ME
                disp(getReport(ME));
                disp("Error saving: " + filepath)
        end
    else
        save(filePath,...
            "-struct",...
            "data")
    end
end