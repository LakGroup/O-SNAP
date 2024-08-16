function save_voronoi_data(filepath, data)
    if exist(filepath,'file')
        try
            save(filepath,...
                "-struct",...
                "data",...
                "-append")
        catch ME
                disp(getReport(ME));
                disp("Error saving: " + filepath)
        end
    else
        try
            save(filepath,...
                "-struct",...
                "data")
        catch ME
            disp(getReport(ME));
                disp("Error saving: " + filepath)
        end
    end
end
