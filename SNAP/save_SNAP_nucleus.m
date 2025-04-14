function save_SNAP_nucleus(path, data)
    try
        save(path,...
            "-struct",...
            "data")
    catch ME
        disp(getReport(ME));
        disp("Error saving: " + path)
    end
end

