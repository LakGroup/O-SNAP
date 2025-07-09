function save_OSNAP_sample(path, data)
    try
        save(path,...
            "-struct",...
            "data")
    catch ME
        disp(getReport(ME));
        disp("Error saving: " + path)
    end
end

