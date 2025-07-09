% load specified variables from file
function data = load_variables_OSNAP(file_path, vars_to_load, options)
    arguments
        file_path string
        vars_to_load cell
        options.verbose logical = true
    end
    try
        [file_has_variables_OSNAP, vars_missing] = has_variables_OSNAP(file_path, vars_to_load, "verbose", options.verbose);
        if file_has_variables_OSNAP
            for i=1:length(vars_to_load)
                try
                    vars_value = load(file_path,vars_to_load{i});
                    data.(vars_to_load{i}) = vars_value.(vars_to_load{i});
                catch
                    ME = MException('SNAP:failed_to_load_existing_variable', ...
                        'Failed to load ''%s'' from %s',vars_to_load{i},file_path);
                    throw(ME)
                end
            end
        else
            ME = MException('SNAP:variable_not_in_file', ...
                        'File is missing ''%s''',string(join(vars_missing,''',''')));
            throw(ME)
        end
    catch ME
        throw(ME);
    end
end
