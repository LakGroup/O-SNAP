% Fast way to check file for variables
function [result, vars_missing] = has_variables(file_path, vars_to_check, options)
    arguments
        file_path string
        vars_to_check cell
        options.verbose logical = false
    end
    vars_in_file = who("-file",file_path);
    vars_missing = setdiff(vars_to_check,vars_in_file);
    if ~isempty(vars_missing)
        result = 0;
        if options.verbose
            fprintf('''%s'' not in: %s\n',string(join(vars_missing, ''', ''')),file_path)
        end
    else
        result = 1;
    end
end
