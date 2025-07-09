function top_in_search_path_OSNAP(function_name,sub_string)
arguments
    function_name string
    sub_string string = []
end
    filepaths = which(function_name,'-all');
    if contains(filepaths{1},sub_string)
        return
    end
    for i=2:numel(filepaths)
        if contains(filepaths{i},sub_string)
            dir_top = fileparts(filepaths{i});
            break
        end
    end
    rmpath(dir_top);
    addpath(dir_top);
end