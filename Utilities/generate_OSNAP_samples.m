function data = generate_OSNAP_samples(file_path,pixel_size)
    % load data
    locs_loaded = readtable(file_path,'ReadRowNames',0,'VariableNamingRule','preserve');
    % flip to match conventional images
    % y = max(y) - y; 
    % if pixel_size == 117 % pixel size in nm
    %     y = 684 - y;
    % end
    % clean data
    [~,data.name,~] = fileparts(file_path);
    data_to_unique(:,1) = locs_loaded{:,1};
    data_to_unique(:,2) = locs_loaded{:,2};
    data_to_unique = unique(data_to_unique,'rows');
    x = data_to_unique(:,1)*pixel_size;
    y = data_to_unique(:,2)*pixel_size;
    data.x = x;
    data.y = y;
    data.pixel_size_nm = pixel_size;
end

