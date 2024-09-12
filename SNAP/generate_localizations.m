function data = generate_localizations(data,pixel_size)
    % calculate voronoi segmentation
    x = data.x_data;
    y = data.y_data;
    % flip to match conventional images
    % y = max(y) - y; 
    if pixel_size == 117
        y = 684 - y;
    end
    % clean data
    data_to_unique(:,1) = x;
    data_to_unique(:,2) = y;
    data_to_unique = unique(data_to_unique,'rows');
    x = data_to_unique(:,1)*pixel_size;
    y = data_to_unique(:,2)*pixel_size;
    data.x = x;
    data.y = y;
end

