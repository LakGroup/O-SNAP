% Generate a nucleus
function points = generate_nucleus(axes_lengths,center,angle,int_periphery_cutoff,...
    density_interior,density_periphery,density_hetero_center,density_hetero,radius_hetero,...
    options)
    arguments
        axes_lengths (1,2) double
        center (1,2) double
        angle double
        int_periphery_cutoff double
        density_interior double
        density_periphery double
        density_hetero_center double
        density_hetero double
        radius_hetero double
        options.plot_flag logical = false
    end
    axes_lengths_interior = axes_lengths*sqrt(int_periphery_cutoff);
    points_interior = generate_ellipse(density_interior,center,axes_lengths_interior);
    points_periphery =  generate_ellipse(density_periphery,center,axes_lengths,axes_lengths_interior);
    points_hetero_center = generate_ellipse(density_hetero_center,center,axes_lengths_interior);
    points_hetero = [];
    for i=1:size(points_hetero_center,1)
        a = 360*rand(); % random angle
        l = random("Normal",radius_hetero,25,1,2); % random ellipse axes
        points_cluster = generate_ellipse(density_hetero,points_hetero_center(i,:),l);
        points_cluster = rotate_points(points_cluster, points_hetero_center(i,:), a);
        points_hetero = [points_hetero; points_cluster];
    end
    points_interior=rotate_points(points_interior,center,angle);
    points_periphery=rotate_points(points_periphery,center,angle);
    points_hetero_center=rotate_points(points_hetero_center,center,angle);
    if ~isempty(points_hetero)
        points_hetero=rotate_points(points_hetero,center,angle);
    end
    % plot
    if options.plot_flag && ~isempty(points_hetero)
        plot_nucleus(points_interior,points_periphery,points_hetero_center,points_hetero);
    elseif options.plot_flag
        plot_nucleus(points_interior,points_periphery,points_hetero_center);
    end
    points = [points_interior; points_periphery; points_hetero];
end

% Plot nucleus
function plot_nucleus(points_interior,points_periphery,points_hetero_center,points_hetero)
    plot_idx =  floor(1:length(points_interior)/100:length(points_interior));
    scatter(points_interior(plot_idx,1),points_interior(plot_idx,2),".k");
    hold on
     plot_idx =  floor(1:length(points_periphery)/100:length(points_periphery));
    scatter(points_periphery(plot_idx,1),points_periphery(plot_idx,2),"sr","filled");
    hold on
    plot_idx =  floor(1:length(points_hetero_center)/100:length(points_hetero_center));
    scatter(points_hetero_center(plot_idx,1),points_hetero_center(plot_idx,2),"xb");
    axis equal
    xlim([-15000 15000])
    ylim([-15000 15000])
    hold on
    if nargin == 4
        plot_idx =  floor(1:length(points_hetero)/100:length(points_hetero));
        scatter(points_hetero(plot_idx:end,1),points_hetero(plot_idx:end,2),0.5,".m")
        hold on
    end
end

% Generate points in an ellipse or ring
function points = generate_ellipse(density,center,axes_lengths1,axes_lengths2)
    points = generate_box_points(density,center,axes_lengths1);
    switch nargin
        case 4
            idx = in_ellipse(points,center,axes_lengths1,axes_lengths2);
        case 3
            idx = in_ellipse(points,center,axes_lengths1);
        otherwise
            return
    end 
    points = points(idx,:);
end

% Rotate points about a center
function points_rotated = rotate_points(points,center,angle)
    c = cosd(angle); s = sind(angle); x = center(1); y = center(2);
    rotate_matrix = [c -s -x*c+y*s+x; s c -x*s-y*c+y; 0 0 1];
    points_rotated = (rotate_matrix*[points'; ones(1,size(points,1))])';
    points_rotated = points_rotated(:,1:2);
end

% Get points randomly located within bounding box
function points = generate_box_points(density,center,axes_lengths)
    n_points = floor(density*prod(axes_lengths));
    x = axes_lengths(1) * rand(n_points,1) - axes_lengths(1)/2 + center(1);
    y = axes_lengths(2) * rand(n_points,1) - axes_lengths(2)/2 + center(2);
    points = [x y];
end

% Determine which points are in ellipse
function idx = in_ellipse(points, center, axes_lengths1, axes_lengths2)
    n_points = size(points,1);
    switch nargin
        case 4
            idx = all([((points(:,1) - repmat(center(:,1),n_points,1))/(axes_lengths2(1)/2)).^2 ...
            + ((points(:,2) - repmat(center(:,2),n_points,1))/(axes_lengths2(2)/2)).^2 >= 1 ...
            ((points(:,1) - repmat(center(:,1),n_points,1))/(axes_lengths1(1)/2)).^2 ...
            + ((points(:,2) - repmat(center(:,2),n_points,1))/(axes_lengths1(2)/2)).^2 <= 1
            ],2);
        case 3
            idx = ((points(:,1) - repmat(center(:,1),n_points,1))/(axes_lengths1(1)/2)).^2 ...
            + ((points(:,2) - repmat(center(:,2),n_points,1))/(axes_lengths1(2)/2)).^2 < 1;
        otherwise
            return
    end
end

