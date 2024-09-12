% WARNING: technically overestimates by including points beyond the outermost ring
% WARNING: uses a LOT of memory; can cause crashes
function radial_density = calculate_radial_density_recursive(p, points)
    n_r = length(p);
    idx_split = ceil(n_r/2);
    p_1 = union(p(1:idx_split));
    in_p_1 = insidepoly(points(:,1),points(:,2),p_1.Vertices(:,1),p_1.Vertices(:,2)); % find points in first half
    if n_r == 1 % if there is only 1 ring, all the points that are passed will be inside that ring
        radial_density = size(points,1)/area(p);% calculate density of given ring
    elseif n_r == 2 % if there are 2 rings, split points into either in ring1 or in ring2 (speeds up calculation)
        n_in_p_1 = sum(in_p_1); % find interior points at each r
        radial_density = [n_in_p_1/area(p(1)); (size(points,1)-n_in_p_1)/area(p(2))];% calculate densities, where ring 2 includes all points not in ring 1
    else % halve the number of rings and recursively call function
        radial_density = zeros(n_r,1);
        radial_density(1:idx_split) = calculate_radial_density(p(1:idx_split), points(in_p_1,:)); % first half
        if idx_split > 1
            radial_density(idx_split+1:end) = calculate_radial_density(p(idx_split+1:end), points(~in_p_1,:)); % second half
        end
    end
end

