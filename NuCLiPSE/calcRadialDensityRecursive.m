% WARNING: technically overestimates by including points beyond the outermost ring
% WARNING: uses a LOT of memory; can cause crashes
function radialDensity = calcRadialDensityRecursive(p, points)
    n_r = length(p);
    idxSplit = ceil(n_r/2);
    p_1 = union(p(1:idxSplit));
    in_p_1 = insidepoly(points(:,1),points(:,2),p_1.Vertices(:,1),p_1.Vertices(:,2)); % find points in first half
    if n_r == 1 % if there is only 1 ring, all the points that are passed will be inside that ring
        radialDensity = size(points,1)/area(p);% calculate density of given ring
    elseif n_r == 2 % if there are 2 rings, split points into either in ring1 or in ring2 (speeds up calculation)
        n_in_p_1 = sum(in_p_1); % find interior points at each r
        radialDensity = [n_in_p_1/area(p(1)); (size(points,1)-n_in_p_1)/area(p(2))];% calculate densities, where ring 2 includes all points not in ring 1
    else % halve the number of rings and recursively call function
        radialDensity = zeros(n_r,1);
        radialDensity(1:idxSplit) = calcRadialDensity(p(1:idxSplit), points(in_p_1,:)); % first half
        if idxSplit > 1
            radialDensity(idxSplit+1:end) = calcRadialDensity(p(idxSplit+1:end), points(~in_p_1,:)); % second half
        end
    end
end