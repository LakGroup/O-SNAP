function [c_idx, cColor] = extract_clusters_from_tree(ax,h,tt)
    % get bottom branches in tree
    bottom_lines = h(arrayfun(@(x) ismember(0,x.YData),h,'uni',1));
    % color information
    bottom_linesColors = arrayfun(@(x) x.Color, bottom_lines,'uni',0);
    bottom_linesColors = vertcat(bottom_linesColors{:});
    % using x-position of bottom branches, find which leaf indices 
    % correspond to each bottom branch, removing branch ends that do not
    % reach the x-axis (i.e. are not integers)
    bottom_lines_x = arrayfun(@(x) unique(x.XData),bottom_lines,'uni',0);
    bottom_lines_x = cellfun(@(x) x(round(x(:)) == (x(:))), bottom_lines_x,'uni',0);
    x_tickLabels = str2num(ax.XTickLabel);
    bottom_lines_xLabels = cellfun(@(x) x_tickLabels(x)', bottom_lines_x,'uni',0);
    % find the color mapping of each bottom branch (cColor: color list of clusters)
    [cColor,~,cColor_idx] = unique(bottom_linesColors,'rows');
    % map the cluster index to each cell (c_idx: cluster that each belongs to)
    c_idx = zeros(size(tt,1),1);
    for i=1:length(bottom_lines_xLabels)
        tt_idx = any(tt==bottom_lines_xLabels{i},2);
        c_idx(tt_idx) = cColor_idx(i);
    end
end

