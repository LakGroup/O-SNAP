function [cIdx, cColor] = extractClustersFromTree(ax,h,tt)
    % get bottom branches in tree
    bottomLines = h(arrayfun(@(x) ismember(0,x.YData),h,'uni',1));
    % color information
    bottomLinesColors = arrayfun(@(x) x.Color, bottomLines,'uni',0);
    bottomLinesColors = vertcat(bottomLinesColors{:});
    % using x-position of bottom branches, find which leaf indices 
    % correspond to each bottom branch, removing branch ends that do not
    % reach the x-axis (i.e. are not integers)
    bottomLinesX = arrayfun(@(x) unique(x.XData),bottomLines,'uni',0);
    bottomLinesX = cellfun(@(x) x(round(x(:)) == (x(:))), bottomLinesX,'uni',0);
    xTickLabels = str2num(ax.XTickLabel);
    bottomLinesXLabels = cellfun(@(x) xTickLabels(x)', bottomLinesX,'uni',0);
    % find the color mapping of each bottom branch (cColor: color list of clusters)
    [cColor,~,cColorIdx] = unique(bottomLinesColors,'rows');
    % map the cluster index to each cell (cIdx: cluster that each belongs to)
    cIdx = zeros(size(tt,1),1);
    for i=1:length(bottomLinesXLabels)
        ttIdx = any(tt==bottomLinesXLabels{i},2);
        cIdx(ttIdx) = cColorIdx(i);
    end
end