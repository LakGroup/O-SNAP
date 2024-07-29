function percentileData = generatePercentileData(workDir, groups, voronoiData)
n_groups = length(groups);
% create list of percentiles
percentiles = (1:19)*5;
n_percentiles = size(percentiles,2);
percentileFields = arrayfun(@(elem) num2str(elem), percentiles, 'uni', 0);
percentileData = cell(1,n_groups);
[sampleNames, groupIdx] = sortSampleToGroups(groups,voronoiData);

%% calculate percentiles to y and the return variable
for g=1:n_groups
    percentileData{g}.group = groups{g};
    nSampInGroup = length(groupIdx{g});
    if nSampInGroup > 0
        percentileData{g}.sampleNames = sampleNames(groupIdx{g});
        for p=1:n_percentiles
            percentileData{g}.(['prctile_' percentileFields{p}]) = arrayfun(@(x) prctile(voronoiData{x}.reduced_log_voronoi_density,percentiles(p)), groupIdx{g});
        end
    end
end

save(fullfile(workDir,'percentileData'),'percentileData');
end