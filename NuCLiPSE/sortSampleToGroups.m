function [sampleNames, groupIdx] = sortSampleToGroups(groups,voronoiData)
sampleNames = cellfun(@(x) x.name, voronoiData, 'uni',1);
sampleGroups = cellfun(@(x) x.group, voronoiData, 'uni',1);
groupIdx = cellfun(@(x) find(sampleGroups==x), groups,'uni',0);
end


