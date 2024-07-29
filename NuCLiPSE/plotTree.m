function [h,tt]=plotTree(ax, pos)
    tree = linkage(pos,'ward');
    D = pdist(pos);
    leafOrder = optimalleaforder(tree,D);
    [h,tt]=dendrogram(tree,'Reorder',leafOrder,'ColorThreshold','default','Parent',ax);
    title("Tree")
end