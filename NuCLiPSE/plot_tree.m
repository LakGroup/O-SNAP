function [h,tt]=plot_tree(ax, pos)
    tree = linkage(pos,'ward');
    D = pdist(pos);
    leaf_order = optimalleaforder(tree,D);
    [h,tt]=dendrogram(tree,'Reorder',leaf_order,'ColorThreshold','default','Parent',ax);
    title("Tree")
end
