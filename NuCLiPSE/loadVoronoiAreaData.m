function voronoiData = loadVoronoiAreaData(T,varsToLoad)
    n_samps = size(T,1);
    if ~(n_samps > 0)
        error('No valid files found; Look in another folder or calculate Voronoi segmentations de novo')
    end
    %% load data
    w = waitbar(0, ['Loading file: 0/' num2str(n_samps)]);
    voronoiData = cell(1,n_samps);
    for s=1:n_samps
        waitbar(s/n_samps, w, ['Loading file: ' num2str(s) '/' num2str(n_samps)]);
        voronoiData{s} = load(fullfile(T{s,"filepath"}),varsToLoad{:});
        voronoiData{s}.BioRep = T{s,"rep"};
        voronoiData{s}.Group = T{s,"group"};
    end
    close(w)
end