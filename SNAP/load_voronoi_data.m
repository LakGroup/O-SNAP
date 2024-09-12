function voronoi_data = load_voronoi_data(T,vars_to_load)
    n_samples = size(T,1);
    if ~(n_samples > 0)
        error('No valid files found; Look in another folder or calculate Voronoi segmentations de novo')
    end
    %% load data
    w = waitbar(0, ['Loading file: 0/' num2str(n_samples)]);
    voronoi_data = cell(1,n_samples);
    for s=1:n_samples
        waitbar(s/n_samples, w, ['Loading file: ' num2str(s) '/' num2str(n_samples)]);
        voronoi_data{s} = load(fullfile(T{s,"filepath"}),vars_to_load{:});
        voronoi_data{s}.biological_replicate = T{s,"rep"};
        voronoi_data{s}.group = T{s,"group"};
    end
    close(w)
end

