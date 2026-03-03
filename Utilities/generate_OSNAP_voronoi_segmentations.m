% -------------------------------------------------------------------------
% generate_OSNAP_voronoi_segmentations.m
% -------------------------------------------------------------------------
% Performs Voronoi segmentation on a sample of interest and optionally
% plots the Voronoi density map. Saves information to the MAT analysis file
% specific to that sample.
%
% Example on how to use it:
%   generate_OSNAP_voronoi_segmentations(...
%                       "D:\Analysis\ExperimentA\Rep1\sample_KO_1.mat")
% -------------------------------------------------------------------------
% Input:
%   filepath: File path to specific sample (nucleus) to analyze.
% Output:
%   data: A struct array with information on the Voronoi data for a sample
%       - x: x-coordinates of localiations [nm]
%       - y: y-coordinates of localiations [nm]
%       - voronoi_areas_all: All initially calculated Voronoi polygon areas
%                           following Voronoi segmentation of localizations
%       - voronoi_areas: Voronoi polygon areas, filtered for NaN values and
%                        areas associated with the boundary localixations,
%                        which tend to demonstrate large outlier values
%       - voronoi_neighbors: A list of indices that store information on
%                            the relationship between the locs/Voronoi 
%                            cells that neighbor each other
%       - reduced_log_voronoi_density: The reduced log10 Voronoi density,
%                                      calculated from the filtered 
%                                      voronoi_areas
%       - faces: Information on how the vertices of the Voronoi cells
%                should connect to each other
%       - vertices: Information on the vertex values of the Voronoi cells
% Options:
%   min_log_vor_density: Lower bound for filtering and visualizing reduced
%                        Voronoi density values
%   max_log_vor_density: Upper bound for filtering and visualizing reduced
%                        Voronoi density values
%   min_log_vor_area: Lower bound for filtering and visualizing Voronoi
%                     area values
%   max_log_vor_area: Upper bound for filtering and visualizing Voronoi
%                     area values
%   plot: Flag to plot results
%   overwrite: Flag on whether to ask user before overwriting data
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
%%
function data = generate_OSNAP_voronoi_segmentations(filepath, options)
arguments
    filepath string
    options.min_log_vor_density double = 0;
    options.max_log_vor_density double = 3;
    options.min_log_vor_area double = 1;
    options.max_log_vor_area double = 3;
    options.plot logical = true;
    options.overwrite logical = false;
end
vor_data_vars = {...
    'x',...
    'y',...
    'voronoi_areas_all',...
    'voronoi_neighbors',...
    'voronoi_areas',...
    'reduced_log_voronoi_density',...
    'faces',...
    'vertices'};

% If overwrite is on, rerun the calculation
if options.overwrite
    data = calculate_voronoi_density(filepath);
    if options.plot
        plot_OSNAP_voronoi_map(filepath,data,options)
    end
% Return if all variables already in file
elseif has_variables_OSNAP(filepath,vor_data_vars) && ~options.overwrite
    % Load data and plot if plot flag
    if options.plot
        if ~exist('data','var')
            data = load_variables_OSNAP(filepath,vor_data_vars);
        end
        plot_OSNAP_voronoi_map(filepath,data,options)
    end
% If plot but no data is present in the file, calculate densities and plot
elseif options.plot
    data = calculate_voronoi_density(filepath);
    plot_OSNAP_voronoi_map(filepath,data,options)
end
end

%% Helper functions
% Calculate voronoi density values and save to file
function data = calculate_voronoi_density(filepath)
    data = load(filepath);
    % Calculate voronoi segmentation
    dt = delaunayTriangulation(data.x,data.y);
    [vertices,connections] = voronoiDiagram(dt);
    voronoi_areas_all = nan(numel(connections),1);
    for i = 1:numel(connections)
        tmp = vertices(connections{i},:);
        if ~any(isnan(tmp) | isinf(tmp))
            voronoi_areas_all(i) = polyarea(tmp(:,1),tmp(:,2));
        end
    end
    % Clean data (boundary points that cause outliers and NaN values)
    voronoi_areas = voronoi_areas_all;
    boundary_idx = boundary(data.x,data.y);
    voronoi_areas(boundary_idx) = [];
    voronoi_areas = voronoi_areas(~isnan(voronoi_areas));
    % connections = connections(idx);
    reduced_voronoi_areas = voronoi_areas ./ mean(voronoi_areas);
    reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
    % Find neighbors
    connectivity_list = dt.ConnectivityList;
    attached_triangles = vertexAttachments(dt);
    neighbors = cell(numel(attached_triangles),1);
    for i = 1:numel(attached_triangles)
        neighbors{i} = connectivity_list(attached_triangles{i},:);
        neighbors{i} = unique(neighbors{i});
        neighbors{i}(neighbors{i}==i) = [];
    end
    % Return relevant info
    data.voronoi_areas_all = voronoi_areas_all;
    data.voronoi_areas = voronoi_areas;
    data.voronoi_neighbors = uneven_cell2mat(neighbors);
    data.reduced_log_voronoi_density = reduced_log_voronoi_density;
    data.faces = uneven_cell2mat(connections)';
    data.vertices = vertices;
    save_OSNAP_sample(filepath, data);
end

% Change structure of neighbors from cell array to matrix to reduce memory
function arr_new = uneven_cell2mat(arr)
    arr_new = nan(max(cellfun('length',arr)),length(arr));
    for i=1:length(arr)
        value = arr{i};
        arr_new(1:numel(value),i) = value;
    end
end

