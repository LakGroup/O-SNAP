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
%   min_log_reduced_vor_density: Lower bound for filtering and visualizing reduced
%                        Voronoi density values
%   max_log_reduced_vor_density: Upper bound for filtering and visualizing reduced
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
function data = generate_OSNAP_features_sample_voronoi_segmentations(filepath, options)
arguments
    filepath string
    options.min_log_reduced_vor_density double = 0;
    options.max_log_reduced_vor_density double = 3;
    options.min_log_vor_area double = 1;
    options.max_log_vor_area double = 3;
    options.plot logical = true;
    options.overwrite logical = false;
end
vor_data_vars = {...
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
        plot_OSNAP_voronoi_map(filepath,data,...
                        "min_log_reduced_vor_density", options.min_log_reduced_vor_density, ...
                        "max_log_reduced_vor_density", options.max_log_reduced_vor_density,...
                        "min_log_vor_area", options.min_log_vor_area, ...
                        "max_log_vor_area", options.max_log_vor_area);
    end
% Return if all variables already in file
elseif has_variables_OSNAP(filepath,vor_data_vars)
    if ~options.overwrite
        [dir_name,name,~] = fileparts(filepath);
        voronoi_density_map_dir = replace(dir_name,"_nucleus_data","_voronoi_density_maps");
        reduced_voronoi_density_map_dir = replace(dir_name,"_nucleus_data","_reduced_voronoi_density_maps");
        % Check plot is not present
        if ~exist(fullfile(voronoi_density_map_dir,name+".png"), 'file') || ~exist(fullfile(reduced_voronoi_density_map_dir,name+"_reduced.png"), 'file')
            % Load data and plot if plot flag
            if options.plot
                if ~exist('data','var')
                    data = load_variables_OSNAP(filepath,vor_data_vars);
                end
                plot_OSNAP_voronoi_map(filepath,data,...
                    "min_log_reduced_vor_density", options.min_log_reduced_vor_density, ...
                    "max_log_reduced_vor_density", options.max_log_reduced_vor_density,...
                    "min_log_vor_area", options.min_log_vor_area, ...
                    "max_log_vor_area", options.max_log_vor_area);
            end
            % Otherwise continue
        end
    end
% If plot but no data is present in the file, calculate densities and plot
elseif options.plot
    data = calculate_voronoi_density(filepath);
    plot_OSNAP_voronoi_map(filepath,data,...
                        "min_log_reduced_vor_density", options.min_log_reduced_vor_density, ...
                        "max_log_reduced_vor_density", options.max_log_reduced_vor_density,...
                        "min_log_vor_area", options.min_log_vor_area, ...
                        "max_log_vor_area", options.max_log_vor_area);
end
end

%% Helper functions
% Calculate voronoi density values and save to file
function data = calculate_voronoi_density(filepath)
    data = load(filepath);
    % Calculate voronoi segmentation
    dt = delaunayTriangulation(data.x,data.y);
    [vertices,connections] = voronoiDiagram(dt);
    % Find Voronoi cells
    voronoi_cells = cellfun(@(x) vertices(x,:),connections,'uni',0);
    % Calculate Voronoi areas
    voronoi_areas_all = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'uni',0);
    idx = cell2mat(cellfun(@(x) isnan(x) | x == Inf,voronoi_areas_all,'uni',0));
    voronoi_areas_all(idx) = {nan};
    voronoi_areas_all = cell2mat(voronoi_areas_all);
    % Clean data (boundary points that cause outliers and NaN values)
    voronoi_areas = voronoi_areas_all;
    boundary_idx = boundary(data.x,data.y);
    voronoi_areas(boundary_idx) = [];
    voronoi_areas = voronoi_areas(~isnan(voronoi_areas));
    % connections = connections(idx);
    reduced_voronoi_areas = voronoi_areas ./ mean(voronoi_areas);
    reduced_log_voronoi_density = log10(1./reduced_voronoi_areas);
    % Find Voronoi neighbors
    connectivity_list = dt.ConnectivityList;
    attached_triangles = vertexAttachments(dt);
    neighbors = cellfun(@(x) connectivity_list(x,:),attached_triangles,'UniformOutput',false);
    neighbors = cellfun(@(x) unique(x),neighbors,'uniformoutput',false);
    for i = 1:length(neighbors)
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
function data_arr = uneven_cell2mat(data_cell)
arguments
    data_cell cell
end
    data_arr = nan(max(cellfun('length',data_cell)),length(data_cell));
    for i=1:length(data_cell)
        value = data_cell{i};
        data_arr(1:numel(value),i) = value;
    end
end

