% -------------------------------------------------------------------------
% plot_OSNAP_DBSCAN_LAD.m
% -------------------------------------------------------------------------
% Calculates features related to morphometrics, DBSCAN heterochromatin-
% associated localization clustering, and radial features for a given
% sample (nucleus). Saves information to the MAT analysis file specific to 
% that sample.
%
% Example on how to use it:
%   extract_OSNAP_morphometrics_dbscan_cluster_radial_features(...
%                   "D:\Analysis\ExperimentA\Rep1\sample_KO_1.mat",...
%                   70)
% -------------------------------------------------------------------------
% Input:
%   filepath: File path to specific sample (nucleus) to analyze.
% Options:
%   data: A struct array with information on the Voronoi data for a sample.
%         If left empty, the code will attempt to load relevant information
%         from the provided file.
%       - boundary: Boundary points of the nucleus
%       - components: PCA components on localization data to normalize
%       - polygon: Polyshape representation of the nucleus boundary
%       - locs_dbscan_cluster_periphery_labels: [Nx3] array containing 
%                                               information on the compact 
%                                               domains at the nucleus 
%                                               periphery. The coordinates 
%                                               of the locs in the first 
%                                               two columns and the index 
%                                               of the compact domain they 
%                                               belong to in the third 
%                                               column
%       - locs_dbscan_cluster_interior_labels: [Nx3] array containing 
%                                              information on the compact 
%                                              domains at the nucleus 
%                                              interior. The coordinates 
%                                              of the locs in the first 
%                                              two columns and the index 
%                                              of the compact domain they 
%                                              belong to in the third 
%                                              column
%       - model_lad_segment_locs: Localizations found to belong to clusters
%                                 defining modeled LADs
%       - model_lad_segment_normals: Vectors indicating the orientation of
%                                    the model LAD orthogonal to the
%                                    tangent
%       - model_lad_segment_thickness: Distance from nucleus boundary of
%                                      the model LAD segments
%   periphery_thresh: Cutoff for what percentage of the polygon to consider
%                     as the periphery. Ex. 0.15 results in a periphery
%                     region with an area equivalent to 15% of the original
%                     polygon shape
% -------------------------------------------------------------------------
% Modified from code (v4.0, 04/05/2024) provided by the Shenoy Lab
% Code written by:
%   Aayush Kant         Shenoy lab, University of Pennsylvania (USA)
%          Current: Kant lab, Indian Institute of Technology Delhi (India)
% Code edited by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   A. Kant, Z. Guo, V. Vinayak, M. V. Neguembor, W. S. Li, V. Agrawal, E. 
%   Pujadas, L. Almassalha, V. Backman, M. Lakadamyali, M. P. Cosma, V. B. 
%   Shenoy, Active transcription and epigenetic reactions synergistically 
%   regulate meso-scale genomic organization. Nature Communications 15, 
%   4338 (2024).
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------

%% Plot DBSCAN and model LAD data
function plot_OSNAP_DBSCAN_LAD(filepath,options)
arguments
    filepath string
    options.data struct = [];
    options.periphery_thresh double = 0.15;
end
    if ~isempty(options.data)
        data = options.data;
    else
        vars_to_load = {...
            'boundary',...
            'components',...
            'polygon',...
            'locs_dbscan_cluster_periphery_labels',...
            'locs_dbscan_cluster_interior_labels',...
            'model_lad_segment_locs',...
            'model_lad_segment_normals',...
            'model_lad_segment_thickness'};
        data = load_variables_OSNAP(filepath,vars_to_load);
    end
    warning('off','MATLAB:polyshape:repairedBySimplify')
    [dir_name,name,~] = fileparts(filepath);
    map_dir = replace(dir_name,"_nucleus_data","_nucleus_analysis_images");
    if ~exist(map_dir, 'dir')
        mkdir(map_dir)
    end
    save_path = fullfile(map_dir,name+"_nucleus_analysis.png");
    % Setup scale bar 
    scale_bar_size_x = 2; % 2 um scale bar
    scale_bar_size_y = 0.5;
    nm_to_um_x = scale_bar_size_x*1000;
    nm_to_um_y = scale_bar_size_y*1000;
    % Plot
    f = figure('visible','off','Position',[200 200 900 750]);
    tiledlayout(1,1,"tilespacing","none");
    nexttile();
    polygon_outer = translate(data.polygon,-(max(data.polygon.Vertices) - range(data.polygon.Vertices)/2));
    p1 = plot(polygon_outer,'FaceColor',[0.6484    0.8047    0.9062],'FaceAlpha',1); hold on
    polygon_inner = scale(polygon_outer,(1-options.periphery_thresh));
    p2 = plot(polygon_inner,'FaceColor',[0.7188    0.6367    0.6250],'FaceAlpha',1); hold on
    % plot(data.boundary(:,1),data.boundary(:,2),'r-','LineWidth',2); hold on
    locs_dbscan_cluster_periphery_labels = normalize_points(data.locs_dbscan_cluster_periphery_labels(:,1:2),data.components,data.boundary);
    locs_dbscan_cluster_interior_labels = normalize_points(data.locs_dbscan_cluster_interior_labels(:,1:2),data.components,data.boundary);
    model_lad_segment_locs = normalize_points(data.model_lad_segment_locs(:,[1 2]),data.components,data.boundary);
    model_lad_segment_normals = data.model_lad_segment_normals*data.components;
    s1 = scatter(locs_dbscan_cluster_periphery_labels(:,1),locs_dbscan_cluster_periphery_labels(:,2),0.5,'.g'); hold on
    s2 = scatter(locs_dbscan_cluster_interior_labels(:,1),locs_dbscan_cluster_interior_labels(:,2),0.5,'.k'); hold on
    for i = 1:length(model_lad_segment_locs(:,1))-1
        point1 = model_lad_segment_locs(i,:);
        point2 = model_lad_segment_locs(i+1,:);
        point4 = point1 + data.model_lad_segment_thickness(i)*model_lad_segment_normals(i,:);
        point3 = point2 + data.model_lad_segment_thickness(i)*model_lad_segment_normals(i+1,:);
        vertices = unique([[point1(1);point2(1);point3(1);point4(1);point1(1)],...
            [point1(2);point2(2);point3(2);point4(2);point1(2)]],'rows');
        pgon = polyshape(vertices);
        if i==1
            p3 = plot(pgon,'FaceColor','k');hold on
        else
            plot(pgon,'FaceColor','k');hold on
        end
    end
    axis equal
    x_bounds=xlim();
    y_bounds=ylim();
    rectangle('Position',[x_bounds(1)+nm_to_um_y y_bounds(1)+nm_to_um_y nm_to_um_x nm_to_um_y],...
        'facecolor','w','edgecolor','none')
    ax = gca();
    set(ax,'Color','k','XTick',[],'YTick',[])
    lgd=legend([p1,p2,p3,s1,s2],...
        "Nucleus Periphery","Nucleus Interior", "Modeled Lamina-associated domains",...
        "Chromatin packing domains","Lamina-associated packing domains",...
        "Location","north","NumColumns",2,"Color","w","FontSize",12);
    lgd.Layout.Tile = 'north';
    drawnow()
    exportgraphics(f,save_path,"Resolution",300,"BackgroundColor","w");
    close(f);
    warning('on','MATLAB:polyshape:repairedBySimplify')
end

%% Normalize localizations
function points_norm = normalize_points(points,components,boundary)
    rotated_boundary = boundary*components;
    points_norm = points*components - min(boundary*components) - range(rotated_boundary)/2;
end