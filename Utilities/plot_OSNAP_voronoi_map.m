% -------------------------------------------------------------------------
% plot_OSNAP_voronoi_map.m
% -------------------------------------------------------------------------
% Creates the Voronoi density maps following Voronoi segmentation. They are
% automatically saved to a directory "OSNAP_nucleus_images". Data is passed
% from OSNAP files through the data struct, which is assumed to be well
% formatted to include the necessary fields related to the nucleus:
%   - vertices:  the coordinates of each unique vertex 
%   - faces: how to connect vertices
%   - voronoi_areas_all: area associated the Voronoi cell of each loc
% Two versions of the Voronoi density maps are generated:
%   - Voronoi area (no suffix - ".png"): Color-coding based upon the raw
%                                        values of the areas for the
%                                        Voronoi cells in nm^2
%   - Reduced Voronoi density ("_reduced.png"): Color-coding based upon the 
%                                        values of the areas for the
%                                        Voronoi cells in nm^2 divided by
%                                        the average Voronoi area of the
%                                        entire nucleus. Take the inverse
%                                        to get reduced density.
% -------------------------------------------------------------------------
% Input:
%   filepath: Table with details on the samples of interest
%                           to generate feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
%   data: A struct array with information on the Voronoi data for a sample
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
% Example on how to use it:
%   plot_OSNAP_voronoi_map("D:\Analysis\ExperimentA\nucleus_X.mat",...
%                          data_nucleus_X,...
%                          0,3,...
%                          1,3)
% Input:
%   OSNAP_sample_file_list: Table with details on the samples of interest
%                           to generate feature information from, including 
%                           the file path, replicate, phenotype, and sample
%                           identifier
% Options:
%   n_processes: Number of cores to run parallel processes on
%   min_log_reduced_vor_density: Lower bound for filtering and visualizing reduced 
%                        Voronoi density values
%   max_log_reduced_vor_density: Upper bound for filtering and visualizing reduced
%                        Voronoi density values
%   min_log_vor_area: Lower bound for filtering and visualizing Voronoi
%                     area values
%   max_log_vor_area: Upper bound for filtering and visualizing Voronoi
%                     area values
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

% Plot the voronoi maps
function plot_OSNAP_voronoi_map(filepath,data,options)
arguments
    filepath string
    data struct
    options.min_log_reduced_vor_density double = 0;
    options.max_log_reduced_vor_density double = 3;
    options.min_log_vor_area double = 1;
    options.max_log_vor_area double = 3;
    options.CT;
end
    [dir_name,name,~] = fileparts(filepath);
    voronoi_density_map_dir = replace(dir_name,"_nucleus_data","_voronoi_density_maps");
    if ~exist(voronoi_density_map_dir, 'dir')
        mkdir(voronoi_density_map_dir)
    end
    reduced_voronoi_density_map_dir = replace(dir_name,"_nucleus_data","_reduced_voronoi_density_maps");
    if ~exist(reduced_voronoi_density_map_dir, 'dir')
        mkdir(reduced_voronoi_density_map_dir)
    end
    % Setup limits
    [x_min,x_max] = bounds(data.x);
    [y_min,y_max] = bounds(data.y);
    % Setup scale bar 
    scale_bar_size_x = 2;
    scale_bar_size_y = 0.5;
    nm_to_um_x = scale_bar_size_x*1000;
    nm_to_um_y = scale_bar_size_y*1000;
    % % Color maps
    CT = interp1([0 2 4 6 8 9]/9,[0 0 0; 0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.3 0 0],linspace(0,1,256),'linear','extrap'); % jet - not recommended
    % % A
    % CT = interp1(linspace(0,1,7),[0.0 0.0 0.0; 0.1216 0.2275 0.2902; 0.3098 0.4118 0.4902; 0.502 0.5686 0.6275; 0.6784 0.7098 0.7569; 0.8431 0.8471 0.8784; 1.0 0.9765 0.9922],linspace(0,1,256),'linear','extrap');
    % % B
    % CT = interp1(linspace(0,1,7),[0.0 0.0 0.0; 0.1412 0.1804 0.5059; 0.1804 0.4157 0.6824; 0.2275 0.6196 0.749; 0.4627 0.7765 0.7412; 0.749 0.898 0.7373; 1.0 1.0 0.8],linspace(0,1,256),'linear','extrap');
    % % C
    % CT = interp1(linspace(0,1,6),[0.0 0.0 0.0; 0.4784 0.0039 0.4667; 0.7725 0.1059 0.5412; 0.9686 0.4078 0.6314; 0.9843 0.7059 0.7255; 0.9961 0.9216 0.8863],linspace(0,1,256),'linear','extrap');
    % % D
    % CT = interp1(linspace(0,1,7),[0.0 0.0 0.0; 0.0627 0.251 0.1412; 0.0667 0.4706 0.2471; 0.251 0.651 0.3451; 0.5255 0.7843 0.4745; 0.7647 0.902 0.6118; 1.0 1.0 0.8],linspace(0,1,256),'linear','extrap');

    %% Plot Voronoi area
    plot_faces(data.vertices,data.faces,data.voronoi_areas_all,flipud(CT),...
        x_min,x_max,y_min,y_max,options.min_log_vor_area,options.max_log_vor_area,...
        nm_to_um_x,nm_to_um_y,...
        'Log Voronoi Area Map','$\textrm{Log}_{10}$(Voronoi polygon area) [$\textrm{nm}^{2}$]',...
        fullfile(voronoi_density_map_dir,name+".png"));
    %% Plot reduced voronoi density
    plot_faces(data.vertices,data.faces,1./(data.voronoi_areas_all/mean(data.voronoi_areas_all,'omitnan')),CT,...
        x_min,x_max,y_min,y_max,options.min_log_reduced_vor_density,options.max_log_reduced_vor_density,...
        nm_to_um_x,nm_to_um_y,...
        'Reduced Log Voronoi Density Map','Reduced $\textrm{Log}_{10}$(Voronoi polygon density) [$\textrm{nm}^{-2}$]',...
        fullfile(reduced_voronoi_density_map_dir,name+"_reduced.png"));
end

% helper function to plot
function plot_faces(vertices,faces,values,CT,...
    x_min,x_max,y_min,y_max,c_min,c_max,...
    nm_to_um_x,nm_to_um_y,...
    f_label,c_label,...
    savepath)
    c_values = log10(values);
    idx = values >= c_min;
    c_values = c_values(idx);
    c_values = min(c_max,c_values);
    c_ticks = c_min:0.5:c_max;
    f = figure('visible','off','name',f_label,'NumberTitle','off','color','k','units','normalized',...
        'position',[0.25 0.15 0.6 0.7]);
    ax = axes(f);
    set(ax,'colormap',CT,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
    set(f, 'InvertHardCopy', 'off'); 
    patch('Vertices',vertices,'Faces',faces(idx,:),'FaceVertexCData',c_values,'FaceColor','flat','edgecolor','none')    
    axis equal
    xlim(ax,[x_min x_max])
    ylim(ax,[y_min y_max])    
    h = colorbar('Position',[0.9 0.1 0.02 0.8],...
        'TickLabelInterpreter','latex',...
        'FontSize',14,...
        'Color','w',...
        'Ticks',c_ticks,...
        'TickLabels',arrayfun(@(elem) sprintf('%.2f',elem), c_ticks,'uni',0));    
    pbaspect([1 1 1])
    ylabel(h,c_label,'FontSize',14,'interpreter','latex');
    caxis([c_min c_max]);
    rectangle('Position',[x_min y_min nm_to_um_x nm_to_um_y],'facecolor','w')
    exportgraphics(ax,savepath,"Resolution",300,"BackgroundColor","k");
    close(f)
end