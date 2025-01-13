% This code works as a standard code for the params extraction of 
% STORM images (@Shenoy_lab) Version 4.0 (last update April 5, 2024)

% Notes: This file is used to extract key params of nucleus as follows:
% 1. Inner Heterochromatin Size 2. LADs thickness 3. Nuclear a
% All variables are saved into a cell 'data', which is written to the
% system

% CHANGELOG (for v4)
% Saving the boundary of domains into the data class
% Option to overwrite the global threshold for individual nuclei
% Made to run in parallel from loading voronoi analysis data

% density storage
% TSA: 0.017972
% GSK dss: 1 threshold: 0.030912 --> 0.4
% Stiffness dss: 1 threshold: 0.021197
% teno stiffness
% If density threshold has been previously calculated, in the next line it can be used directly to save time%
% density_threshold = 0.014431999613957;
% OWflag = false;          % if true overwrites the density threshold for each nuclei <NOT RECOMMENDED>

function data = generate_morphometrics_dbscan_radial_data(file_path, density_threshold, options)
arguments
    file_path string
    density_threshold double
    options.eps double = 20;
    options.min_num double = 3;
    options.ellipse_inc double = 0.1;
    options.periphery_thresh double = 0.15;
    options.plot logical = true;
end

    data = load_variables(file_path, {'x','y'});
    locs = [data.x data.y];

    %% find nuclear boundary
    if ~has_variables(file_path,{'boundary'})
        bd_old = locs(boundary(locs, 0.5),:);
        bd = smoothdata(bd_old(1:length(bd_old)-1,:),'gaussian',10);
        clearvars bd_old
        bd(end+1,:) = bd(1,:);
        data.boundary = bd;
    else
        tmp = load(file_path,'boundary');
        bd = tmp.boundary;
        clearvars tmp
    end

    %% number of localizations
    if ~has_variables(file_path,{'locs_number'})
        data.locs_number = length(locs(:,1));
    end

    %% approximate size of nucleus
    if ~has_variables(file_path,{'nucleus_radius'})
        nucleus_radius = sqrt(polyarea(bd(:,1),bd(:,2))/pi);
        data.nucleus_radius = nucleus_radius;
    end

    %% localization density with respect to nucleus area
    if ~has_variables(file_path,{'locs_density'})
        if ~exist('nucleus_radius','var') && has_variables(file_path,'nucleus_radius')
            load(file_path,'nucleus_radius');
        end
        data.locs_density = length(locs(:,1))/(pi*nucleus_radius^2);
    end

    %% dbscan
    if ~has_variables(file_path,{'locs_dbscan_cluster_labels'})
        if has_variables(file_path,{'voronoi_areas_all'})
            load(file_path,'voronoi_areas_all');
            locs_to_cluster = locs((1./voronoi_areas_all)>=density_threshold,:);
            labels = dbscan(locs_to_cluster,options.eps,options.min_num);
            hetero_flt = removerows(locs_to_cluster,'ind',find(labels == -1)); % filter background noise
            labels_flt = removerows(labels,'ind',find(labels == -1)); 
            locs_dbscan_cluster_labels = [hetero_flt,labels_flt];
            data.locs_dbscan_cluster_labels = locs_dbscan_cluster_labels;
        end
    else
        load(file_path,'locs_dbscan_cluster_labels');
        hetero_flt = locs_dbscan_cluster_labels(:,1:2);
        labels_flt = locs_dbscan_cluster_labels(:,3);
    end
    clearvars density hetero labels locs_dbscan_cluster_labels
    
    %% define dbscan clusters at the periphery (LADs) and interior (non-LADs)
    if ~has_variables(file_path,{...
            'locs_dbscan_cluster_periphery_labels',...
            'locs_dbscan_cluster_interior_labels'})
        labels_flt_shuffle = shuffle_label(labels_flt); % shuffle label
        locs_dbscan_cluster_periphery_is = [];
        locs_dbscan_cluster_periphery_index = [];
        locs_dbscan_cluster_interior_labels_is = [];
        locs_dbscan_cluster_interior_labels_index = [];
        cnt = 1;
        nl_cnt = 1;
        for i=1:max(labels_flt_shuffle)
            grp = hetero_flt(labels_flt == i,:);
            % sampling from grp
            mid_point = [mean(grp(:,1)), mean(grp(:,2))];
            rid = randi([1 length(grp(:,1))],1,ceil(length(grp)/15));
            rnd_point = grp(rid,:);
            sam_points = [mid_point; rnd_point];
            dist_mat = pdist2(sam_points, bd);
            min_dist = min(min(dist_mat));
            if min_dist <= 0.05*data.nucleus_radius && length(grp(:,1))>=35
                locs_dbscan_cluster_periphery_is = [locs_dbscan_cluster_periphery_is; find(labels_flt == i)]; % old index
                locs_dbscan_cluster_periphery_index = [locs_dbscan_cluster_periphery_index; cnt*ones(length(grp(:,1)),1)]; % new index
                cnt = cnt + 1;
            elseif length(grp(:,1))>=35
                temp_bd = grp(boundary(grp,0.3),:);
                locs_dbscan_cluster_interior_labels_is = [locs_dbscan_cluster_interior_labels_is; find(labels_flt == i)];
                locs_dbscan_cluster_interior_labels_index = [locs_dbscan_cluster_interior_labels_index; nl_cnt*ones(length(grp(:,1)),1)];
                nl_cnt = nl_cnt + 1;
                clear temp_bd
            end
        end
        locs_dbscan_cluster_interior_labels = hetero_flt(locs_dbscan_cluster_interior_labels_is,:);
        locs_dbscan_cluster_periphery_labels = hetero_flt(locs_dbscan_cluster_periphery_is,:);
        data.locs_dbscan_cluster_interior_labels = [locs_dbscan_cluster_interior_labels,locs_dbscan_cluster_interior_labels_index];
        data.locs_dbscan_cluster_periphery_labels = [locs_dbscan_cluster_periphery_labels,locs_dbscan_cluster_periphery_index];
        clearvars locs_dbscan_cluster_labels hetero_flt labels_flt labels_flt_shuffle
    else
        load(file_path,'locs_dbscan_cluster_periphery_labels');
        load(file_path,'locs_dbscan_cluster_interior_labels');
        locs_dbscan_cluster_periphery_index = locs_dbscan_cluster_periphery_labels(:,3);
        locs_dbscan_cluster_periphery_labels = locs_dbscan_cluster_periphery_labels(:,1:2);
        locs_dbscan_cluster_interior_labels_index = locs_dbscan_cluster_interior_labels(:,3);
        locs_dbscan_cluster_interior_labels = locs_dbscan_cluster_interior_labels(:,1:2);
    end

    %% analysis of lads chromatin domain size
    if ~has_variables(file_path,{...
            'periphery_dbscan_cluster_n_locs',...
            'periphery_dbscan_cluster_center',...
            'periphery_dbscan_cluster_density'})
        if has_variables(file_path,{'periphery_dbscan_cluster_radius'})
            load(file_path,'periphery_dbscan_cluster_radius');
            n_locs_dbscan_cluster_periphery = max(locs_dbscan_cluster_periphery_labels(:,3));
            data.periphery_dbscan_cluster_n_locs = zeros(n_locs_dbscan_cluster_periphery,1);
            for i=1:n_locs_dbscan_cluster_periphery
                data.periphery_dbscan_cluster_n_locs(i) = sum(locs_dbscan_cluster_periphery_index==i);
            end
            data.periphery_dbscan_cluster_density = data.periphery_dbscan_cluster_n_locs./(periphery_dbscan_cluster_radius.^2*pi);
        else
            if ~exist('locs_dbscan_cluster_periphery_labels','var') && has_variables(file_path,'locs_dbscan_cluster_periphery_labels')
                load(file_path,'locs_dbscan_cluster_periphery_labels');
                locs_dbscan_cluster_periphery_index = locs_dbscan_cluster_periphery_labels(:,3);
                locs_dbscan_cluster_periphery_labels = locs_dbscan_cluster_periphery_labels(:,1:2);
            end
            n_locs_dbscan_cluster_periphery = max(locs_dbscan_cluster_periphery_index);
            periphery_dbscan_cluster_n_locs = zeros(n_locs_dbscan_cluster_periphery,1);
            periphery_dbscan_cluster_center = zeros(n_locs_dbscan_cluster_periphery,2);
            periphery_dbscan_cluster_area = zeros(n_locs_dbscan_cluster_periphery,1);
            periphery_dbscan_cluster_density = zeros(n_locs_dbscan_cluster_periphery,1);
            for i=1:n_locs_dbscan_cluster_periphery
                grp = locs_dbscan_cluster_periphery_labels(locs_dbscan_cluster_periphery_index==i,:);
                periphery_dbscan_cluster_n_locs(i) = length(grp);
                periphery_dbscan_cluster_center(i,:) = [mean(grp(:,1)),mean(grp(:,2))];
                [k,periphery_dbscan_cluster_area_i] = boundary(grp,0.3);
                periphery_dbscan_cluster_area(i) = periphery_dbscan_cluster_area_i;
                periphery_dbscan_cluster_density(i) = periphery_dbscan_cluster_n_locs(i)./periphery_dbscan_cluster_area_i;
                data.domain_boundary{i,1}=grp(k,:);
                clear grp 
            end
            data.periphery_dbscan_cluster_n_locs = periphery_dbscan_cluster_n_locs;
            data.periphery_dbscan_cluster_center = periphery_dbscan_cluster_center;
            data.periphery_dbscan_cluster_radius = sqrt(periphery_dbscan_cluster_area/pi);
            data.periphery_dbscan_cluster_density = periphery_dbscan_cluster_density;
        end
    end
    
    %% analysis of inner heterochromatin domain size
    if ~has_variables(file_path,{...
            'interior_dbscan_cluster_n_locs',...
            'interior_dbscan_cluster_center',...
            'interior_dbscan_cluster_density'})
        if has_variables(file_path,{'interior_dbscan_cluster_radius'})
            load(file_path,'interior_dbscan_cluster_radius');
            n_hetero = max(locs_dbscan_cluster_interior_labels(:,3));
            data.interior_dbscan_cluster_n_locs = zeros(n_hetero,1);
            for i=1:n_hetero
                data.interior_dbscan_cluster_n_locs(i) = sum(locs_dbscan_cluster_interior_labels_index==i);
            end
            data.interior_dbscan_cluster_density = data.interior_dbscan_cluster_n_locs./(pi*interior_dbscan_cluster_radius.^2);
        else
            if ~exist('locs_dbscan_cluster_interior_labels','var') && has_variables(file_path,'locs_dbscan_cluster_interior_labels')
                load(file_path,'locs_dbscan_cluster_interior_labels');
                locs_dbscan_cluster_interior_labels_index = locs_dbscan_cluster_interior_labels(:,3);
                locs_dbscan_cluster_interior_labels = locs_dbscan_cluster_interior_labels(:,1:2);
            end
            n_hetero = max(locs_dbscan_cluster_interior_labels_index);
            interior_dbscan_cluster_n_locs = zeros(n_hetero,1);
            interior_dbscan_cluster_center = zeros(n_hetero,2);
            interior_dbscan_cluster_radius = zeros(n_hetero,1);
            interior_dbscan_cluster_density = zeros(n_hetero,1);
            for i=1:n_hetero
                grp = locs_dbscan_cluster_interior_labels(locs_dbscan_cluster_interior_labels_index==i,:);
                interior_dbscan_cluster_n_locs(i) = length(grp);
                interior_dbscan_cluster_center(i,:) = [mean(grp(:,1)),mean(grp(:,2))];
                [k,hetero_area]=boundary(grp,0.3);
                interior_dbscan_cluster_radius(i) = sqrt(hetero_area/pi);
                interior_dbscan_cluster_density(i) = interior_dbscan_cluster_n_locs(i)./hetero_area;
                data.domain_boundary{i,1}=grp(k,:);
                clear grp 
            end
            data.interior_dbscan_cluster_n_locs = interior_dbscan_cluster_n_locs;
            data.interior_dbscan_cluster_center = interior_dbscan_cluster_center;
            data.interior_dbscan_cluster_radius = interior_dbscan_cluster_radius;
            data.interior_dbscan_cluster_density = interior_dbscan_cluster_density;
        end
    else
        load(file_path,'interior_dbscan_cluster_center');
    end

    %% find  lad segment spacing
    if ~has_variables(file_path,{'model_lad_segment_spacing'})
        N_neighbour = 5;
        dist_mat = mink(pdist2(interior_dbscan_cluster_center,interior_dbscan_cluster_center),N_neighbour+1,2);
        model_lad_segment_spacing = sum(dist_mat,2)/N_neighbour;
        data.model_lad_segment_spacing = model_lad_segment_spacing;
        clearvars dist_mat N_neighbour
    end
    
    %% compute lads thicknesses
    if ~has_variables(file_path,{...
            'model_lad_thickness',...
            'model_lad_segment_normals',...
            'model_lad_segment_locs',...
            'model_lad_segment_length',...
            'model_lad_segment_thickness'})
        % create accumulate center location
        % bd_acc = cumsum([0;vecnorm(diff(bd),2,2)]);
        % create unit vector
        tangents = diff(bd);
        u_tangents=tangents./sqrt(sum(tangents.^2,2));
        normals = ([0 -1;1 0]*u_tangents(1:end,:)')';
        normals(end+1,:) = normals(1,:);
        % define LADs segments
        seg_num = 50;
        seg_len = floor(size(bd,1)/seg_num);
        % calculate lads thickness within each segment
        bd_idx = 1;
        indent_len = 2000;
        seg_len_store = [];
        seg_area = [];
        model_lad_segment_thickness = [];
        model_lad_segment_locs = [];
        model_lad_segment_normals = [];
        if mod(size(bd,1),seg_num)==0
            num_iters = seg_num;
        else
            num_iters = seg_num+1;
        end
        for seg_idx = 1:num_iters
            point1 = bd(bd_idx,:);
            point4 = point1 + indent_len*normals(bd_idx,:);
            model_lad_segment_normals = [model_lad_segment_normals;normals(bd_idx,:)];
            if seg_idx<num_iters
                bd_idx = bd_idx +seg_len;
            else
                bd_idx = bd_idx + size(bd,1)-seg_num*seg_len-1;
            end
            point2 = bd(bd_idx,:);
            point3 = point2 + indent_len*normals(bd_idx,:);
            in = inpolygon(locs_dbscan_cluster_periphery_labels(:,1),locs_dbscan_cluster_periphery_labels(:,2),...
                [point1(1);point2(1);point3(1);point4(1);point1(1)],...
                [point1(2);point2(2);point3(2);point4(2);point1(2)]);
            grp = locs_dbscan_cluster_periphery_labels(in,:);
            bd_grp = grp(boundary(grp,0.5),:);
            area = polyarea(bd_grp(:,1),bd_grp(:,2));
            seg_area = [seg_area;area];
            seg_len_store = [seg_len_store;norm(point1-point2)];
            model_lad_segment_thickness = [model_lad_segment_thickness;area/norm(point1-point2)];
            model_lad_segment_locs = [model_lad_segment_locs;point1];
            
        end
        model_lad_segment_locs(end+1,:) = model_lad_segment_locs(1,:);
        model_lad_segment_normals(end+1,:) = model_lad_segment_normals(1,:);
        model_lad_segment_thickness(model_lad_segment_thickness <= 0) = 1;
        data.model_lad_segment_normals = model_lad_segment_normals;
        data.model_lad_segment_locs = model_lad_segment_locs;
        data.model_lad_segment_length = seg_len_store;
        data.model_lad_segment_thickness = model_lad_segment_thickness;
        data.model_lad_thickness = sum(seg_area)/sum(seg_len_store);
    else
        load(file_path, 'model_lad_segment_normals','model_lad_segment_locs','model_lad_segment_thickness');
    end

    %% Perform a Principal Component Analysis to rotate the data (maximizing the length)
    if ~all(isfield(data,{'locs_norm','x_length','y_length','eigenvalues','components'}))
        [components,rotated_boundary,eigenvalues] = pca([bd(:,1) bd(:,2)]); % Uses the built-in Matlab pca routine.
        locs_norm = normalize_points(locs,components,bd,rotated_boundary);        
        x_length = max(locs_norm(:,1))-min(locs_norm(:,1)); % Determine the length of the cloud in x.
        y_length = max(locs_norm(:,2))-min(locs_norm(:,2)); % Determine the length of the cloud in y.
        data.locs_norm = locs_norm;
        data.eigenvalues = eigenvalues;
        data.x_length = x_length;
        data.y_length = y_length;
        data.components = components;
    else
        [~,rotated_boundary,~] = pca([bd(:,1) bd(:,2)]);
        load(file_path,'x_length','y_length','locs_norm','components');
    end

    %% Generate polygon representation of data
    if ~has_variables(file_path,{'polygon'})
        polygon = polyshape(rotated_boundary(1:end-1,1),rotated_boundary(1:end-1,2),'Simplify',false); % Create a polygon of the boundary points (which will be completely filled).
        [Cx, Cy] = centroid(polygon);
        polygon = translate(polygon, -[Cx Cy]);
        data.polygon = polygon;
    else
        load(file_path,'polygon');
    end

    %% calculate the boundary descriptors
    if ~has_variables(file_path,{'elastic_energy','bending_energy','border_curvature'})
        [data.elastic_energy, data.bending_energy] = Energy(rotated_boundary); % Calculate the elastic energy and bending energy for each border.
        data.border_curvature = Curvature(rotated_boundary); % Only calculate the curvature for the outer boundary of the point cloud. 
    end

    %% calculate radial localization density
    if ~has_variables(file_path,{'radial_loc_density'})
        radial_loc_density = calculate_radial_loc_density(locs_norm, [x_length y_length], options.ellipse_inc);
        data.radial_loc_density = radial_loc_density;
    else
        load(file_path,'radial_loc_density')
        if size(radial_loc_density,1) < floor(1/options.ellipse_inc)
            radial_loc_density = calculate_radial_loc_density(locs_norm, [x_length y_length], options.ellipse_inc);
            data.radial_loc_density = radial_loc_density;
        end
    end

    %% calculate radial dbscan cluster density
    if ~has_variables(file_path,{'radial_dbscan_cluster_density'})
        radial_dbscan_cluster_density = calculate_radial_loc_density(normalize_points(interior_dbscan_cluster_center,components,bd,rotated_boundary),...
            [x_length y_length], options.ellipse_inc);
        data.radial_dbscan_cluster_density = radial_dbscan_cluster_density;
    else
        load(file_path,'radial_dbscan_cluster_density')
        if size(radial_dbscan_cluster_density,1) ~= floor(1/options.ellipse_inc)
            radial_dbscan_cluster_density = calculate_radial_loc_density(normalize_points(interior_dbscan_cluster_center,components,bd,rotated_boundary),... 
            [x_length y_length], options.ellipse_inc);
            data.radial_dbscan_cluster_density = radial_dbscan_cluster_density;
        end
    end

    %% calculate periphery/interior density
    if ~has_variables(file_path,{'interior_loc_density','periphery_loc_density'})
        [interior_loc_density, periphery_loc_density]  = calculate_periphery_loc_density(locs_norm,polygon,...
            options.periphery_thresh);
        data.interior_loc_density = interior_loc_density;
        data.periphery_loc_density = periphery_loc_density;
    end

    %% calculate periphery/interior density of Heterochromatin clusters
    if ~has_variables(file_path,{'interior_cluster_density','periphery_cluster_density'})
        [interior_cluster_density, periphery_cluster_density]  = calculate_periphery_loc_density(normalize_points(interior_dbscan_cluster_center,components,bd,rotated_boundary),...
            polygon, options.periphery_thresh);
        data.interior_cluster_density = interior_cluster_density;
        data.periphery_cluster_density = periphery_cluster_density;
    end

    if options.plot
        [dir_name,name,~] = fileparts(file_path);
        map_dir = replace(dir_name,"SNAP_nucleus_data","SNAP_nucleus_images");
        if ~exist(map_dir, 'dir')
            mkdir(map_dir)
        end
        plot_file = fullfile(map_dir,name+"_nucleus_analysis.png");
        f = figure('visible','off');
        polygon_outer = translate(polygon,-(max(polygon.Vertices) - range(polygon.Vertices)/2));
        plot(polygon_outer); hold on
        polygon_inner = scale(polygon_outer,(1-options.periphery_thresh));
        plot(polygon_inner); hold on
        % plot(bd(:,1),bd(:,2),'r-','LineWidth',2); hold on
        locs_dbscan_cluster_periphery_labels = normalize_points(locs_dbscan_cluster_periphery_labels,components,bd,rotated_boundary);
        locs_dbscan_cluster_interior_labels = normalize_points(locs_dbscan_cluster_interior_labels,components,bd,rotated_boundary);
        model_lad_segment_locs = normalize_points(model_lad_segment_locs,components,bd,rotated_boundary);
        model_lad_segment_normals = model_lad_segment_normals*components;
        scatter(locs_dbscan_cluster_periphery_labels(:,1),locs_dbscan_cluster_periphery_labels(:,2),0.5,'g'); hold on
        scatter(locs_dbscan_cluster_interior_labels(:,1),locs_dbscan_cluster_interior_labels(:,2),0.5,'k'); hold on
        legend('off')
        for i = 1:length(model_lad_segment_locs(:,1))-1
            point1 = model_lad_segment_locs(i,:);
            point2 = model_lad_segment_locs(i+1,:);
            point4 = point1 + model_lad_segment_thickness(i)*model_lad_segment_normals(i,:);
            point3 = point2 + model_lad_segment_thickness(i)*model_lad_segment_normals(i+1,:);
            pgon = polyshape([point1(1);point2(1);point3(1);point4(1);point1(1)],...
                [point1(2);point2(2);point3(2);point4(2);point1(2)]);
            plot(pgon,'FaceColor','k');hold on
            axis equal
        end
        drawnow()
        saveas(f,plot_file);
        close(f);
    end
    clearvars -except file_path data
    save_SNAP_nucleus(file_path, data);
end

function points_norm = normalize_points(points,components,boundary,rotated_boundary)
    points_norm = points*components - min(boundary*components) - range(rotated_boundary)/2;
end

