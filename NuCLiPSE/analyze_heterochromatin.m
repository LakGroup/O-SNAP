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
%If density threshold has been previously calculated, in the next line it can be used directly to save time%
%density_threshold = 0.014431999613957;
% OWflag = false;          % if true overwrites the density threshold for each nuclei <NOT RECOMMENDED>

function data = analyze_heterochromatin(filepath, density_threshold, eps, min_num, ellipse_inc, periphery_thresh)

    data = load_variables(filepath, {'x','y'});
    locs = [data.x data.y];

    if ~has_variables(filepath,{'boundary'})
        % find nuclear boundary
        bd_old = locs(boundary(locs, 0.5),:);
        bd = smoothdata(bd_old(1:length(bd_old)-1,:),'gaussian',10); clear bd_old
        bd(end+1,:) = bd(1,:);
        data.boundary = bd;
        clearvars bd_old
    else
        tmp = load(filepath,'boundary');
        bd = tmp.boundary;
        clearvars tmp
    end

    % number of localizations
    if ~has_variables(filepath,{'locs_number'})
        data.locs_number = length(locs(:,1));
    end

    % approximate size of nucleus
    if ~has_variables(filepath,{'nucleus_radius'})
        nucleus_radius = sqrt(polyarea(bd(:,1),bd(:,2))/pi);
        data.nucleus_radius = nucleus_radius;
    end

    if ~has_variables(filepath,{'locs_density'})
        if ~exist('nucleus_radius','var') && has_variables(filepath,'nucleus_radius')
            load(filepath,'nucleus_radius');
        end
        data.locs_density = length(locs(:,1))/(pi*nucleus_radius^2);
    end

    % dbscan
    if ~has_variables(filepath,{'hetero_flt_with_label'})
        if has_variables(filepath,{'voronoi_areas_all'})
            load(filepath,'voronoi_areas_all');
        end
        hetero = locs((1./voronoi_areas_all)>=density_threshold,:);
        labels = dbscan(hetero,eps,min_num);
        hetero_flt = removerows(hetero,'ind',find(labels == -1)); % filter background noise
        labels_flt = removerows(labels,'ind',find(labels == -1)); 
        hetero_flt_with_label = [hetero_flt,labels_flt];
        data.hetero_flt_with_label = hetero_flt_with_label;
    else
        load(filepath,'hetero_flt_with_label');
        hetero_flt = hetero_flt_with_label(:,1:2);
        labels_flt = hetero_flt_with_label(:,3);
    end
    clearvars density hetero labels hetero_flt_with_label
    
    % define LADs and non-LADs domain
    if ~has_variables(filepath,{'lads','non_lads'})
        labels_flt_shuffle = shuffleLabel(labels_flt); % shuffle label
        num_groups = length(unique(labels_flt_shuffle));
        lads_is = [];
        lads_index = [];
        non_lads_is = [];
        non_lads_index = [];
        cnt = 1;
        nl_cnt = 1;
        for i=1:num_groups
            grp = hetero_flt(labels_flt == i,:);
            % sampling from grp
            mid_point = [mean(grp(:,1)), mean(grp(:,2))];
            rid = randi([1 length(grp(:,1))],1,ceil(length(grp)/15));
            rnd_point = grp(rid,:);
            sam_points = [mid_point; rnd_point];
            dist_mat = pdist2(sam_points, bd);
            min_dist = min(min(dist_mat));
            if min_dist <= 0.05*data.nucleus_radius && length(grp(:,1))>=35
                % min_dist <= 0.05*approximate_radius && length(grp(:,1))>=35
                lads_is = [lads_is; find(labels_flt == i)]; % old index
                lads_index = [lads_index; cnt*ones(length(grp(:,1)),1)]; % new index
                cnt = cnt + 1;
            elseif length(grp(:,1))>=35
                temp_bd = grp(boundary(grp,0.3),:);
                % temp_area = polyarea(temp_bd(:,1),temp_bd(:,2));
                % temp_r = sqrt(temp_area/pi);
                non_lads_is = [non_lads_is; find(labels_flt == i)];
                non_lads_index = [non_lads_index; nl_cnt*ones(length(grp(:,1)),1)];
                nl_cnt = nl_cnt + 1;
                clear temp_r temp_area temp_bd
            end
        end
        non_lads = hetero_flt(non_lads_is,:);
        lads = hetero_flt(lads_is,:);
        data.lads = [lads,lads_index];
        data.non_lads = [non_lads,non_lads_index];
        data.lads2total = length(lads)/(length(lads)+length(non_lads));
        clearvars hetero_flt_with_label hetero_flt labels_flt labels_flt_shuffle
    else
        load(filepath,'lads');
        load(filepath,'non_lads');
        lads_index = lads(:,3);
        lads = lads(:,1:2);
        non_lads_index = non_lads(:,3);
        non_lads = non_lads(:,1:2);
    end

    if ~has_variables(filepath,{'lads_n_locs','lads_center','lads_density'})
        if has_variables(filepath,{'lads_area'})
            load(filepath,'lads_area');
            num_of_lads = length(unique(lads_index));
            data.lads_n_locs = arrayfun(@(x) length(lads(lads_index==x,:)), 1:num_of_lads,'uni',1)';
            data.lads_density = data.lads_n_locs./lads_area;
        else
            if ~exist('lads','var') && has_variables(filepath,'lads')
                load(filepath,'lads');
                lads_index = lads(:,3);
                lads = lads(:,1:2);
            end
            % analysis of inner ladschromatin domain size
            num_of_lads = length(unique(lads_index));
            lads_n_locs = zeros(num_of_lads,1);
            lads_center = zeros(num_of_lads,2);
            lads_area = zeros(num_of_lads,1);
            lads_density = zeros(num_of_lads,1);
            for i=1:num_of_lads
                grp = lads(lads_index==i,:);
                lads_n_locs(i) = length(grp);
                lads_center(i,:) = [mean(grp(:,1)),mean(grp(:,2))];
                [k,lads_area_i]=boundary(grp,0.3);
                lads_area(i) = lads_area_i;
                lads_density(i) = lads_n_locs(i)./lads_area_i;
                data.domain_boundary{i,1}=grp(k,:);
                clear grp 
            end
            data.lads_n_locs = lads_n_locs;
            data.lads_center = lads_center;
            data.lads_area = lads_area;
            data.lads_density = lads_density;
        end
    end
    
    if ~has_variables(filepath,{'hetero_n_locs','hetero_center','hetero_density'})
        if has_variables(filepath,{'hetero_radius'})
            load(filepath,'hetero_radius');
            num_of_hetero = length(unique(non_lads_index));
            data.hetero_n_locs = arrayfun(@(x) length(non_lads(non_lads_index==x,:)), 1:num_of_hetero,'uni',1)';
            data.hetero_density = data.hetero_n_locs./(pi*hetero_radius.^2);
        else
            if ~exist('non_lads','var') && has_variables(filepath,'non_lads')
                load(filepath,'non_lads');
                non_lads_index = non_lads(:,3);
                non_lads = non_lads(:,1:2);
            end
            % analysis of inner Heterochromatin domain size
            num_of_hetero = length(unique(non_lads_index));
            hetero_n_locs = zeros(num_of_hetero,1);
            hetero_center = zeros(num_of_hetero,2);
            hetero_radius = zeros(num_of_hetero,1);
            hetero_density = zeros(num_of_hetero,1);
            for i=1:num_of_hetero
                grp = non_lads(non_lads_index==i,:);
                hetero_n_locs(i) = length(grp);
                hetero_center(i,:) = [mean(grp(:,1)),mean(grp(:,2))];
                [k,hetero_area]=boundary(grp,0.3);
                hetero_radius(i) = sqrt(hetero_area/pi);
                hetero_density(i) = hetero_n_locs(i)./hetero_area;
                data.domain_boundary{i,1}=grp(k,:);
                clear grp 
            end
            data.hetero_n_locs = hetero_n_locs;
            data.hetero_center = hetero_center;
            data.hetero_radius = hetero_radius;
            data.hetero_density = hetero_density;
        end
    else
        load(filepath,'hetero_center');
    end

    % find spacing
    if ~has_variables(filepath,{'spacing'})
        N_neighbour = 5;
        dist_mat = mink(pdist2(hetero_center,hetero_center),N_neighbour+1,2);
        spacing = sum(dist_mat,2)/N_neighbour;
        data.spacing = spacing;
        clearvars dist_mat N_neighbour
    end
    
    %% compute lads thickness
    if ~has_variables(filepath,{'seg_normals','lad_thickness','seg_locs','lads_seg_thickness'})
        % create accumulate center location
        bd_acc = zeros(length(bd(:,1)),1);
        for i=2:length(bd(:,1))
            bd_acc(i) = bd_acc(i-1) + norm(bd(i,:)-bd(i-1,:));
        end
        % create unit vector
        tangents=diff(bd);
        u_tangents=tangents./sqrt(sum(tangents.^2,2));
        normals = zeros(length(bd(:,1))-1,2);
        for i=1:length(bd(:,1))-1
            normals(i,:)=[0 -1;1 0]*u_tangents(i,:)';
        end
        normals(end+1,:) = normals(1,:);
        % define LADs segments
        seg_num = 50;
        seg_len = floor(length(bd(:,1))/seg_num);
        % calculate lads thickness within each segment
        bd_idx = 1;
        indent_len = 2000;
        seg_len_store = [];
        seg_area = [];
        lads_seg_thickness = [];
        seg_locs = [];
        seg_normals = [];
        if mod(length(bd(:,1)),seg_num)==0
            num_iters = seg_num;
        else
            num_iters = seg_num+1;
        end
        for seg_idx = 1:num_iters
            point1 = bd(bd_idx,:);
            point4 = point1 + indent_len*normals(bd_idx,:);
            seg_normals = [seg_normals;normals(bd_idx,:)];
            if seg_idx<num_iters
                bd_idx = bd_idx +seg_len;
            else
                bd_idx = bd_idx + length(bd(:,1))-seg_num*seg_len-1;
            end
            point2 = bd(bd_idx,:);
            point3 = point2 + indent_len*normals(bd_idx,:);
            in = inpolygon(lads(:,1),lads(:,2),...
                [point1(1);point2(1);point3(1);point4(1);point1(1)],...
                [point1(2);point2(2);point3(2);point4(2);point1(2)]);
            grp = lads(in,:);
            bd_grp = grp(boundary(grp,0.5),:);
            area = polyarea(bd_grp(:,1),bd_grp(:,2));
            seg_area = [seg_area;area];
            seg_len_store = [seg_len_store;norm(point1-point2)];
            lads_seg_thickness = [lads_seg_thickness;area/norm(point1-point2)];
            seg_locs = [seg_locs;point1];
            
        end
        seg_locs(end+1,:) = seg_locs(1,:);
        seg_normals(end+1,:) = seg_normals(1,:);
        lads_seg_thickness(lads_seg_thickness <= 0) = 1;
        data.seg_normals = seg_normals;
        data.seg_locs = seg_locs;
        data.lads_seg_len = seg_len_store;
        data.lads_seg_thickness = lads_seg_thickness;
        data.lad_thickness = sum(seg_area)/sum(seg_len_store);
    else
        load(filepath, 'seg_normals','seg_locs','lads_seg_thickness');
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
        load(filepath,'x_length','y_length','locs_norm','components');
    end

    %% Generate polygon representation of data
    if ~has_variables(filepath,{'polygon'})
        polygon = polyshape(rotated_boundary(1:end-1,1),rotated_boundary(1:end-1,2),'Simplify',false); % Create a polygon of the boundary points (which will be completely filled).
        [Cx, Cy] = centroid(polygon);
        polygon = translate(polygon, -[Cx Cy]);
        data.polygon = polygon;
    else
        load(filepath,'polygon');
    end

    %% Calculate the boundary descriptors.
    if ~has_variables(filepath,{'elastic_energy','bending_energy','border_curvature'})
        [data.elastic_energy, data.bending_energy] = Energy(rotated_boundary); % Calculate the elastic energy and bending energy for each border.
        data.border_curvature = Curvature(rotated_boundary); % Only calculate the curvature for the outer boundary of the point cloud. 
    end

    %% calculate radius density
    if ~has_variables(filepath,{'radial_density'})
        radial_density = calculate_radial_density(locs_norm, [x_length y_length], ellipse_inc);
        data.radial_density = radial_density;
    else
        load(filepath,'radial_density')
        if size(radial_density,1) < floor(1/ellipse_inc)
            radial_density = calculate_radial_density(locs_norm, [x_length y_length], ellipse_inc);
            data.radial_density = radial_density;
        end
    end

    %% repeat radial_density for radial_hetero_density
    if ~has_variables(filepath,{'radial_hetero_density'})
        radial_hetero_density = calculate_radial_density(normalize_points(hetero_center,components,bd,rotated_boundary),...
            [x_length y_length], ellipse_inc);
        data.radial_hetero_density = radial_hetero_density;
    else
        load(filepath,'radial_hetero_density')
        if size(radial_hetero_density,1) ~= floor(1/ellipse_inc)
            radial_hetero_density = calculate_radial_density(normalize_points(hetero_center,components,bd,rotated_boundary),... 
            [x_length y_length], ellipse_inc);
            data.radial_hetero_density = radial_hetero_density;
        end
    end

    %% calculate periphery/interior density
    if ~has_variables(filepath,{'interior_density','periphery_density'})
        [interior_density, periphery_density]  = calculate_periphery_density(locs_norm,polygon,...
            periphery_thresh);
        data.interior_density = interior_density;
        data.periphery_density = periphery_density;
    end

    %% calculate periphery/interior density of Heterochromatin clusters
    if ~has_variables(filepath,{'interior_hetero_density','periphery_hetero_density'})
        [interior_hetero_density, periphery_hetero_density]  = calculate_periphery_density(normalize_points(hetero_center,components,bd,rotated_boundary),...
            polygon, periphery_thresh);
        data.interior_hetero_density = interior_hetero_density;
        data.periphery_hetero_density = periphery_hetero_density;
    end

    [dir_name,name,~] = fileparts(filepath);
    plot_file = fullfile(replace(dir_name,"voronoi_data","nuclei_images"),name+"_nucleus_analysis.png");
    % if ~exist(plot_file,'file')
        f = figure('visible','off');
        polygon_outer = translate(polygon,-(max(polygon.Vertices) - range(polygon.Vertices)/2));
        plot(polygon_outer); hold on
        polygon_inner = scale(polygon_outer,(1-periphery_thresh));
        plot(polygon_inner); hold on
        % plot(bd(:,1),bd(:,2),'r-','LineWidth',2); hold on
        lads = normalize_points(lads,components,bd,rotated_boundary);
        non_lads = normalize_points(non_lads,components,bd,rotated_boundary);
        seg_locs = normalize_points(seg_locs,components,bd,rotated_boundary);
        seg_normals = seg_normals*components;
        scatter(lads(:,1),lads(:,2),0.5,'g'); hold on
        scatter(non_lads(:,1),non_lads(:,2),0.5,'k'); hold on
        legend('off')
        for i = 1:length(seg_locs(:,1))-1
            point1 = seg_locs(i,:);
            point2 = seg_locs(i+1,:);
            point4 = point1 + lads_seg_thickness(i)*seg_normals(i,:);
            point3 = point2 + lads_seg_thickness(i)*seg_normals(i+1,:);
            pgon = polyshape([point1(1);point2(1);point3(1);point4(1);point1(1)],...
                [point1(2);point2(2);point3(2);point4(2);point1(2)]);
            plot(pgon,'FaceColor','k');hold on
            axis equal
        end
        drawnow()
        saveas(f,plot_file);
        close(f);
    % end
    clearvars -except filepath data
    save_voronoi_data(filepath, data);
end

function points_norm = normalize_points(points,components,boundary,rotated_boundary)
    points_norm = points*components - min(boundary*components) - range(rotated_boundary)/2;
end
