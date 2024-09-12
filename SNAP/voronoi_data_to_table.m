function T = voronoi_data_to_table(data_info_table)
    %% parameters
    percentile_thresh = [40 70];
    %% define variables of interest
    vars_string = ["group","biological_replicate","name"];
    vars_morphological = ["area","major_axis","minor_axis","convex_area",...
        "gyration_radius","eigenvalues_ratio","eigenentropy",...
        "perimeter","convex_perimeter",...
        "aspect_ratio","form_ratio","rectangularity","circularity","convex_circularity",...
        "convexity","solidity","equivalent_diameter","fiber_length","fiber_width","curl",...
        "major_axis_norm_area","minor_axis_norm_area","bending_energy_norm_area"];
    vars_value = [...
        "nucleus_radius",...
        "lads2total",...
        "locs_density","locs_number",...
        "elastic_energy","bending_energy",...
        "periphery_density","interior_density",...
        "periphery_hetero_density","interior_hetero_density"];
    vars_chromatin = "voronoi_density_" + percentile_thresh;
    vars_distribution_names = [...
        "hetero_radius","hetero_density","hetero_n_locs",...
        "lads_area","lads_density","lads_n_locs","lads_seg_len","lads_seg_thickness",...
        "reduced_log_voronoi_density","spacing","border_curvature"];
    vars_distribution = reshape(vars_distribution_names + ["_median";"_std";"_skewness"],1,[]);
    vars_n_names = ["lads","non_lads","cluster_area","hetero_radius"];
    vars_n = vars_n_names + "_n";
    vars_radial_density = compose("radial_density_ring_%02d", 1:10);
    vars_radial_heterochromatin_density = compose("radial_heterochromatin_density_ring_%02d", 1:10);
    vars_radial_delta = ["radial_density_ring_gradient_major_axis","radial_density_ring_gradient_minor_axis",...
        "radial_heterochromatin_density_ring_gradient_major_axis","radial_heterochromatin_density_ring_gradient_minor_axis"];
    load(data_info_table{1,"filepath"},'area_thresholds');
    vars_cluster = reshape("voronoi_cluster" +compose("_%03.0fnm_",area_thresholds) + ["area";"density";"gyration_R";"n_locs"]+"_",1,[]);
    vars_cluster = reshape(vars_cluster + ["median";"std";"skewness"],1,[]);
    vars = [vars_string,...
        vars_morphological,...
        vars_value,...
        vars_chromatin,...
        vars_distribution,...
        vars_n,...
        vars_radial_density,...
        vars_radial_heterochromatin_density,...
        vars_radial_delta,...
        vars_cluster];
    vars_to_load = cellstr([vars_value vars_distribution_names vars_n_names,...
        ["area_thresholds"...
        "voronoi_areas","locs_norm","polygon",...
        "x_length","y_length",...
        "eigenvalues",...
        "radial_density","radial_hetero_density",...
        "cluster_area","cluster_density","cluster_gyration_R","cluster_n_locs"]]);
    % group voronoi data
    n_samp = size(data_info_table,1);
    T = table('Size',[n_samp length(vars)], ...
        'VariableTypes',[repmat("string",length(vars_string),1); ...
                         repmat("double",length(vars)-length(vars_string),1)]',...
        'VariableNames',vars);
    for s=1:n_samp
        %% string variables
        T{s,"group"} = data_info_table{s,"group"};
        T{s,"biological_replicate"} = data_info_table{s,"rep"};
        [~,a]= regexp(data_info_table{s,"name"}, ".*-([x0-9]+).*","split","tokens");
        T{s,"name"} = a{1};
        try
            voronoi_data = load_variables(data_info_table{s,"filepath"},vars_to_load);
            %% value variables
            for j=1:length(vars_value)
                if ~isempty(voronoi_data.(vars_value(j)))
                    T{s,vars_value(j)} = voronoi_data.(vars_value(j));
                else
                    T{s,vars_value(j)} = NaN;
                end
            end
            %% morphological variables
            locs_norm = voronoi_data.locs_norm;
            % Generate convex hull
            locs_norm_downsize = locs_norm(1:1000:end,:);
            conv_hull_points = locs_norm_downsize(boundary(locs_norm_downsize,0),:); % Retrieve the actual boundary points.
            polygon_conv_hull_cluster = polyshape(conv_hull_points(1:end-1,1),conv_hull_points(1:end-1,2)); % Create a polygon of the convex hull points (which will be completely filled).
            % Store the shape descriptors: localization descriptors.
            T{s,"area"} = area(voronoi_data.polygon);
            T{s,"major_axis"} = max([voronoi_data.x_length,voronoi_data.y_length]);
            T{s,"minor_axis"} = min([voronoi_data.x_length,voronoi_data.y_length]);
            T{s,"convex_area"} = area(polygon_conv_hull_cluster);
            T{s,"gyration_radius"} = sqrt(mean((locs_norm(:,1)-mean(locs_norm(:,1))).^2 + (locs_norm(:,2)-mean(locs_norm(:,2))).^2));
            eigenvalues = voronoi_data.eigenvalues;
            T{s,"eigenvalues_ratio"} = eigenvalues(1) / eigenvalues(2);
            T{s,"eigenentropy"} = -eigenvalues(1)/sum(eigenvalues)*log(eigenvalues(1)/sum(eigenvalues)) - eigenvalues(2)/sum(eigenvalues)*log(eigenvalues(2)/sum(eigenvalues));
            % Store the shape descriptors: boundary descriptors.
            T{s,"perimeter"} = perimeter(voronoi_data.polygon);
            T{s,"convex_perimeter"} = perimeter(polygon_conv_hull_cluster);
            % Store the dependent shape descriptors (coordinate-based).
            T{s,"aspect_ratio"} = T{s,"major_axis"} / T{s,"minor_axis"};
            T{s,"form_ratio"} = T{s,"area"} / (T{s,"major_axis"}^2);
            T{s,"rectangularity"} = T{s,"area"} / (T{s,"major_axis"}*T{s,"minor_axis"});
            T{s,"circularity"} = 4 * pi * T{s,"area"} / (T{s,"perimeter"}^2);
            T{s,"convex_circularity"} = 4 * pi*T{s,"convex_area"} / (T{s,"convex_perimeter"}^2);
            T{s,"convexity"} = T{s,"perimeter"} / T{s,"convex_perimeter"};
            T{s,"solidity"} = T{s,"area"} / T{s,"convex_area"};
            T{s,"equivalent_diameter"} = sqrt(4 * T{s,"area"} / pi);
            T{s,"fiber_length"} = (T{s,"perimeter"} + realsqrt(abs((T{s,"perimeter"})^2 - 16 * T{s,"area"}))) / 4;
            T{s,"fiber_width"} = T{s,"area"} / T{s,"fiber_length"};
            T{s,"curl"} = T{s,"major_axis"} / T{s,"fiber_length"};
            T{s,"major_axis_norm_area"} = T{s,"major_axis"} ./ T{s,"area"};
            T{s,"minor_axis_norm_area"} = T{s,"minor_axis"} ./ T{s,"area"};
            T{s,"bending_energy_norm_area"} = T{s,"bending_energy"}./ T{s,"area"};            
            %% chromatin variables
            v_density = 1./(voronoi_data.voronoi_areas);
            for j=1:length(vars_chromatin)
                T{s,vars_chromatin(j)} = prctile(v_density,percentile_thresh(j));
            end
            %% distribution variables
            for j=1:length(vars_distribution_names)
                if ~isempty(voronoi_data.(vars_distribution_names(j)))
                    T{s,vars_distribution_names(j)+"_median"} = median(voronoi_data.(vars_distribution_names(j)),"omitmissing");
                    T{s,vars_distribution_names(j)+"_std"} = std(voronoi_data.(vars_distribution_names(j)),"omitmissing");
                    T{s,vars_distribution_names(j)+"_skewness"} = skewness(voronoi_data.(vars_distribution_names(j)),1);
                else
                    T{s,vars_distribution_names(j)+["_median","_std","_skewness"]} = [NaN NaN NaN];
                end
            end
            %% N variables
            for j=1:size(vars_n_names,2)
                if ~isempty(voronoi_data.(vars_n_names(j)))
                    T{s,vars_n(j)} = length(voronoi_data.(vars_n_names(j)));
                else
                    T{s,vars_n(j)} = NaN;
                end
            end
            %% radial features
            for j=1:length(vars_radial_density)
                T{s,vars_radial_density(j)} = voronoi_data.radial_density(j);
            end
            for j=1:length(vars_radial_heterochromatin_density)
                T{s,vars_radial_heterochromatin_density(j)} = voronoi_data.radial_hetero_density(j);
            end
            n_ring = floor(length(vars_radial_density)*0.8);
            T{s,"radial_density_ring_gradient_major_axis"} = mean(voronoi_data.radial_density(1:n_ring))/(0.8*T{s,"major_axis"});
            T{s,"radial_density_ring_gradient_minor_axis"} = mean(voronoi_data.radial_density(1:n_ring))/(0.8*T{s,"minor_axis"});
            T{s,"radial_heterochromatin_density_ring_gradient_major_axis"} = mean(voronoi_data.radial_hetero_density(1:n_ring))/(0.8*T{s,"major_axis"});
            T{s,"radial_heterochromatin_density_ring_gradient_minor_axis"} = mean(voronoi_data.radial_hetero_density(1:n_ring))/(0.8*T{s,"minor_axis"});
            %% cluster features
            metrics = ["area","density","gyration_R","n_locs"];
            for j=1:length(metrics)
                for k=1:5
                    prefix = sprintf("voronoi_cluster_%03.0fnm_%s_",voronoi_data.area_thresholds(k),metrics(j));
                    T{s,prefix+"median"} = median(voronoi_data.("cluster_"+metrics(j)){k},"omitmissing");
                    T{s,prefix+"std"} = std(voronoi_data.("cluster_"+metrics(j)){k},"omitmissing");
                    T{s,prefix+"skewness"} = skewness(voronoi_data.("cluster_"+metrics(j)){k},1);
                end
            end
            clearvars voronoi_data
        catch ME
            disp(data_info_table{s,"filepath"})
            disp(getReport(ME))
            T{s,4:end} = NaN(1,size(T,2)-3);
        end
    end
end
