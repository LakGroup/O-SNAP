% Generate multiple nucleiload_vars
function sim_data = generate_nucleus_set(save_dir,n,options)
    arguments
        save_dir string;
        n double = 1;
        options.major_axis (1,2) double = [2e4,2.5e3];
        options.minor_axis (1,2) double = [1.3e4,2e3];
        options.int_periphery_cutoff (1,2) double = [0.85,0.05];
        options.density_interior (1,2) double = [7.4e-3,2e-3];
        options.density_periphery (1,2) double = [7.5e-3,2.3e-3];
        options.density_hetero_center (1,2) double = [1e-5,3e-6];
        options.density_hetero (1,2) double = [2.3e-2,4.8e-3];
        options.radius_hetero (1,2) double = [30,4.5];
        options.plot_flag logical = false;
        options.n_processes double = 12;
    end
    if options.plot_flag
        figure();
        tiledlayout('flow');
    end
    % reduce broadcast overhead
    major_axis = repmat(options.major_axis,options.n_processes,1);
    major_axis_mean = major_axis(:,1);
    major_axis_std = major_axis(:,2);
    minor_axis = repmat(options.minor_axis,options.n_processes,1);
    minor_axis_mean = minor_axis(:,1);
    minor_axis_std = minor_axis(:,2);
    int_periphery_cutoff = repmat(options.int_periphery_cutoff,options.n_processes,1);
    int_periphery_cutoff_mean = int_periphery_cutoff(:,1);
    int_periphery_cutoff_std = int_periphery_cutoff(:,2);
    density_interior = repmat(options.density_interior,options.n_processes,1);
    density_interior_mean = density_interior(:,1);
    density_interior_std = density_interior(:,2);
    density_periphery = repmat(options.density_periphery,options.n_processes,1);
    density_periphery_mean = density_periphery(:,1);
    density_periphery_std = density_periphery(:,2);
    density_hetero_center = repmat(options.density_hetero_center,options.n_processes,1);
    density_hetero_center_mean = density_hetero_center(:,1);
    density_hetero_center_std = density_hetero_center(:,2);
    density_hetero = repmat(options.density_hetero,options.n_processes,1);
    density_hetero_mean = density_hetero(:,1);
    density_hetero_std = density_hetero(:,2);
    radius_hetero = repmat(options.radius_hetero,options.n_processes,1);
    radius_hetero_mean = radius_hetero(:,1);
    radius_hetero_std = radius_hetero(:,2);
    plot_flag = repmat(options.plot_flag,options.n_processes,1);
    % split data for parallel
    [sim_data, n_processes] = split_data_for_parallel(cell(1,n),...
        options.n_processes, false);
    idx = split_data_for_parallel(num2cell(1:n), options.n_processes,false);
    % create nuclei
    tic
    parfor p=1:n_processes
        for s=1:length(sim_data{p})
            % create random values based on options
            while 1
                major_axis_val = random('Normal',major_axis_mean(p),major_axis_std(p));
                minor_axis_val = random('Normal',minor_axis_mean(p),minor_axis_std(p));
                int_periphery_cutoff_val = random('Normal',int_periphery_cutoff_mean(p),int_periphery_cutoff_std(p));
                density_interior_val = random('Normal',density_interior_mean(p),density_interior_std(p));
                density_periphery_val = random('Normal',density_periphery_mean(p),density_periphery_std(p));
                density_hetero_center_val = random('Normal',density_hetero_center_mean(p),density_hetero_center_std(p));
                density_hetero_val = random('Normal',density_hetero_mean(p),density_hetero_std(p)) - density_interior_mean(p);
                radius_hetero_val = random('Normal',radius_hetero_mean(p),radius_hetero_std(p)); % random avg cluster size per nucleus
                if any([major_axis_val, minor_axis_val,...
                    int_periphery_cutoff_val, density_interior_val,...
                    density_periphery_val, density_hetero_center_val,...
                    density_hetero_val, radius_hetero_val])
                    break
                end
            end
            angle_val = 360*rand();
            if plot_flag(p)
                nexttile
            end
            points = generate_nucleus([major_axis_val minor_axis_val],[0 0],angle_val,int_periphery_cutoff_val,...
            density_interior_val,density_periphery_val,density_hetero_center_val,density_hetero_val,radius_hetero_val,...
            "plot_flag",plot_flag(p));
            sim_data{p}{s}.x_data = points(:,1);
            sim_data{p}{s}.y_data = points(:,2);
            %% save localization data for subsequent analysis
            sim_data{p}{s} = generate_loc_data(sim_data{p}{s}, 1);
            sim_data{p}{s}.name = ['simulated_nucleus-group1-storm-' sprintf('%02.0f',idx{p}{s})];
            save_voronoi_analysis_data(fullfile(save_dir,[sim_data{p}{s}.name '.mat']), sim_data{p}{s});
            % disp("   " + sim_data{p}{s}.name + ": " + string(datetime))
        end
    end
    sim_data = [sim_data{:}];
    toc
end
