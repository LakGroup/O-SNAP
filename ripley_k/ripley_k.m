function ripley_k(work_dir,groups,reps)
dir_path_split = regexp(work_dir,'\','split');
analysis_name = dir_path_split(end);
n_steps = 10;
n_batch = 10;
max_ten = 10^(ceil(log10(100000))-1);
dist = max_ten/n_steps:max_ten/n_steps:max_ten;
Ripley_t = get_valid_voronoi_data(work_dir,groups,reps,{'locs_norm','boundary'});
last_col_old = size(Ripley_t,2);
%% calculate K values
for i=1:n_steps
    Ripley_t.("K_"+num2str(dist(i))) = zeros(size(Ripley_t,1),1);
end
K_exp = cell(1,length(groups));
for g=1:length(groups)
    idx = strcmpi(Ripley_t.group,groups(g));
    filepaths = Ripley_t{idx,'filepath'};
    % measure observed K
    K_obs_g = arrayfun(@(f) ripley_KObs(f, dist), filepaths,'uni',0);
    Ripley_t{idx,last_col_old+1:last_col_old+n_steps} = cell2mat(K_obs_g);
    clearvars K_g
    % measure expected K
    K_exp{g} = cell2mat(arrayfun(@(f) ripley_KExp(f, dist,n_batch), filepaths,'uni',0));
end
K_exp = vertcat(K_exp{:});
%% plot
figure();
% expected
K_exp_mean = mean(K_exp);
K_exp_std = std(K_exp);
p = plot(dist,K_exp_mean,'k--');
hold on
c = fill([dist dist(end:-1:1)], [K_exp_mean+K_exp_std K_exp_mean(end:-1:1)-K_exp_std(end:-1:1)],p.Color);
c.FaceAlpha = 0.2;
c.EdgeColor = "none";
hold on
% observed
K_obs_mean = zeros(length(groups),n_steps);
K_obs_std = zeros(length(groups),n_steps);
for g=1:length(groups)
    idx = strcmpi(Ripley_t.group,groups(g));
    K_obs_g = Ripley_t{idx,last_col_old+1:last_col_old+n_steps};
    K_obs_mean(g,:) = mean(K_obs_g);
    K_obs_std(g,:) = std(K_obs_g);
    p = plot(dist,K_obs_mean(g,:));
    hold on
    % plot observed confidence interval
    c = fill([dist dist(end:-1:1)],...
        [K_obs_mean(g,:)+K_obs_std(g,:) K_obs_mean(g,end:-1:1)-K_obs_std(g,end:-1:1)],p.Color);
    c.FaceAlpha = 0.2;
    c.EdgeColor = "none";
    hold on
end
title('Ripleys K');
labels = [groups;repelem({''},1,length(groups))];
labels = labels(:)';
legend([{'Random pattern',''} labels]);
savefig(fullfile(work_dir,"ripley_K.fig"));
save(fullfile(work_dir,analysis_name+"_ripley_K.mat"),"Ripley_t");
end
function K = ripley_K_calc(points,dist)
    [min_val,max_val] = bounds(points);
    K = arrayfun(@(d) Ripleys_k(points,d,[min_val(1) max_val(1) min_val(2) max_val(2)],false), dist,'uni',1);
end
function K = ripley_KObs(filepath, dist)
    data = load_variables(filepath,{'locs_norm'});
    K = ripley_K_calc(data.locs_norm(1:500:end,:),dist);
end
% creates random set of XY points in defined boundary, where X is stored in
% odd columns while Y is in even columns
function K = ripley_KExp(filepath,dist,n_batch)
    K = zeros(n_batch,length(dist));
    % set up
    data = load_variables(filepath,{'boundary','locs_number'});
    n_points = floor(data.locs_number/500);
    % create triangulation of polygon
    XY = [data.boundary];
    tri = delaunay(XY);
    % use a 2x2 determinant, in a vectorized form
    v1 = XY(tri(:,1),:);
    v2 = XY(tri(:,2),:);
    v3 = XY(tri(:,3),:);
    % translate
    v1 = v1-v3;
    v2 = v2-v3;
    % vectorized determinant
    % divide by factorial(2) for the area
    areas = (v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;
    % normalize the areas to sum to 1
    areas = areas/sum(areas);
    % genereate random points
    for i=1:n_batch
        R = rand(n_points,1);
        tind = discretize(R,cumsum([0;areas]));
        v1 = XY(tri(tind,1),:);
        v2 = XY(tri(tind,2),:);
        v3 = XY(tri(tind,3),:);
        R1 = rand(n_points,1);
        points_batch = v1.*repmat(R1,[1 2]) + v2.*repmat(1-R1,[1 2]);
        R2 = sqrt(rand(n_points,1));
        points_batch = points_batch.*repmat(R2,[1 2]) + v3.*repmat(1-R2,[1 2]);
        K(i,:) = ripley_K_calc(points_batch,dist);
    end
end

