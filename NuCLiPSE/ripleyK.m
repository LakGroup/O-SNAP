function ripleyK(workDir,groups,reps)
dirPathSplit = regexp(workDir,'\','split');
analysisName = dirPathSplit(end);
nSteps = 10;
nBatch = 10;
maxTen = 10^(ceil(log10(100000))-1);
dist = maxTen/nSteps:maxTen/nSteps:maxTen;
RipleyT = getValidVoronoiAreaData(workDir,groups,reps,{'locs_norm','boundary'});
lastColOld = size(RipleyT,2);
%% calculate K values
for i=1:nSteps
    RipleyT.("K_"+num2str(dist(i))) = zeros(size(RipleyT,1),1);
end
K_exp = cell(1,length(groups));
for g=1:length(groups)
    idx = strcmpi(RipleyT.group,groups(g));
    filepaths = RipleyT{idx,'filepath'};
    % measure observed K
    K_obs_g = arrayfun(@(f) ripleyKObs(f, dist), filepaths,'uni',0);
    RipleyT{idx,lastColOld+1:lastColOld+nSteps} = cell2mat(K_obs_g);
    clearvars K_g
    % measure expected K
    K_exp{g} = cell2mat(arrayfun(@(f) ripleyKExp(f, dist,nBatch), filepaths,'uni',0));
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
K_obs_mean = zeros(length(groups),nSteps);
K_obs_std = zeros(length(groups),nSteps);
for g=1:length(groups)
    idx = strcmpi(RipleyT.group,groups(g));
    K_obs_g = RipleyT{idx,lastColOld+1:lastColOld+nSteps};
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
savefig(fullfile(workDir,"RipleyK.fig"));
save(fullfile(workDir,analysisName+"_RipleyK.mat"),"RipleyT");
end
function K = ripleyKCalc(points,dist)
    [minVal,maxVal] = bounds(points);
    K = arrayfun(@(d) RipleysK(points,d,[minVal(1) maxVal(1) minVal(2) maxVal(2)],false), dist,'uni',1);
end
function K = ripleyKObs(filepath, dist)
    data = loadVars(filepath,{'locs_norm'});
    K = ripleyKCalc(data.locs_norm(1:500:end,:),dist);
end
% creates random set of XY points in defined boundary, where X is stored in
% odd columns while Y is in even columns
function K = ripleyKExp(filepath,dist,nBatch)
    K = zeros(nBatch,length(dist));
    % set up
    data = loadVars(filepath,{'boundary','locs_number'});
    nPoints = floor(data.locs_number/500);
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
    for i=1:nBatch
        R = rand(nPoints,1);
        tind = discretize(R,cumsum([0;areas]));
        v1 = XY(tri(tind,1),:);
        v2 = XY(tri(tind,2),:);
        v3 = XY(tri(tind,3),:);
        R1 = rand(nPoints,1);
        points_batch = v1.*repmat(R1,[1 2]) + v2.*repmat(1-R1,[1 2]);
        R2 = sqrt(rand(nPoints,1));
        points_batch = points_batch.*repmat(R2,[1 2]) + v3.*repmat(1-R2,[1 2]);
        K(i,:) = ripleyKCalc(points_batch,dist);
    end
end