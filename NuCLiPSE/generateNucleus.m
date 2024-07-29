% Generate a nucleus
function points = generateNucleus(axesLengths,center,angle,intPeripheryCutoff,...
    densityInterior,densityPeriphery,densityHeteroCenter,densityHetero,radiusHetero,...
    options)
    arguments
        axesLengths (1,2) double
        center (1,2) double
        angle double
        intPeripheryCutoff double
        densityInterior double
        densityPeriphery double
        densityHeteroCenter double
        densityHetero double
        radiusHetero double
        options.plotFlag logical = false
    end
    axesLengthsInterior = axesLengths*sqrt(intPeripheryCutoff);
    pointsInterior = generateEllipse(densityInterior,center,axesLengthsInterior);
    pointsPeriphery =  generateEllipse(densityPeriphery,center,axesLengths,axesLengthsInterior);
    pointsHeteroCenter = generateEllipse(densityHeteroCenter,center,axesLengthsInterior);
    pointsHetero = [];
    for i=1:size(pointsHeteroCenter,1)
        a = 360*rand(); % random angle
        l = random("Normal",radiusHetero,25,1,2); % random ellipse axes
        pointsCluster = generateEllipse(densityHetero,pointsHeteroCenter(i,:),l);
        pointsCluster = rotatePoints(pointsCluster, pointsHeteroCenter(i,:), a);
        pointsHetero = [pointsHetero; pointsCluster];
    end
    pointsInterior=rotatePoints(pointsInterior,center,angle);
    pointsPeriphery=rotatePoints(pointsPeriphery,center,angle);
    pointsHeteroCenter=rotatePoints(pointsHeteroCenter,center,angle);
    if ~isempty(pointsHetero)
        pointsHetero=rotatePoints(pointsHetero,center,angle);
    end
    % plot
    if options.plotFlag && ~isempty(pointsHetero)
        plotNucleus(pointsInterior,pointsPeriphery,pointsHeteroCenter,pointsHetero);
    elseif options.plotFlag
        plotNucleus(pointsInterior,pointsPeriphery,pointsHeteroCenter);
    end
    points = [pointsInterior; pointsPeriphery; pointsHetero];
end

% Plot nucleus
function plotNucleus(pointsInterior,pointsPeriphery,pointsHeteroCenter,pointsHetero)
    plotIdx =  floor(1:length(pointsInterior)/100:length(pointsInterior));
    scatter(pointsInterior(plotIdx,1),pointsInterior(plotIdx,2),".k");
    hold on
     plotIdx =  floor(1:length(pointsPeriphery)/100:length(pointsPeriphery));
    scatter(pointsPeriphery(plotIdx,1),pointsPeriphery(plotIdx,2),"sr","filled");
    hold on
    plotIdx =  floor(1:length(pointsHeteroCenter)/100:length(pointsHeteroCenter));
    scatter(pointsHeteroCenter(plotIdx,1),pointsHeteroCenter(plotIdx,2),"xb");
    axis equal
    xlim([-15000 15000])
    ylim([-15000 15000])
    hold on
    if nargin == 4
        plotIdx =  floor(1:length(pointsHetero)/100:length(pointsHetero));
        scatter(pointsHetero(plotIdx:end,1),pointsHetero(plotIdx:end,2),0.5,".m")
        hold on
    end
end

% Generate points in an ellipse or ring
function points = generateEllipse(density,center,axesLengths1,axesLengths2)
    points = generateBoxPoints(density,center,axesLengths1);
    switch nargin
        case 4
            idx = inEllipse(points,center,axesLengths1,axesLengths2);
        case 3
            idx = inEllipse(points,center,axesLengths1);
        otherwise
            return
    end 
    points = points(idx,:);
end

% Rotate points about a center
function pointsRotated = rotatePoints(points,center,angle)
    c = cosd(angle); s = sind(angle); x = center(1); y = center(2);
    rotateMatrix = [c -s -x*c+y*s+x; s c -x*s-y*c+y; 0 0 1];
    pointsRotated = (rotateMatrix*[points'; ones(1,size(points,1))])';
    pointsRotated = pointsRotated(:,1:2);
end

% Get points randomly located within bounding box
function points = generateBoxPoints(density,center,axesLengths)
    nPoints = floor(density*prod(axesLengths));
    x = axesLengths(1) * rand(nPoints,1) - axesLengths(1)/2 + center(1);
    y = axesLengths(2) * rand(nPoints,1) - axesLengths(2)/2 + center(2);
    points = [x y];
end

% Determine which points are in ellipse
function idx = inEllipse(points, center, axesLengths1, axesLengths2)
    nPoints = size(points,1);
    switch nargin
        case 4
            idx = all([((points(:,1) - repmat(center(:,1),nPoints,1))/(axesLengths2(1)/2)).^2 ...
            + ((points(:,2) - repmat(center(:,2),nPoints,1))/(axesLengths2(2)/2)).^2 >= 1 ...
            ((points(:,1) - repmat(center(:,1),nPoints,1))/(axesLengths1(1)/2)).^2 ...
            + ((points(:,2) - repmat(center(:,2),nPoints,1))/(axesLengths1(2)/2)).^2 <= 1
            ],2);
        case 3
            idx = ((points(:,1) - repmat(center(:,1),nPoints,1))/(axesLengths1(1)/2)).^2 ...
            + ((points(:,2) - repmat(center(:,2),nPoints,1))/(axesLengths1(2)/2)).^2 < 1;
        otherwise
            return
    end
end