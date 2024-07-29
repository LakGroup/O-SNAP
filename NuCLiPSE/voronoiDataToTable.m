function T = voronoiDataToTable(dataInfoTable)
    % define variables of interest
    varsString = {'Group','BioRep','name'};
    varsString_N = size(varsString,2);
    varsDouble = {'nucleus_radius','lads2total',...
        'locs_density','locs_number','ElasticEnergy','BendingEnergy',...
        'peripheryDensity','interiorDensity','peripheryHeteroDensity','interiorHeteroDensity'};
    varsDouble_N = size(varsDouble,2);
    varsMedian = {'hetero_radius','hetero_density','hetero_nLocs',...
        'lads_area','lads_density','lads_nLocs',...
        'lads_seg_len','lads_seg_thickness',...
        'reduced_log_voronoi_density','spacing','BorderCurvature'};
    varsMedian_N = size(varsMedian,2);
    varsN = {'lads','non_lads'};
    varsN_N = size(varsN,2);
    varsGeometric = {'Area','MajorAxis','MinorAxis','ConvexArea',...
    'GyrationRadius','EigenvaluesRatio','Eigenentropy',...
    'Perimeter','ConvexPerimeter',...
    'Density','AspectRatio','FormRatio','Rectangularity','Circularity','ConvexCircularity',...
    'Convexity','Solidity','EquivalentDiameter','FiberLength','FiberWidth','Curl',...
    'MajorAxisNormArea','MinorAxisNormArea','BendingEnergyNormArea'};
    varsGeometric_N = size(varsGeometric,2);
    radDenCol = arrayfun(@(x) sprintf('radialDensity%02d',x), 1:10,'uni',0);
    radHetDenCol = arrayfun(@(x) sprintf('radialHeteroDensity%02d',x), 1:10,'uni',0);
    clusterAreaCol = arrayfun(@(x) sprintf('clusterArea%02d',x), 1:5,'uni',0);
    clusterDensityCol = arrayfun(@(x) sprintf('clusterDensity%02d',x), 1:5,'uni',0);
    clusterGyrationRCol = arrayfun(@(x) sprintf('clusterGyrationR%02d',x), 1:5,'uni',0);
    clusterNLocsCol = arrayfun(@(x) sprintf('clusterNLocs%02d',x), 1:5,'uni',0);
    varsChromatin = [{'voronoiDensity40','voronoiDensity70'},...
        radDenCol,radHetDenCol,...
        clusterAreaCol,clusterDensityCol,clusterGyrationRCol,clusterNLocsCol,...
        {'DeltaRadialDensityMajor','DeltaHeteroRadialDensityMajor',...
        'DeltaRadialDensityMinor','DeltaHeteroRadialDensityMinor'}];
    varsChromatin_N = size(varsChromatin,2);
    varsToLoad = [varsDouble varsMedian varsN,...
        {'name',...
        'voronoi_areas','locs_norm','polygon',...
        'x_length','y_length',...
        'Eigenvalues',...
        'radialDensity','radialHeteroDensity',...
        'clusterArea','clusterDensity','clusterGyrationR','clusterNLocs'}];
    vars = cellfun(@(x) replace(x,'_',''), [varsString ...
        varsDouble ...
        cellfun(@(x) [x '_median'], varsMedian,'uni',0) ...
        cellfun(@(x) [x '_N'], varsN,'uni',0),...
        varsGeometric, ...
        varsChromatin ...
        ],'uni',0);
    vars_N = size(vars,2);
    % % ignore samples that aren't valid
    % idx = validVarIdx(voronoiData,[varsDouble varsMedian varsN]);
    % voronoiData = voronoiData(all(idx,1));
    % group voronoi data
    nSamp = size(dataInfoTable,1);
    T = table('size',[nSamp vars_N], ...
    'VariableTypes',[repmat("string",varsString_N,1); repmat("double",varsDouble_N+varsMedian_N+varsN_N+varsGeometric_N+varsChromatin_N,1)],...
    'VariableNames',vars);
    % w = waitbar(0, ['Loading file: 0/' num2str(nSamp)]);
    for s=1:nSamp
        % waitbar(s/nSamp, w, ['Loading file: ' num2str(s) '/' num2str(nSamp)]);
        voronoiData = load(fullfile(dataInfoTable{s,"filepath"}),varsToLoad{:});
        T{s,"Group"} = dataInfoTable{s,"group"};
        T{s,"BioRep"} = dataInfoTable{s,"rep"};
        [~,a]= regexp(voronoiData.name, '.*-([x0-9]+).*','split','tokens');
        T{s,"name"} = a{1};
        for j=1:varsDouble_N
            if ~isempty(voronoiData.(varsDouble{j}))
                T{s,varsString_N+j} = voronoiData.(varsDouble{j});
            else
                T{s,varsString_N+j} = NaN;
            end
        end
        for j=1:varsMedian_N
            if ~isempty(voronoiData.(varsMedian{j}))
                T{s,varsString_N+varsDouble_N+j} = median(voronoiData.(varsMedian{j}),"omitmissing");
            else
                T{s,varsString_N+j} = NaN;
            end
        end
        for j=1:varsN_N
            if ~isempty(voronoiData.(varsN{j}))
                T{s,varsString_N+varsDouble_N+varsMedian_N+j} = max(voronoiData.(varsN{j})(:,3));
            else
                T{s,varsString_N+j} = NaN;
            end
        end
        v_density = 1./(voronoiData.voronoi_areas);
        T{s,"voronoiDensity40"} = prctile(v_density,40);
        T{s,"voronoiDensity70"} = prctile(v_density,70);
        locs_norm = voronoiData.locs_norm;
        % Generate convex hull
        locs_norm_downsize = locs_norm(1:1000:end,:);
        ConvHullPoints = locs_norm_downsize(boundary(locs_norm_downsize,0),:); % Retrieve the actual boundary points.
        PolygonConvHullCluster = polyshape(ConvHullPoints(1:end-1,1),ConvHullPoints(1:end-1,2)); % Create a polygon of the convex hull points (which will be completely filled).
        % Store the shape descriptors: localization descriptors.
        T{s,"Area"} = area(voronoiData.polygon);
        T{s,"MajorAxis"} = max([voronoiData.x_length,voronoiData.y_length]);
        T{s,"MinorAxis"} = min([voronoiData.x_length,voronoiData.y_length]);
        T{s,"ConvexArea"} = area(PolygonConvHullCluster);
        T{s,"GyrationRadius"} = sqrt(mean((locs_norm(:,1)-mean(locs_norm(:,1))).^2 + (locs_norm(:,2)-mean(locs_norm(:,2))).^2));
        Eigenvalues = voronoiData.Eigenvalues;
        T{s,"EigenvaluesRatio"} = Eigenvalues(1) / Eigenvalues(2);
        T{s,"Eigenentropy"} = -Eigenvalues(1)/sum(Eigenvalues)*log(Eigenvalues(1)/sum(Eigenvalues)) - Eigenvalues(2)/sum(Eigenvalues)*log(Eigenvalues(2)/sum(Eigenvalues));
        % Store the shape descriptors: boundary descriptors.
        T{s,"Perimeter"} = perimeter(voronoiData.polygon);
        T{s,"ConvexPerimeter"} = perimeter(PolygonConvHullCluster);
        % Store the dependent shape descriptors (coordinate-based).
        T{s,"Density"} = T{s,"locsnumber"} / T{s,"Area"};
        T{s,"AspectRatio"} = T{s,"MajorAxis"} / T{s,"MinorAxis"};
        T{s,"FormRatio"} = T{s,"Area"} / (T{s,"MajorAxis"}^2);
        T{s,"Rectangularity"} = T{s,"Area"} / (T{s,"MajorAxis"}*T{s,"MinorAxis"});
        T{s,"Circularity"} = 4 * pi * T{s,"Area"} / (T{s,"Perimeter"}^2);
        T{s,"ConvexCircularity"} = 4 * pi*T{s,"ConvexArea"} / (T{s,"ConvexPerimeter"}^2);
        T{s,"Convexity"} = T{s,"Perimeter"} / T{s,"ConvexPerimeter"};
        T{s,"Solidity"} = T{s,"Area"} / T{s,"ConvexArea"};
        T{s,"EquivalentDiameter"} = sqrt(4 * T{s,"Area"} / pi);
        T{s,"FiberLength"} = (T{s,"Perimeter"} + realsqrt(abs((T{s,"Perimeter"})^2 - 16 * T{s,"Area"}))) / 4;
        T{s,"FiberWidth"} = T{s,"Area"} / T{s,"FiberLength"};
        T{s,"Curl"} = T{s,"MajorAxis"} / T{s,"FiberLength"};
        T{s,"MajorAxisNormArea"} = T{s,"MajorAxis"} / T{s,"Area"};
        T{s,"MinorAxisNormArea"} = T{s,"MinorAxis"} / T{s,"Area"};
        T{s,"BendingEnergyNormArea"} = T{s,"BendingEnergy"} / T{s,"Area"};
        for j=1:length(radDenCol)
            T{s,radDenCol{j}} = voronoiData.radialDensity(j);
            T{s,radHetDenCol{j}} = voronoiData.radialHeteroDensity(j);
        end
        nRing = floor(length(radDenCol)*0.8);
        T{s,"DeltaRadialDensityMajor"} = mean(voronoiData.radialDensity(1:nRing))/(0.8*T{s,"MajorAxis"});
        T{s,"DeltaRadialDensityMinor"} = mean(voronoiData.radialDensity(1:nRing))/(0.8*T{s,"MinorAxis"});
        T{s,"DeltaHeteroRadialDensityMajor"} = mean(voronoiData.radialHeteroDensity(1:nRing))/(0.8*T{s,"MajorAxis"});
        T{s,"DeltaHeteroRadialDensityMinor"} = mean(voronoiData.radialHeteroDensity(1:nRing))/(0.8*T{s,"MinorAxis"});
        for j=1:length(clusterAreaCol)
            T{s,clusterAreaCol{j}} = median(voronoiData.clusterArea{j},"omitmissing");
            T{s,clusterDensityCol{j}} = median(voronoiData.clusterDensity{j},"omitmissing");
            T{s,clusterGyrationRCol{j}} = median(voronoiData.clusterGyrationR{j},"omitmissing");
            T{s,clusterNLocsCol{j}} = median(voronoiData.clusterNLocs{j},"omitmissing");
        end
        clearvars voronoiData
    end
    % close(w)
end
% 
% function idx = validVarIdx(dataStruct,vars)
%     idx =  cell2mat(cellfun(@(j) cellfun(@(i) any(strcmp(fieldnames(i),j)), dataStruct,'uni',1), vars, 'uni',0)');
% end