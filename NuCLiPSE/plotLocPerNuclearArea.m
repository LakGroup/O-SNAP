function plotLocPerNuclearArea(workDir, groups, voronoiData)
n_groups = length(groups);
[~, groupIdx] = sortSampleToGroups(groups,voronoiData);
%% get data
nucleusAreas = cell(1,n_groups);
locPerArea = cell(1,n_groups);
for g=1:n_groups
    nSampInGroup = length(groupIdx{g});
    if nSampInGroup > 0
        nucleusAreas{g} =  cellfun(@(x) x.nucleus_area*10^-6, voronoiData(groupIdx{g}),'uni',1);
        locPerArea{g} = cellfun(@(x) x.nucleus_totalLocs./(x.nucleus_area*10^-6), voronoiData(groupIdx{g}),'uni',1);
    end
end
%% process data
validGroups = logical(cellfun(@(x) length(x), groupIdx,'uni',1));
x = categorical(groups(validGroups));
y = (cellfun(@(x) mean(x), locPerArea(validGroups),'uni',1))'; % in μm
y_err = (cellfun(@(x) std(x), locPerArea(validGroups),'uni',1))'; % in μm
z = (cellfun(@(x) mean(x), nucleusAreas(validGroups),'uni',1))'; % in μm
z_err = (cellfun(@(x) std(x), nucleusAreas(validGroups),'uni',1))'; % in μm
nil = zeros(length(y),1);
%% plot
figure('visible','off');
title({'Localizations per nuclear area','Nuclear area'})
b = bar(x, [y nil], 'grouped');
ylabel('Localizations per Nuclear Area $(\textrm{loc}/{\mu m}^{-2})$','interpreter','latex');
hold on
yyaxis right
bar(x, [nil z], 'grouped');
ylabel('Nuclear Area $({\mu m}^{2})$','interpreter','latex');
hold on
x_pos = zeros(2,length(b(1).XEndPoints));
for j = 1:length(b(1).XEndPoints)
    x_pos(1,j) = sort(b(1).XEndPoints(j));
    x_pos(2,j) = sort(b(2).XEndPoints(j));
end
yyaxis left
errorbar(x_pos(1,:),y,y_err,'k','linestyle','none');
ylim([0 14000])
yyaxis right
errorbar(x_pos(2,:), z, z_err,'k','linestyle','none');
ylim([0 300])
saveas(gcf,fullfile(workDir, 'LocPerArea_NuclearArea.png'));
end