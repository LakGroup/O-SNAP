
% n_r = size(T_radialDensity,2);
% c = jet(n_r);
% 
% figure(); 
% h = cell(1,n_r);
% for i=1:n_r
%     h{i} = histogram((T_radialDensity{:,i}),'FaceColor',c(i,:),'BinLimits',[0 maxDensity],'NumBins',50);
%     hold on
%     if i>1
%         for j=1:i-1
%             h{j}.FaceAlpha = h{j}.FaceAlpha * 0.7;
%             h{j}.EdgeAlpha = h{j}.EdgeAlpha * 0.7;
%         end
%     end
%     xlim([0 maxDensity])
%     ylim([0 180])
%     pause(1/n_r)
% end
% for j=1:n_r
%     h{j}.FaceAlpha = 1;
%     h{j}.EdgeAlpha = 1;
% end
% pause(0.2)
% 
% figure();
% for i=1:length(groups)
%     plot(1:n_r,median(X(strcmpi(T.Group,groups{i}),:)));
%     hold on
% end
% xlabel('Ring N')
% ylabel('Density');
% legend(groups);
% 
% figure();
% for i=1:length(groups)
%     plot(1:n_r,median(Y(strcmpi(T.Group,groups{i}),:)));
%     hold on
% end
% xlabel('Ring N')
% ylabel('Density');
% legend(groups);