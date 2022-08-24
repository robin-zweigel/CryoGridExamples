function plotAL(tvec,depthofFrozenGroud)
% Function to plot the yearly active layer based on the depth of frozen
% ground. 
% Cas Renette, Léo Martin, Oslo University 2021.

[t_y,d_y]=timeAverages(tvec,depthofFrozenGroud,'yearly', 'max');
plot(t_y,d_y,'k','LineWidth',2)
hold on
set(gca, 'YDir','reverse')
plot(t_y,d_y,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
title('Active layer trend')
ylabel('Active layer (m)','FontWeight','bold')
xlabel('Year','FontWeight','bold');
datetick
hold off

end

