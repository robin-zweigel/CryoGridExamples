function plotTemp(name,startdate,enddate, time,Temp, grid)
% Function plotting temperature results from CryoGrid simulations.
% Cas Renette, Léo Martin, Oslo University 2021.

load('toolbox\colormaps\vik.mat'); %load oleron scientific colormap

cmap=vik; %select colormap
set(imagesc(time, grid, Temp, 'AlphaData',~isnan(Temp)));
% colormap(ax2,cmap);
colormap(cmap);
axis xy;
datetick;
title([name, ': T'],'Interpreter','none');
a=colorbar;
ylabel(a,'T[^{\circ}C]','FontSize',10,'Rotation',90,'FontWeight','bold');
caxis([-50, 50]);
cRange = caxis;
xlim([startdate,enddate]);
hold on;
contour(time,grid,Temp,[-50 0 50],'LineColor','k','LineWidth',1.25);
caxis(cRange);
%annotation('textbox',[0.15, 0.15, 0.2, 0.2],'String',str1,'FitBoxToText','on');
xlabel('Year','FontWeight','bold');
ylabel('Elevation [m.a.s.l.]','FontWeight','bold');

end

