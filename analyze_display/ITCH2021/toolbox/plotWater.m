function plotWater(name,startdate,enddate, time, water, grid)
% Function plotting Water results from CryoGrid simulations.
% Cas Renette, Léo Martin, Oslo University 2021.

load('toolbox\colormaps\oslo.mat'); %load colormap
 
cmap=oslo; %select colormap
set(imagesc(time, grid, water, 'AlphaData',~isnan(water))); 
colormap(flipud(cmap));
axis xy;
datetick;
title([name, ': Water'],'Interpreter','none');
a=colorbar;
ylabel(a,'volumetric water content','FontSize',10,'Rotation',90,'FontWeight','bold');
caxis([0 0.8]);
xlim([startdate,enddate]);
xlabel('Year','FontWeight','bold');
ylabel('Elevation [m.a.s.l.]','FontWeight','bold');

end