function plotIce(name,startdate,enddate, time, ice, grid)
% Function plotting ice results from CryoGrid simulations.
% Cas Renette, Léo Martin, Oslo University 2021.

load('toolbox\colormaps\oslo.mat'); %load colormap

cmap=oslo; %select colormap
set(imagesc(time, grid, ice, 'AlphaData',~isnan(ice))); 
colormap(flipud(cmap));
axis xy;
datetick;
title([name, ': Ice'],'Interpreter','none');
a=colorbar;
ylabel(a,'volumetric ice content','FontSize',10,'Rotation',90,'FontWeight','bold');
caxis([0 0.8]);
xlim([startdate,enddate]);
xlabel('Year','FontWeight','bold');
ylabel('Elevation [m.a.s.l.]','FontWeight','bold');

end