clear all

meta.folder_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\';
meta.source_extention = 'CryoGrid\source\';
meta.result_extention = 'CryoGridExamplesMongolia\results\';
% forcing_extention = 'CryoGridExamplesMongolia\forcing\Terelj';
addpath(genpath([meta.folder_path meta.source_extention])); % Required to read classes correcly


% ------- EDIT BY USER ----------
meta.run_name = 'Terelj_test';
meta.run_number = '20160901';

grid.snow_height    = 0;.25;
grid.snow_spacing   = .025;
grid.ground_depth   = 2;
grid.ground_spacing = .05;

depths = [grid.snow_height:-grid.snow_spacing:0 -grid.ground_spacing:-grid.ground_spacing:-grid.ground_depth]';

%% Plot temperatures
result= read_ground_temperatures(meta, depths);

figure
pcolor(result.time,result.depths,result.T)
title(meta.run_name)
shading flat
hold on
contour(result.time,result.depths,result.T,[0 0],'k')
datetick
c = colorbar;
ylabel(c,'Temperature [{\circ}C]')
xlim([result.time(1) result.time(end)])
caxis([-41 35])