clear all

folder_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\';
source_extention = 'CryoGrid\source\';
result_extention = 'CryoGridExamplesMongolia\results\';
forcing_extention = 'CryoGridExamplesMongolia\forcing\Ny_Ålesund\';
addpath(genpath([folder_path source_extention])); % Required to read classes correcly


% =========================
%   start EDIT by user
%
run_name = 'NyÅ_simpleSnow';
run_number = '20190901';

grid.snow_height    = .25;
grid.snow_spacing   = .05;
grid.ground_depth   = 2.5;
grid.ground_spacing = .05;

%
%   end EDIT by user
% =========================

% Result data
load([folder_path result_extention run_name '\' run_name '_' run_number '.mat'])

%% Interpolate relevant values to grid
for i = 1:size(out.STRATIGRAPHY{1, 1},1)
    class_name = class(out.STRATIGRAPHY{1, 1}{i, 1});
    if strcmp(class_name(1:6),'GROUND')
        grid.z_surface = out.STRATIGRAPHY{1, 1}{i, 1}.STATVAR.upperPos;
    end
end

grid.z_base = grid.z_surface - grid.ground_depth;
grid.z_top = grid.z_surface + grid.snow_height;

result.grid = [grid.z_top:-grid.snow_spacing:grid.z_surface (grid.z_surface-grid.ground_spacing):-grid.ground_spacing:grid.z_base]';

for i = 1:size(out.STRATIGRAPHY,2)
    temp.grid = [];
    temp.T = [];
    temp.waterIce = []; 
    temp.saturation = [];
    
    for j = 1:size(out.STRATIGRAPHY{1,i},1)
        T           = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.T;
        waterIce    = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.waterIce;
        layerThick  = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick;
        area        = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.area;
        mineral     = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.mineral;
        organic     = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.organic;
        saturation  = waterIce ./ (layerThick.* area - mineral - organic);
        
        temp.grid       = [temp.grid; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos-cumsum(layerThick) + layerThick./2];
        temp.T          = [temp.T; T(1); T];
        temp.waterIce   = [temp.waterIce; waterIce(1); waterIce];
        temp.saturation = [temp.saturation; saturation(1); saturation];
        
    end

    result.waterIce(:,i)    = interp1(temp.grid,temp.waterIce,result.grid);
    result.saturation(:,i)  = interp1(temp.grid,temp.saturation,result.grid);
    result.T(:,i)           = interp1(temp.grid,temp.T,result.grid);

end

%% Plot water content

figure
pcolor(out.TIMESTAMP,result.grid,result.saturation)
shading flat
hold on
contour(out.TIMESTAMP,result.grid,result.T,[0 0],'k')
datetick
title(run_name)
xlim([out.TIMESTAMP(1) out.TIMESTAMP(end)])
c = colorbar;
ylabel(c,'Saturation [-]')