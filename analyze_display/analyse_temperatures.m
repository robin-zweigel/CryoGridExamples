% Plot and compare temperatures for two runs
% R. B. Zweigel, January 2022
clear all
% close all

result_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGridExamples\results\';
source_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGrid\source\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
run1 = 'Terelj_NS_CLM5_newForcing';
run2 = 'Terelj_NS_CLM5_newForcing2';

grid_vec = [.05 3.5]; % specify grid to interpolate to - [Spacing Depth]

variable = 'T';
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
temp1 = read_aggregate_layer([result_path run1], variable);
temp2 = read_aggregate_layer([result_path run2], variable);

% Define grid
grid = 0:grid_vec(1):grid_vec(2);

% Ground
depth1 = [zeros(1,size(temp1.GROUND.layerThick,2)); cumsum(temp1.GROUND.layerThick)];
midpoints1 = (depth1(1:end-1,:)+depth1(2:end,:))./2;
result1 = interp2(temp1.time,midpoints1(:,1),temp1.GROUND.(variable),temp1.time,grid');
result1(1,:) = result1(2,:); % assign values to top row (depth = 0) 

depth2 = [zeros(1,size(temp2.GROUND.layerThick,2)); cumsum(temp2.GROUND.layerThick)];
midpoints2 = (depth2(1:end-1,:)+depth2(2:end,:))./2;
result2 = interp2(temp2.time,midpoints2(:,1),temp2.GROUND.(variable),temp2.time,grid');
result2(1,:) = result2(2,:); % assign values to top row (depth = 0) 

% Snow

% Vegetation
if size(temp1.VEGETATION.layerThick) > 0
    Tv1 = temp1.VEGETATION.T;
end
if size(temp2.VEGETATION.layerThick) > 0
    Tv2 = temp2.VEGETATION.T;
end

%% 
figure

if exist('Tv1')
    subplot(3,2,1)
    hold on
    plot(temp1.time,Tv1)
    plot([temp1.time(1) temp1.time(end)],[0 0],'k--')
    datetick
    xlim([temp2.time(1) temp2.time(end)])
    ylim([-50 70])
    title('Canopy Temperature')
end

if exist('Tv2')
    subplot(3,2,2)
    hold on
    plot(temp2.time,Tv2)
    plot([temp2.time(1) temp2.time(end)],[0 0],'k--')
    datetick
    xlim([temp2.time(1) temp2.time(end)])
    title('Canopy Temperature')
    ylim([-50 70])
end

subplot(3,2,3)
plot(temp1.time,result1(1,:))
datetick
xlim([temp2.time(1) temp2.time(end)])
ylim([-30 40])
title('Ground Surface Temperature')

subplot(3,2,4)
plot(temp2.time,result2(1,:))
datetick
xlim([temp2.time(1) temp2.time(end)])
ylim([-30 40])
title('Ground Surface Temperature')

subplot(3,2,5)
plotTemp(run1,temp2.time(1),temp2.time(end), temp1.time,result1, -grid);

subplot(3,2,6)
plotTemp(run2,temp2.time(1),temp2.time(end), temp2.time,result2, -grid);
