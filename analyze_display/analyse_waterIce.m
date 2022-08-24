% Plot and compare water/ice for two runs
% R. B. Zweigel, January 2022
clear all
% close all

result_path = 'C:\Users\robinbz\Documents\GitHub\CryoGridExamplesMongolia\results\';
source_path = 'C:\Users\robinbz\Documents\GitHub\CryoGrid\source\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
run1 = 'Terelj_flat_CLM5_snow';
run2 = 'Terelj_flat_CLM5';

grid_vec = [.025 2]; % specify grid to interpolate to - [Spacing Depth]

variables = {'water', 'waterIce', 'ice', 'energy'};
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
data1 = read_aggregate_layer([result_path run1], variables);
data2 = read_aggregate_layer([result_path run2], variables);

% Define grid
grid = 0:grid_vec(1):grid_vec(2);

% Ground
depth1 = [zeros(1,size(data1.GROUND.layerThick,2)); cumsum(data1.GROUND.layerThick)];
midpoints1 = (depth1(1:end-1,:)+depth1(2:end,:))./2;
for var_i = 1:length(variables)
    temp1 = data1.GROUND.(string(variables(var_i))) ./ (data1.GROUND.area.*data1.GROUND.layerThick);
    result1.(string(variables(var_i))) = interp2(data1.time,midpoints1(:,1),temp1,data1.time,grid');
    result1.(string(variables(var_i)))(1,:) = result1.(string(variables(var_i)))(2,:); % assign values to top row (depth = 0)
end

depth2 = [zeros(1,size(data2.GROUND.layerThick,2)); cumsum(data2.GROUND.layerThick)];
midpoints2 = (depth2(1:end-1,:)+depth2(2:end,:))./2;
for var_i = 1:length(variables)
    temp2 = data2.GROUND.(string(variables(var_i))) ./ (data2.GROUND.area.*data2.GROUND.layerThick);
    result2.(string(variables(var_i))) = interp2(data2.time,midpoints1(:,1),temp2,data2.time,grid');
    result2.(string(variables(var_i)))(1,:) = result2.(string(variables(var_i)))(2,:); % assign values to top row (depth = 0)
end

%% Plot the waterIce

figure
subplot(2,2,1)
plotWaterIce(run1, data1.time(1), data1.time(end), data1.time, result1.waterIce, -grid)

subplot(2,2,3)
plotWater(run1, data1.time(1), data1.time(end), data1.time, result1.water, -grid)

subplot(2,2,2)
plotWaterIce(run2, data1.time(1), data1.time(end), data1.time, result2.waterIce, -grid)

subplot(2,2,4)
plotWater(run2, data1.time(1), data1.time(end), data1.time, result2.water, -grid)
