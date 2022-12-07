% Plot and compare daily temperatures simulated temperatures to observations
% R. B. Zweigel, September 2022
clear all
% close all

result_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGridExamples\results\';
source_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGrid\source\';
obs_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\PhD stuff\Data\Field data\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
run = 'Terelj_NS_CLM5_newForcing2';

variable = 'T';

grid_vec = [.05 3.5]; % specify grid to interpolate to - [Spacing Depth]

obs_filename = 'FA_daily_soil.mat';
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
temp = read_aggregate_layer([result_path run], variable);

% Observations
load([obs_path obs_filename])

% Interpolate to obs. depths
depth = [zeros(1,size(temp.GROUND.layerThick,2)); cumsum(temp.GROUND.layerThick)];
midpoints = (depth(1:end-1,:)+depth(2:end,:))./2;
temp2 = interp2(temp.time,midpoints(:,1),temp.GROUND.(variable),temp.time,(data.z_Tg')./100);
temp2(data.z_Tg == 0,:) = interp2(temp.time,midpoints(:,1),temp.GROUND.(variable),temp.time,.05);
result = squeeze(mean(reshape(temp2,size(temp2,1),[],size(temp2,2)/4),2));
result_time = temp.time(1):temp.time(end);


%% Plotting
    
figure

for i = 1:6
    subplot(6,1,i)
    plot(data.time,data.Tg(:,i))
    hold on
    plot(result_time,result(i,:))
    datetick
    title(['Ground temperatures at ' num2str(data.z_Tg(i)) ' cm'])
    xlim([result_time(1) result_time(end)])
end
legend('Observations','Simulations')

