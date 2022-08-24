% Plot and compare surface energy balance for two runs
% R. B. Zweigel, January 2022
clear all
% close all

result_path = 'C:\Users\robinbz\Documents\GitHub\CryoGridExamples\results\';
source_path = 'C:\Users\robinbz\Documents\GitHub\CryoGrid\source\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
% run1 = 'Terelj_NS_CLM5_snow';
run2 = 'Terelj_NS_CLM5_base';

variables = {'Lin','Lout','Sin','Sout','Qe','Qh'};
%
% == == == END EDIT == == == ==

% data1 = read_aggregate_surface([result_path run1], variables);
data2 = read_aggregate_surface([result_path run2], variables);

%% Plot the results

figure
subplot(2,2,1)
hold on
plot(data1.time,data1.Lin(1,:)-data1.Lout(1,:))
plot(data2.time,data2.Lin(1,:)-data2.Lout(1,:),'--')
datetick
xlim([data1.time(1) data1.time(end)])
title('Absorbed Longwave (Lin-Lout)')
ylabel('Radiation [W/m^2]')
legend({string(run1),string(run2)},'Location','southwest')

subplot(2,2,2)
hold on
plot(data1.time,data1.Sin(1,:)-data1.Sout(1,:))
plot(data2.time,data2.Sin(1,:)-data2.Sout(1,:),'--')
datetick
xlim([data1.time(1) data1.time(end)])
title('Absorbed Shortwave (Sin-Sout)')
ylabel('Radiation [W/m^2]')

subplot(2,2,3)
hold on
plot(data1.time,sum(data1.Qe))
plot(data2.time,sum(data2.Qe),'--')
datetick
xlim([data1.time(1) data1.time(end)])
title('Latent energy (Qe)')
ylabel('Flux [W/m^2]')

subplot(2,2,4)
hold on
plot(data1.time,sum(data1.Qh))
plot(data2.time,sum(data2.Qh),'--')
datetick
xlim([data1.time(1) data1.time(end)])
title('Sensible heat (Qe)')
ylabel('Flux [W/m^2]')