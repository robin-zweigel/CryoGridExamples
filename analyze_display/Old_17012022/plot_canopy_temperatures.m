% Plot vegetation temperatures for VEGETATION_simpleShading and
% VEGETATION_CLM5 classes

% Robin B. Zweigel, October 2021

clear all

meta.folder_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\';
meta.source_extention = 'CryoGrid\source\';
meta.result_extention = 'CryoGridExamplesMongolia\results\';
meta.forcing_extention = 'CryoGridExamplesMongolia\forcing\Terelj';
addpath(genpath([meta.folder_path meta.source_extention])); % Required to read classes correcly

run_name = {'Terelj_test'};
meta.run_number = '20150901';

%%
figure
hold on
for n = 1:length(run_name)
    meta.run_name = run_name{n};
    result = read_canopy_temperatures(meta);
    plot(result.time,result.Tv)
    
    drawnow
    disp(['MACT = ' num2str(mean(result.Tv(1,:))) '{\circ}C'])
end

datetick

title('Canopy temperature')
ylabel('Temperature [{\circ}C]')
legend(run_name)
xlim([result.time(1) result.time(end)])