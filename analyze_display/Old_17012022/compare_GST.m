clear all

meta.folder_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\';
meta.source_extention = 'CryoGrid\source\';
meta.result_extention = 'CryoGridExamplesMongolia\results\';
% forcing_extention = 'CryoGridExamplesMongolia\forcing\Terelj';
addpath(genpath([meta.folder_path meta.source_extention])); % Required to read classes correcly


% ------- EDIT BY USER ----------
run_name = { 'Terelj_test'};
meta.run_number = '20150901';

depths = -0.025;

%%
figure
hold on
for i = 1:length(run_name)
    meta.run_name = cell2mat(run_name(i));
    result = read_ground_temperatures(meta, depths);
    
    plot(result.time,result.T) 
    drawnow
    disp(['MAGST = ' num2str(mean(result.T(1,:))) '{\circ}C'])
end

datetick

title(['Ground surface temperature: ' num2str(depths) 'm'])
ylabel('Temperature [{\circ}C]')
legend(run_name)
xlim([result.time(1) result.time(end)])