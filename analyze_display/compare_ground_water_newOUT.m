% Plot and compare daily ground water using new OUT class
% R. B. Zweigel, October 2022
% clear all
% close all

result_path = 'C:\Users\robinbz\Documents\GitHub\CryoGridExamples\results\';
source_path = 'C:\Users\robinbz\Dokumenter\GitHub\CryoGrid\source\';
obs_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\PhD stuff\Data\Field data\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
run_name = 'Terelj_Xice_no_seepage';

obs_filename = 'FA_daily_soil.mat';
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
files = dir(fullfile([result_path run_name],'*0901.mat')); % select only regular out files
result.time = [];
result.height = [];
result.water = [];
result.ice = [];
for i = 1:length(files)
    load([result_path run_name '\' files(i).name])
    result.time = [result.time out.TIMESTAMP];
    result.height = out.HEIGHTS; % does not change temporally, but needs to be stored
    result.water = [result.water out.STATVAR.water];
    result.ice = [result.ice out.STATVAR.ice];
end

% Observations
load([obs_path obs_filename])
%%

z_mois = max(0.05,data.z_mois./100); % from cm depth to m
[~, unique_t, ~] = unique(result.time); % Some save date is duplicated in some out
mois_simulated = interp2(result.time(unique_t),result.height,result.water(:,unique_t),result.time(unique_t)',1500-z_mois);

mois_simulated = squeeze(mean(reshape(mois_simulated,6,4,[]),2));
time = result.time(unique_t(1:4:end));

figure
for i = 1:length(z_mois)
    subplot(length(z_mois),1,i)
    plot(time,mois_simulated(i,:))
    datetick
    hold on
%     plot(data.time,data.mois(:,i))
    xlim([data.time(1) result.time(end)])
    title(['Ground temperature at ' num2str(data.z_mois(i)) 'cm depth'])
%     if mod(i,2) == 0 
%         ylabel('Temperature [^{\circ}C]')
%     end
%     h = gca;
%     if h.YLim(2)-h.YLim(1) < minYLim
%         diffY = minYLim - (h.YLim(2)-h.YLim(1));
%         ylim([h.YLim(1)-diffY/2 h.YLim(2)+diffY/2])
%     end
end
legend('Simulated','Measured')