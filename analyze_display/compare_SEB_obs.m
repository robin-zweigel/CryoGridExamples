% Plot and compare daily SEB
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

TempVars = {'Sin','Lin','Lout','Sout','LAI','Qe','Qh','evap','sublim','transp'};

obs_filename = 'FA_daily_data.mat';
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
temp = read_aggregate_surfaceTEMP([result_path run], TempVars);

% Observations
load([obs_path obs_filename])

%% Plot

figure
plot(data.time,data.Sin)
hold on
plot(temp.time(1:4:end),squeeze(mean(reshape(temp.Sin(2,:),[],length(temp.time)/4),1)))
