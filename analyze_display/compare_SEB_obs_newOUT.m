% Plot and compare daily SEB using new OUT class
% R. B. Zweigel, October 2022
clear all
close all

result_path = 'C:\Users\robinbz\Documents\GitHub\CryoGridExamples\results\';
source_path = 'C:\Users\robinbz\Dokumenter\GitHub\CryoGrid\source\';
obs_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\PhD stuff\Data\Field data\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
run_name = 'Terelj_scaled_snow_denser_ground';

obs_filename = 'FA_daily_data.mat';
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
files = dir(fullfile([result_path run_name],'*0901.mat')); % select only regular out files
result.time = [];
result.Sin  = [];
result.Sout = [];
result.Lin  = [];
result.Lout = [];
result.Qe = [];
result.Qh = [];

for i = 1:length(files)
    load([result_path run_name '\' files(i).name])
    result.time = [result.time out.TIMESTAMP];
    result.Sin  = [result.Sin out.STATVAR.Sin];
    result.Sout = [result.Sout out.STATVAR.Sout];
    result.Lin  = [result.Lin out.STATVAR.Lin];
    result.Lout = [result.Lout out.STATVAR.Lout];
    result.Qe = [result.Qe out.STATVAR.Qe];
    result.Qh = [result.Qh out.STATVAR.Qh];
end
Sin = sum(reshape(sum(result.Sin(2:3,1:floor(length(result.Sin)/4)*4)),4,[]))./(24*60*60);
Lin = sum(reshape(sum(result.Lin(2:3,1:floor(length(result.Lin)/4)*4)),4,[]))./(24*60*60);
Sout = sum(reshape(sum(result.Sout(2:3,1:floor(length(result.Sout)/4)*4)),4,[]))./(24*60*60);
Lout = sum(reshape(sum(result.Lout(2:3,1:floor(length(result.Lout)/4)*4)),4,[]))./(24*60*60);
Qe_ground = sum(reshape(sum(result.Qe(2:3,1:floor(length(result.Qe)/4)*4)),4,[]))./(24*60*60);
Qh_ground = sum(reshape(sum(result.Qh(2:3,1:floor(length(result.Qh)/4)*4)),4,[]))./(24*60*60);
Qe_veg = sum(reshape(result.Qe(1,1:floor(length(result.Qe)/4)*4),4,[]))./(24*60*60);
Qh_veg = sum(reshape(result.Qh(1,1:floor(length(result.Qh)/4)*4),4,[]))./(24*60*60);
time = result.time(1:4:floor(length(result.Lin)/4)*4);

% Observations
load([obs_path obs_filename])

%%
figure
subplot(2,1,1)
plot(data.time,data.Sin)
hold on
plot(time,Sin)
datetick
legend('Measured','Simulated')
title('Incoming Shortwave below canopy')
ylabel('Radiation [W/m^2]')
xlim([datenum(2004,1,1) datenum(2009,9,1)])

subplot(2,1,2)
plot(data.time,data.Lin)
hold on
plot(time,Lin)
datetick
legend('Measured','Simulated')
title('Incoming Longwave below canopy')
ylabel('Radiation [W/m^2]')
xlim([datenum(2004,1,1) datenum(2009,9,1)])

%% Net radiation
data.Rin    = data.Sin + data.Lin;
data.Rout   = data.Sout + data.Lout;
data.Rnet   = data.Rin - data.Rout;

Rin     = Sin + Lin;
Rout    = Sout + Lout;
Rnet    = Rin - Rout;

figure
plot(data.time,data.Rin)
hold on
plot(result.time(1:4:floor(length(result.Lin)/4)*4),Rin)
datetick
ylabel('Incoming radiation total [W/m^2]')

figure
plot(data.time,data.Rout)
hold on
plot(result.time(1:4:floor(length(result.Lin)/4)*4),Rout)
datetick
ylabel('Outgoing radiation total [W/m^2]')

figure
plot(data.time,data.Rnet)
hold on
plot(time,Rnet)
ylabel('Net radiation [W/m^2]')
datetick