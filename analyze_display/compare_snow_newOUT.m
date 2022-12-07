% Compare snow data from measurements and simulations

% R. B. Zweigel, October 2022
clear all
% close all

result_path = 'C:\Users\robinbz\Documents\GitHub\CryoGridExamples\results\';
source_path = 'C:\Users\robinbz\Dokumenter\GitHub\CryoGrid\source\';
obs_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\PhD stuff\Data\Field data\';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%
run_name = 'Terelj_scaled_snow';

obs_filename = 'FA_daily_data.mat';
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
files = dir(fullfile([result_path run_name],'*0901.mat')); % select only regular out files
result.time = [];
result.d_snow = [];
result.swe = [];
result.snowfall = [];
result.rainfall = [];
for i = 1:length(files)
    load([result_path run_name '\' files(i).name])
    result.time = [result.time out.TIMESTAMP];
    result.d_snow = [result.d_snow out.STATVAR.d_snow]; % does not change temporally, but needs to be stored
    result.swe = [result.swe out.STATVAR.swe];
    result.snowfall = [result.snowfall out.STATVAR.snowfall];
    result.rainfall = [result.rainfall out.STATVAR.rainfall];
end


% Observations
load([obs_path obs_filename])

%%
figure
time = floor(result.time(1:4:end));
snowfall = sum(reshape(result.snowfall,4,[]));
rainfall = sum(reshape(result.rainfall,4,[]));

subplot('position',[.05 .65 .9 .30])
bar(time,snowfall+rainfall,'k')
hold on
bar(time,rainfall,'b')
datetick
xlim([min(data.time_precip) max(result.time)])
legend('Snowfall','Rainfall','Location','northwest')
ylabel('Precipitation [m/day]')

% Snow depths through time
subplot('position',[.05 .25 .9 .35])
plot(result.time,result.d_snow)
datetick
hold on
plot(data.time_precip,data.snow./100,'.')
xlim([min(data.time_precip) max(result.time)])
legend('Reference simulation','Measurements','Location','northwest')
ylabel('Snow depth [m]')

% Linear fit between 
% Simulated snowdepth at measurement times
d_snow = interp1(result.time,result.d_snow,data.time_precip);

lims = [0 .32];
no_yrs = 7;

for i = 1:no_yrs
    subplot('position',[.05+.9*(i-1)/no_yrs .05 .9/(no_yrs+1) .15])
    I = find(data.time_precip > datenum(2002+i,9,1) & data.time_precip < datenum(2003+i,9,1));
    ratio = max(data.snow(I)./100)/max(d_snow(I));
    plot(lims,lims,'k')
    hold on
    scatter(d_snow(I),data.snow(I)./100)
    best_fit = fitlm(d_snow(I),data.snow(I)./100,'Intercept',false);
    plot(lims,best_fit.Coefficients.Estimate.*lims)
    xlim(lims); ylim(lims)
    text(.1,.02,['Coefficient = ' num2str(best_fit.Coefficients.Estimate,4)])
    text(.165,.045,['ratio = ' num2str(ratio,4)])
    text(.02,.3,[num2str(2002+i) '/' num2str(2003+i)])
    if i==1
        ylabel('Measured snow depth [m]')
    elseif i == 4
        xlabel('Reference simulation snow depth [m]')
    end
end

    