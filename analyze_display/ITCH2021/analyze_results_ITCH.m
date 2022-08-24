%% 0. Set up the path for CryoGrid\source
clearvars -except results;
source_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGrid\source'; % <- here locate the location of CryoGrid\source to give knowledge of the class 
addpath(genpath(source_path));
addpath('toolbox\'); addpath('toolbox\colormaps\')

%% 1. Import, process and aggregate your results
    % 1.1   Indicate locations of your results and the simulations you look at
CGsim.path='results';
CGsim.run='Terelj_flat_reference_no_snow';
    % 1.2   Interpolation (input interpVect=NaN if not used)
interpVect=[0.05 3]; % first value: cell thickness, second value depth interpolated from surface, everything in meters
    % 1.3   Extract data from run
if ~exist('results','var')
    results = aggregateResults(fullfile(CGsim.path,CGsim.run),interpVect);
    save(['results_' datestr(now,'yymmdd')],'results');
    fprintf('results saved\n')
end

%% 2. Plot your results
close all
    % 2.0   Define dates
startdate = datenum('01-09-2000','dd-mm-yyyy');
enddate = datenum('01-09-2020','dd-mm-yyyy');

    % 2.1   Tempearture plot along depth and time
fh = figure(1); fh.WindowState = 'maximized';
plotTemp(CGsim.run,startdate,enddate, results.TIMESTAMP,results.T, results.z);

    % 2.2   Ice plot along depth and time
fh = figure(2); fh.WindowState = 'maximized';
plotIce(CGsim.run,startdate,enddate, results.TIMESTAMP, results.ice, results.z)

    % 2.3   Water plot along depth and time
fh = figure(3); fh.WindowState = 'maximized';
plotWater(CGsim.run,startdate,enddate, results.TIMESTAMP, results.water, results.z)

    % 2.4   Water and ice plot along depth and time
fh = figure(4); fh.WindowState = 'maximized';
plotWaterIce(CGsim.run,startdate,enddate, results.TIMESTAMP, results.waterIce, results.z)

    % 2.5   Active layer plot
fh = figure(5); fh.WindowState = 'maximized';
plotAL(results.TIMESTAMP,results.frozenGround_depth)

    % 2.6   Working with a given variable at a given depth
variable = results.T; % Here choose which variable
depth4seasonalPlot=0.5; % depth at wich the data is read in meters
[~,ind4seasonalPlot]=min(abs((results.z(1)-results.z)-depth4seasonalPlot));
        % 2.6.1     Seasonality of a given variable at a given depth
fh = figure(6); fh.WindowState = 'maximized';
seasonalPlot(results.TIMESTAMP,variable(ind4seasonalPlot,:),'monthly','mean',1); % here possible to replace T by ice or water and monthly by daily
title('Temperature seasonality'); ylabel('Temperature at a certain depth (°C)','FontWeight','bold'); xlabel('Year','FontWeight','bold');
        % 2.6.2     Yeraly trend
fh = figure(7); fh.WindowState = 'maximized';
yearlyPlot(results.TIMESTAMP,variable(ind4seasonalPlot,:));
title('Temperature long term'); ylabel('Mean annual temperature at a certain depth (°C)','FontWeight','bold'); xlabel('Year','FontWeight','bold');
