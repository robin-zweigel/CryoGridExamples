% Read vars
% R. B. Zweigel, January 2022
clear all
% close all

result_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGridExamples\results\';
source_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGrid\source';
addpath(genpath(source_path));
addpath(genpath('toolbox\'))

% == == == EDIT BY USER == == ==
%

run = 'Terelj_NS_CLM5_drainage';

%grid_vec = [.025 2]; % specify grid to interpolate to - [Spacing Depth]

TempVars = {'Sin','Lin','Lout','Sout','LAI','Qe','Qh','evap','sublim','transp'};
%
% == == == END EDIT == == == ==

% extract and aggregate specified variable
data = read_aggregate_surfaceTEMP([result_path run], TempVars);