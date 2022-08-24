clear all

folder_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\';
source_extention = 'CryoGrid\source\';
result_extention = 'CryoGridExamplesMongolia\results\';
forcing_extention = 'CryoGridExamplesMongolia\forcing\Terelj\';
addpath(genpath([folder_path source_extention])); % Required to read classes correcly


% =========================
%   start EDIT by user
%
run_name = 'Terelj_test';
run_number = '20160901';

%
%   end EDIT by user
% =========================

% Result data
load([folder_path result_extention run_name '\' run_name '_' run_number '.mat'])

%%

for i = 1:size(out.STRATIGRAPHY,2)
    Qe_canopy(i) = out.STRATIGRAPHY{1, i}{1, 1}.TEMP.Qe_canopy;
    Qt_sun(i) = out.STRATIGRAPHY{1, i}{1, 1}.TEMP.Qt_sun;
    Qt_sha(i) = out.STRATIGRAPHY{1, i}{1, 1}.TEMP.Qt_sha;
    Qe(i)   = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.Qe;
    Qh(i)   = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.Qh;
    T(i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T;
    Ts(i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.Ts;
    Tg(i) = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(1);
end

time = out.TIMESTAMP;
%%
figure
hold on
% plot(time,Qe_canopy)
% plot(time,Qt_sun)
% plot(time,Qt_sha)
plot(time,T)
plot(time,Tg)
plot(time,Ts)
% plot(time,Qh)
% plot(time,Qe)
datetick