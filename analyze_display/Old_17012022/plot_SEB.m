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
run_number = '20140504';
class = 1;

%
%   end EDIT by user
% =========================

% Result data
load([folder_path result_extention run_name '\' run_name '_' run_number '.mat'])

%%

for i = 1:size(out.STRATIGRAPHY,2)
    Lout(i) = out.STRATIGRAPHY{1, i}{class, 1}.STATVAR.Lout;
    Sout(i) = out.STRATIGRAPHY{1, i}{class, 1}.STATVAR.Sout;
	Lin(i)  = out.STRATIGRAPHY{1, i}{class, 1}.STATVAR.Lin;
    Sin(i)  = out.STRATIGRAPHY{1, i}{class, 1}.STATVAR.Sin;
    Labs(i) = out.STRATIGRAPHY{1, i}{class, 1}.TEMP.L_abs;
    Sabs(i) = out.STRATIGRAPHY{1, i}{class, 1}.TEMP.S_abs;
    Qe(i)   = out.STRATIGRAPHY{1, i}{class, 1}.STATVAR.Qe;
    Qh(i)   = out.STRATIGRAPHY{1, i}{class, 1}.STATVAR.Qh;
end

% Lout = mean(reshape(Lout,4,[]));
% Sout = mean(reshape(Sout,4,[]));
% Lin  = mean(reshape(Lin,4,[]));
% Sin  = mean(reshape(Sin,4,[]));
% Labs = mean(reshape(Labs,4,[]));
% Sabs = mean(reshape(Sabs,4,[]));
% Qe   = mean(reshape(Qe,4,[]));
% Qh   = mean(reshape(Qh,4,[]));
% time = out.TIMESTAMP(2:4:end);
%%
figure
hold on
plot(out.TIMESTAMP,Labs)

plot(out.TIMESTAMP,Sabs)
plot(out.TIMESTAMP,Qe)
plot(out.TIMESTAMP,Qh)
datetick