clear all
close all

folder_path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\';
source_extention = 'CryoGrid\source\';
result_extention = 'CryoGridExamplesMongolia\results\';
forcing_extention = 'CryoGridExamplesMongolia\forcing\Terelj';
addpath(genpath([folder_path source_extention])); % Required to read classes correcly


% =========================
%   start EDIT by user
%
run_name = 'Terelj_reference_run_no_snow';
run_number = '20120901';
%
%   end EDIT by user
% =========================

% Result data
load([folder_path result_extention run_name '\' run_name '_' run_number '.mat'])

% Forcing data
temp.time = squeeze(ncread([folder_path forcing_extention '\t2m.nc'], 'time'));
temp.tp = [0; 0; squeeze(ncread([folder_path forcing_extention '\tp.nc'], 'tp'))];
temp.t2m = squeeze(ncread([folder_path forcing_extention '\t2m.nc'], 't2m'))-273.15;
T_snow = 0; T_rain = 0; % Thresholds for snow/rain, need to be the same as in Excel file! Find smoother solution later
temp.rainfall = temp.tp./3600 .* (double(temp.t2m >= T_rain) + double(temp.t2m < T_rain & temp.t2m > T_snow).*(1- (temp.t2m - T_snow)./ max(1e-12,T_rain - T_snow)));
temp.snowfall = temp.tp./3600 .* (double(temp.t2m <= T_snow) + double(temp.t2m < T_rain & temp.t2m > T_snow).*(temp.t2m - T_snow)./ max(1e-12,T_rain - T_snow));

%% Extract the info
% Notation: negative - loss of water, positive - gain of water

d_snow = out.TIMESTAMP.*0; swe = d_snow; area = d_snow; isChild = d_snow;
runoff_subsurface = d_snow; runoff_surface = d_snow; sublim = d_snow;
evap = d_snow; depos = d_snow; cond  = d_snow; groundwater = d_snow;
groundice = d_snow; rainfall = d_snow; snowfall = d_snow;

timestep = out.PARA.output_timestep;
result_time = out.TIMESTAMP;
forcing_time = datenum(1900,1,1) + double(temp.time)./24;

% Water balance terms:
% 1. Precipitation
rainfall = interp1(forcing_time,temp.rainfall,result_time);
snowfall = interp1(forcing_time,temp.snowfall,result_time);

c = struct; % Contains constants
c.rho_w = 1000; % Density of water [kg/m3]
c.L_s = 2.83e+06; % Latent heat of sublimation


for i = 1:size(out.STRATIGRAPHY,2)
    top_class = class(out.STRATIGRAPHY{1, i}{1, 1}); % Uppermos class in the stratigraphy, can be CHILD
    area(i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.area(1);
    isChild(i) = area(i) ~= out.STRATIGRAPHY{1, i}{end, 1}.STATVAR.area(1);
    isSnow = strcmp( top_class(1:4), 'SNOW');
    
    
    % 2. Lateral runoff [m3]
    runoff_surface(i) = out.LATERAL{1, i}{1, 1}.STATVAR.surface_run_off;
    runoff_subsurface(i) = out.LATERAL{1, i}{2, 1}.STATVAR.subsurface_run_off;
    
    % 3. evaporation/transpiration/condensation/sublimation
    % Latent heat fluxes [W] - not per unit area!
    if isChild(i) % SNOW covers only part of ground
        Qe = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.Qe*out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.area(1); % Total latent energy flux is in GROUND class
        Qe_snow = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.Qe.*area(i); % Latent energy flux from snow
        Qe_ground = Qe - Qe_snow;
    elseif isSnow
        Qe = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.Qe;
        Qe_snow = Qe;
        Qe_ground = 0;
    else % GROUND is top_class
        Qe = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.Qe;
        Qe_snow = 0;
        Qe_ground = Qe;
    end
    
    T_ground_surface = out.STRATIGRAPHY{1, i+isSnow}{1, 1}.STATVAR.T(1);
    c.L_w = 1e3.*(2500.8 - 2.36.*T_ground_surface); % Latent heat of vaporization [J/kg K]
    
    %Sublimation and depositon (resublimation) only from snow
    sublim(i) = - Qe_snow ./(c.rho_w .* c.L_s) .* double(Qe_snow > 0);
    depos(i) = - Qe_snow ./(c.rho_w .* c.L_s) .* double(Qe_snow <= 0);
    
    % Evaporation and condensation only on ground
    evap(i) = - Qe_ground ./ (c.L_w.*c.rho_w).*double(Qe_ground > 0);
    cond(i) = - Qe_ground ./ (c.L_w.*c.rho_w).*double(Qe_ground <= 0);
    
    % 4. Storage terms
    if isSnow
        d_snow(i) = sum(out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.layerThick);
        swe(i) = sum(out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.waterIce);
    end
    
    for j = 1+isSnow:size(out.STRATIGRAPHY{1, i},1)
        groundwater(i) = groundwater(i) + sum(out.STRATIGRAPHY{1, i}{j, 1}.STATVAR.water);
        groundice(i) = groundice(i) + sum(out.STRATIGRAPHY{1, i}{j, 1}.STATVAR.ice);
    end
    
end

% Convert from accumulated volume [m^3] to flux [m^3/s]
runoff_surface = - [ 0 diff(runoff_surface)./ (24*3600*timestep) ];
runoff_subsurface = - [ 0 diff(runoff_subsurface)./ (24*3600*timestep) ];

%% Evaluate

water_in = (rainfall + snowfall + cond + depos).*timestep*24*3600; % [m3]
water_out = (evap + sublim + runoff_surface + runoff_subsurface).*timestep*24*3600;  % [m3]
d_storage = -[0 diff(groundwater) + diff(groundice) + diff(swe)]; % [m3]

% water_out = (mean(evap) + mean(sublim) + mean(runoff_surface) + mean(runoff_subsurface)).*365*24*3600
% water_in = (mean(snowfall) + mean(rainfall) + mean(cond) + mean (depos)).*365*24*3600
% d_storage = (groundwater(end) + groundice(end) + swe(end) - groundwater(1) - groundice(1) - swe(1))
