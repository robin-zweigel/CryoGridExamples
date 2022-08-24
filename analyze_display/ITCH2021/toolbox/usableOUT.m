function [result] = usableOUT(out, interpVect)
% Functions that transform the native CG outputs into more usable ones.
%
% Input argument : the out file from CG
% Ouput argument : the tidied output
%
% Note 1 : This functions does not define the class that are used. If
% needed, run the first line of the main CryoGrid program (from the
% beggining until the time integration routine, without including it) so
% that matlab knows them.
%
% Note 2 : The function is for now coded for cases without Xcess ice, with
% a fixed grid (regarding spacing, mineral and organic content) along time.
% See the fprintf for cases discrimintation.
%
% Note 3 : Arbitrary choices were made regarding the variables of interest
% and this function can provide the base for other functions to extract and
% present other variables.
%
% Note 4 : the function offer the possibility to interpolate the result at
% a new resolution. For this, uses the input interpVect variable. It has
% to include 2 element, the resolution of the interpolation and the depth
% to which the interpolation is stopped, both in meters.
%          Example : interpVect = [0.02 3];
%
% Author : Léo Martin, l.c.p.martin@uu.nl, October 2020, Oslo

% Initiate output
nb_dates=length(out.TIMESTAMP);
result.TIMESTAMP=out.TIMESTAMP;
result.T = [];
result.waterIce =[];
result.water = [];
result.ice = [];
result.air=[];
result.thermCond=[];
result.snow(nb_dates).layerThick=[];
result.snow(nb_dates).area=[];
result.snow(nb_dates).distrThick=[];
result.snow(nb_dates).rho=[];
result.snow(nb_dates).T=[];
result.snow(nb_dates).Qh=[];
result.snow(nb_dates).Qe=[];
result.snow(nb_dates).Lstar=[];
result.snow(nb_dates).albedo=[];
result.snow(nb_dates).thermCond=[];
result.snow(nb_dates).water=[];
result.snow(nb_dates).ice=[];
result.snow(nb_dates).air=[];
result.snow(nb_dates).excessWater=[];
result.lateral.LAT_WATER_RESERVOIR(nb_dates).subsurface_run_off=[];
result.SEB(nb_dates).Qh=[];
result.SEB(nb_dates).Qe=[];
result.SEB(nb_dates).Qe_pot=[];
result.SEB(nb_dates).Lout=[];
result.SEB(nb_dates).Sout=[];
result.SEB(nb_dates).Lstar=[];
result.SEB(nb_dates).u_star=[];
result.runoff=nan(nb_dates,1);
result.excessWater=nan(nb_dates,1);

% Find number of ground modules (non snow) % changed RBZ 10.01.2022
nb_ground = 0;
for i=1:length(out.STRATIGRAPHY{1,1})
    className = class(out.STRATIGRAPHY{1,1}{i,1});
    nb_ground = nb_ground + strcmp('GROUND',className(1:6));
end

% Initialize upper and lower position
ground_names=cell(1,nb_ground);
for i=1:nb_ground
    ground_names{i}=['ground' num2str(i)];
end
upperPos  =cell(nb_dates,nb_ground);
lowerPos  =cell(nb_dates,nb_ground);
layerThick=cell(nb_dates,nb_ground);
area      =cell(nb_dates,nb_ground);

% Browse and store dimension and snow
for date_i=1:length(out.TIMESTAMP)
    
    STRAT_i=out.STRATIGRAPHY{1,date_i};
    
    snow_yes_no = 0; % RBZ 10.01.2022, added new way of deriving if snow, compatible with vegetation
    for snow_i = length(STRAT_i)
        className = class(STRAT_i(snow_i));
        if strcmp('SNOW',className(1:4))
            snow_yes_no = 1;
        end
    end
    % Store ground information
    for layer_i=length(STRAT_i)-nb_ground+1:1:length(STRAT_i)
        upperPos{date_i,layer_i-(length(STRAT_i)-nb_ground)}  =STRAT_i{layer_i,1}.STATVAR.upperPos;
        lowerPos{date_i,layer_i-(length(STRAT_i)-nb_ground)}  =STRAT_i{layer_i,1}.STATVAR.lowerPos;
        layerThick{date_i,layer_i-(length(STRAT_i)-nb_ground)}=STRAT_i{layer_i,1}.STATVAR.layerThick;
        area{date_i,layer_i-(length(STRAT_i)-nb_ground)}=STRAT_i{layer_i,1}.STATVAR.area;
    end
    
    %Store snow information % Changed RBZ 10.01.2022, new way of checking
    %if snow, see above
    if ~snow_yes_no
        result.snow(date_i).layerThick=NaN;
        result.snow(date_i).area=NaN;
        result.snow(date_i).distrThick=NaN;
        result.snow(date_i).rho=NaN;
        result.snow(date_i).T=NaN;
        result.snow(date_i).Qh=NaN;
        result.snow(date_i).Qe=NaN;
        result.snow(date_i).Lstar=NaN;
        result.snow(date_i).albedo=NaN;
        result.snow(date_i).thermCond=NaN;
        result.snow(date_i).water=NaN;
        result.snow(date_i).ice=NaN;
        result.snow(date_i).air=NaN;
        result.snow(date_i).excessWater=NaN;
    else
        result.snow(date_i).layerThick=sum(STRAT_i{1,1}.STATVAR.layerThick);
        result.snow(date_i).area=mean(STRAT_i{1,1}.STATVAR.area);
        result.snow(date_i).distrThick=(result.snow(date_i).layerThick * result.snow(date_i).area)/STRAT_i{length(STRAT_i)-nb_ground+1,1}.STATVAR.area(end);
        result.snow(date_i).rho=mean(STRAT_i{1,1}.STATVAR.waterIce./STRAT_i{1,1}.STATVAR.layerThick ./STRAT_i{1,1}.STATVAR.area .*1000);
        result.snow(date_i).T=mean(STRAT_i{1,1}.STATVAR.T);
        result.snow(date_i).Qh=STRAT_i{1,1}.STATVAR.Qh;
        result.snow(date_i).Qe=STRAT_i{1,1}.STATVAR.Qe;
        result.snow(date_i).Lstar=STRAT_i{1,1}.STATVAR.Lstar;
        result.snow(date_i).albedo=STRAT_i{1,1}.STATVAR.albedo;
        result.snow(date_i).thermCond=mean(STRAT_i{1,1}.STATVAR.thermCond);
        result.snow(date_i).water=mean(STRAT_i{1,1}.STATVAR.water./STRAT_i{1,1}.STATVAR.volume);
        result.snow(date_i).ice=mean(STRAT_i{1,1}.STATVAR.ice./STRAT_i{1,1}.STATVAR.volume);
        result.snow(date_i).air=mean(STRAT_i{1,1}.STATVAR.air./STRAT_i{1,1}.STATVAR.volume);
        result.snow(date_i).excessWater=STRAT_i{1,1}.STATVAR.excessWater;
    end
    
    % Find index of first ground
    if length(STRAT_i)==nb_ground
        index2store=1;
    else
        index2store=2;
    end
    % Store SEB info
    result.SEB(date_i).Qh=STRAT_i{index2store,1}.STATVAR.Qh;
    result.SEB(date_i).Qe=STRAT_i{index2store,1}.STATVAR.Qe;
    if isfield(STRAT_i{index2store,1}.STATVAR,'Qe_pot') % See error sometimes no Qe_pot with Richards ?
        result.SEB(date_i).Qe_pot=STRAT_i{index2store,1}.STATVAR.Qe_pot;
    else
        result.SEB(date_i).Qe_pot=NaN;
    end
    result.SEB(date_i).Lout=STRAT_i{index2store,1}.STATVAR.Lout;
    result.SEB(date_i).Sout=STRAT_i{index2store,1}.STATVAR.Sout;
    result.SEB(date_i).Lstar=STRAT_i{index2store,1}.STATVAR.Lstar;
    result.SEB(date_i).u_star=STRAT_i{index2store,1}.STATVAR.u_star;
    
    % Store scalar info from STATVAR
    if isfield(STRAT_i{index2store,1}.STATVAR,'runoff') % No runoff for VEGETATION classes
        result.runoff(date_i)=STRAT_i{index2store,1}.STATVAR.runoff;
    else
        result.runoff(date_i)=NaN;
    end
    result.excessWater(date_i)=STRAT_i{index2store,1}.STATVAR.excessWater;
    
    % Sotre Lateral information (to be adjusted to the number of lateral module and the data they produce)
%     result.lateral.LAT_WATER_RESERVOIR(date_i).subsurface_run_off=out.LATERAL{1,date_i}{1,1}.STATVAR.subsurface_run_off;
end

% Check for trivial case with no subsidence and no grid modification
if max(cell2mat(upperPos(:,1)))<= min(cell2mat(upperPos(:,1)))
    
    % fprintf('usableOUT : No subsidence, simple depth processing\n')
    
    % Create z axis
    result.z.thick=vertcat(layerThick{1,:});
    result.z.TopCell=upperPos{1,1}-[0; cumsum(vertcat(layerThick{1,:}))];
    result.z.TopCell(end)=[];
    result.z.MidCell=result.z.TopCell-0.5.*vertcat(layerThick{1,:});
    variable_i=1:nb_ground;
    snow_daily=layerThick(1,:);
    for i=1:nb_ground
        snow_daily{i}=variable_i(i).*snow_daily{i}./snow_daily{i};
    end
    result.z.grounds=vertcat(snow_daily{:});
    
    % Fill matrixes
    result.T = nan(length(result.z.grounds),length(result.TIMESTAMP));
    result.waterIce =result.T;
    result.water = result.T;
    result.ice = result.T;
    result.air = result.T;
    result.thermCond=result.T;
    
    for date_i=1:length(out.TIMESTAMP)
        
        STRAT_i=out.STRATIGRAPHY{1,date_i};
        
        % Store ground information
        for layer_i=length(STRAT_i)-nb_ground+1:1:length(STRAT_i)
            % fprintf('\t%1.0f\n',layer_i)
            result.T(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.T;
            result.waterIce(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.waterIce;
            result.water(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.water;
            result.ice(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.ice;
            result.air(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.air;
            result.thermCond(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.thermCond;
        end
        
    end
    
    % Fill soil data
    STRAT_i=out.STRATIGRAPHY{1,1};
    result.soil.organic=nan(length(result.z.grounds),1);
    result.soil.mineral=result.soil.organic;
    % Store ground information
    for layer_i=length(STRAT_i)-nb_ground+1:1:length(STRAT_i)
        result.soil.organic(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),1)=STRAT_i{layer_i,1}.STATVAR.organic;
        result.soil.mineral(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),1)=STRAT_i{layer_i,1}.STATVAR.mineral;
    end
    
    % Finalize dimensions and snow
    ground_names=cell(1,nb_ground);
    for i=1:nb_ground
        ground_names{i}=['ground' num2str(i)];
    end
    result.dimensions.upperPos  =cell2struct(upperPos,ground_names,2);
    result.dimensions.lowerPos  =cell2struct(lowerPos,ground_names,2);
    result.dimensions.layerThick=cell2struct(layerThick,ground_names,2);
    % result.dimensions.area=cell2struct(area,ground_names,2);
    result.dimensions.area=area{1,1}(1);
    
    % Streamline dimensions
    result.dimensions.upperPos(2:end)=[];
    result.dimensions.lowerPos(2:end)=[];
    result.dimensions.layerThick(2:end)=[];
    % result.dimensions.area(2:end)=[];
    
    % Compute volumetric fractions
    [result.ice]=[result.ice]                  ./([result.z.thick]*[result.dimensions.area]);
    [result.water]=[result.water]              ./([result.z.thick]*[result.dimensions.area]);
    [result.waterIce]=[result.waterIce]        ./([result.z.thick]*[result.dimensions.area]);
    [result.air]=[result.air]                  ./([result.z.thick]*[result.dimensions.area]);
    [result.soil.organic]=[result.soil.organic]./([result.z.thick]*[result.dimensions.area]);
    [result.soil.mineral]=[result.soil.mineral]./([result.z.thick]*[result.dimensions.area]);
    
    % Interpolation
    if nargin > 1
        % Check interpolation window
        domainDepth=result.z.TopCell(1)-result.z.TopCell(end) + result.z.thick(end);
        if interpVect(2)> domainDepth
            fprintf('UsableOUT : Interpolation window deeper than modelled domain,\n            -> narrowed to model domain.\n')
            interpVect(2)=domainDepth;
        end
        % Define depth vector
        zInterp=(result.z.TopCell(1):(-1)*interpVect(1):result.z.TopCell(1)-interpVect(2))';
        % interpolate everybody
        goodZin=[result.z.TopCell(1);result.z.MidCell];
        result.T=interp1(goodZin,[ result.T(1,:) ;result.T],zInterp);
        result.water=interp1(goodZin,[ result.water(1,:) ;result.water],zInterp);
        result.air=interp1(goodZin,[ result.air(1,:) ;result.air],zInterp);
        result.ice=interp1(goodZin,[ result.ice(1,:) ;result.ice],zInterp);
        result.waterIce=interp1(goodZin,[ result.waterIce(1,:) ;result.waterIce],zInterp);
        result.soil.organic=interp1(goodZin,[ result.soil.organic(1,:) ;result.soil.organic],zInterp);
        result.soil.mineral=interp1(goodZin,[ result.soil.mineral(1,:) ;result.soil.mineral],zInterp);
        result.thermCond=interp1(goodZin,[ result.thermCond(1,:) ;result.thermCond],zInterp);
        
        result.z=zInterp;
        
    end
    
    % Add permafrost elevation info
    if nargin == 1
        [ result.frozenGround.ALT_ind, index_alongTime, result.frozenGround.ALT_val, depth_alongTime, z_alongTime ] = pfTable(result.T,result.z.TopCell(1),result.z.MidCell);
    else
        [ result.frozenGround.ALT_ind, index_alongTime, result.frozenGround.ALT_val, depth_alongTime, z_alongTime ] = pfTable(result.T,result.z(1),result.z);
    end
    result.frozenGround.yearly(nb_dates).index=[];
    add=num2cell(index_alongTime);
    [result.frozenGround.yearly.index]=add{:};
    add=num2cell(depth_alongTime);
    [result.frozenGround.yearly.depth]=add{:};
    add=num2cell(z_alongTime);
    [result.frozenGround.yearly.z]=add{:};
    
    if isfield(out.PARA,'time_average')
        % Compute daily average
        
        % Average ground data
        datetime_out=datetime(result.TIMESTAMP,'ConvertFrom','datenum');
        
        [result.TIMESTAMP,result.T]=timeAverages(datetime_out,result.T        ,out.PARA.time_average,'mean');
        [~,result.water]           =timeAverages(datetime_out,result.water    ,out.PARA.time_average,'mean');
        [~,result.air]             =timeAverages(datetime_out,result.air      ,out.PARA.time_average,'mean');
        [~,result.ice]             =timeAverages(datetime_out,result.ice      ,out.PARA.time_average,'mean');
        [~,result.waterIce]        =timeAverages(datetime_out,result.waterIce ,out.PARA.time_average,'mean');
        [~,result.thermCond]       =timeAverages(datetime_out,result.thermCond,out.PARA.time_average,'mean');
        [~,result.runoff]          =timeAverages(datetime_out,result.runoff   ,out.PARA.time_average,'lastvalue');
        [~,result.excessWater]          =timeAverages(datetime_out,result.excessWater   ,out.PARA.time_average,'lastvalue');
        
        % Average snow data
        snow_daily=struct2cell(result.snow);
        snow_daily=permute(snow_daily,[3 1 2]);
        where_empty=cellfun(@isempty,snow_daily); % Empty values replaced by NaN
        snow_daily(where_empty)={NaN};
        snow_daily=cell2mat(snow_daily);
        [~,averaged_snow] = timeAverages(datetime_out,snow_daily,out.PARA.time_average,'mean');
        [~,averaged_excessWater] = timeAverages(datetime_out,snow_daily,out.PARA.time_average,'lastvalue');
        averaged_snow(:,end)=[];
        averaged_snow=[averaged_snow averaged_excessWater(:,end)]; % Whatch out, excessWater is the only cummulative variable
        snow_daily=num2cell(averaged_snow);
        result.snow=cell2struct(snow_daily,{'layerThick','area','distrThick','rho','T','Qh','Qe','Lstar','albedo','thermCond','water','ice','air','excessWater'},2);
        
        % Average SEB data
        SEB_daily=struct2table(result.SEB);
        SEB_daily=table2array(SEB_daily)';
        [~,averaged_SEB] = timeAverages(datetime_out,SEB_daily,out.PARA.time_average,'mean');
        SEB_daily=num2cell(averaged_SEB');
        result.SEB=cell2struct(SEB_daily,{'Qh','Qe','Qe_pot','Lout','Sout','Lstar','u_star'},2);
        
        % Average Lateral reservoir data
        reservoir_daily=struct2table(result.lateral.LAT_WATER_RESERVOIR);
        reservoir_daily=table2array(reservoir_daily)';
        [~,averaged_reservoir] = timeAverages(datetime_out,reservoir_daily,out.PARA.time_average,'lastvalue');
        reservoir_daily=num2cell(averaged_reservoir');
        result.lateral.LAT_WATER_RESERVOIR=cell2struct(reservoir_daily,{'subsurface_run_off'},2);
        
        % Average frozen ground data
        frozen_daily=struct2table(result.frozenGround.yearly);
        frozen_daily=table2array(frozen_daily)';
        [~,averaged_frozen] = timeAverages(datetime_out,frozen_daily,out.PARA.time_average,'mean');
        frozen_daily=num2cell(averaged_frozen');
        result.frozenGround.yearly=cell2struct(frozen_daily,{'index','depth','z'},2);
    end
    
else
    
    fprintf('usableOUT : this function does not handle simulation results with subsidence\n            Output only partially filled.\n')
    
end

end