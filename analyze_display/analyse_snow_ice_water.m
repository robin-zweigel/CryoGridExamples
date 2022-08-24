%% extract and analyse ice, snow and water content and fluxes
% R. B. Zweigel, August 2022

path = '\\kant\geo-geohyd-u1\robinbz\Dokumenter\GitHub\CryoGridExamples\results\Terelj_NS_CLM5_long\';
files = dir(fullfile(path, '*.mat'));

data.time = [];
data.evap = [];
data.transp = [];
data.sublim = [];
data.SWE = [];
data.d_snow = [];

for k = 1:length(files)
    load([path files(k,1).name]);
    
    for i = 1:length(out.TIMESTAMP)
        for j = 1:2
            temp.evap(j,i) = out.STRATIGRAPHY{1,i}{j,1}.TEMP.evap;
            temp.sublim(j,i) = out.STRATIGRAPHY{1,i}{j,1}.TEMP.sublim;
        end
        temp.transp(i) = out.STRATIGRAPHY{1,i}{1,1}.TEMP.transp;
        
        groundclass = class(out.STRATIGRAPHY{1,i}{2,1});
        if strcmp(groundclass(1:4),'SNOW')
            temp.SWE(i) = sum(out.STRATIGRAPHY{1,i}{2,1}.STATVAR.waterIce);
            temp.d_snow(i) = sum(out.STRATIGRAPHY{1,i}{2,1}.STATVAR.layerThick);
            %     elseif strcmp(groundclass(1:6),'GROUND') && ~isempty(out.STRATIGRAPHY{1,i}{2,1}.CHILD)
            %         temp.SWE(i) = out.STRATIGRAPHY{1,i}{2,1}.CHILD.STATVAR.waterIce;
            %         temp.d_snow(i) = out.STRATIGRAPHY{1,i}{2,1}.CHILD.STATVAR.layerThick./out.STRATIGRAPHY{1,i}{2,1}.CHILD.STATVAR.area;
        else
            temp.SWE(i) = 0;
            temp.d_snow(i) = 0;
        end
        
    end
    
    data.time = [data.time out.TIMESTAMP];
    data.evap = [data.evap temp.evap];
    data.sublim = [data.sublim temp.sublim];
    data.transp = [data.transp temp.transp];
    data.SWE = [data.SWE temp.SWE];
    data.d_snow = [data.d_snow temp.d_snow];
    clear temp
end


