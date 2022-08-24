function result = read_aggregate_out(run_path, variables)
files = dir(fullfile(run_path,'*.mat'));

snowLayers = 20;

result.time = [];
% result.GROUND = struct(string(variable),[]);
result.GROUND.layerThick = [];
result.GROUND.area = [];
result.GROUND.class = {};
% result.VEGETATION = struct(string(variable),[]);
result.VEGETATION.layerThick = [];
result.VEGETATION.area = [];
result.VEGETATION.class = {};
% result.SNOW = struct(string(variable),[]);
result.SNOW.layerThick = [];
result.SNOW.area = [];
result.SNOW.class = {};
for var_i = 1:length(variables)
    result.GROUND.(string(variables(var_i))) = [];
    result.VEGETATION.(string(variables(var_i))) = [];
    result.SNOW.(string(variables(var_i))) = [];
end

% Go through all run files
for run_i = 1:length(files)
    load([run_path '\' files(run_i).name])
    
    result.time = [result.time out.TIMESTAMP];
    
    % Go through all timesteps
    for time_i = 1:length(out.TIMESTAMP)
        I = find(result.time == out.TIMESTAMP(time_i));
        
        % Go through all the stratigrapy
        for strat_i = 1:size(out.STRATIGRAPHY{1, time_i})
            CURRENT = class(out.STRATIGRAPHY{1, time_i}{strat_i, 1});
            pos = 0;
            
            if strcmp(CURRENT(1:10),'VEGETATION')
                for var_i = 1:length(variables)
                    result.VEGETATION.(string(variables(var_i)))(:,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.(string(variables(var_i)));
                end
                result.VEGETATION.layerThick(:,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.layerThick;
                result.VEGETATION.area(:,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.area;
                result.VEGETATION.class(:,I) = {CURRENT};
            end
            
            if strcmp(CURRENT(1:6),'GROUND')
                pos = max(pos) + [1:length(out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.layerThick)]; % Account for several GROUND classes
                for var_i = 1:length(variables)
                    result.GROUND.(string(variables(var_i)))(pos,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.(string(variables(var_i)));
                end
                result.GROUND.layerThick(pos,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.layerThick;
                result.GROUND.area(pos,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.area;
                result.GROUND.class(pos,I) = {CURRENT};
            end
            
            if strcmp(CURRENT(1:4),'SNOW')
                if isempty(result.SNOW.layerThick) % First create space for snow, since it has dynamic size
                    for var_i = 1:length(variables)
                        result.SNOW.(string(variables(var_i)))(:,:) = NaN(snowLayers,length(result.time));
                    end
                    result.SNOW.layerThick(:,:) = NaN(snowLayers,length(result.time));
                    result.SNOW.area(:,:) = NaN(snowLayers,length(result.time));
                end
                
                J = (snowLayers - length(out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.layerThick) + 1):snowLayers;
                for var_i = 1:length(variables)
                    result.SNOW.(string(variables(var_i)))(J,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.(string(variables(var_i)));
                end
                result.SNOW.layerThick(J,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.layerThick;
                result.SNOW.area(J,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.area;
                result.SNOW.class(J,I) = {CURRENT};
            end
            
        end
        
    end
    
end



end