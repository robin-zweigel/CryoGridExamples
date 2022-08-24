function result = read_aggregate_surfaceTEMP(run_path, variables)
files = dir(fullfile(run_path,'*.mat'));

result.time = [];
result.area = [];
result.class = {};
for var_i = 1:length(variables)
    result.(string(variables(var_i))) = [];
end

% go though all the run files
for run_i = 1:length(files)
    load([run_path '\' files(run_i).name])
    
    result.time = [result.time out.TIMESTAMP];
    
    % go through all the timesteps
    for time_i = 1:length(out.TIMESTAMP)
        I = find(result.time == out.TIMESTAMP(time_i));
        
        % go through all the stratigraphy
        for strat_i = 1:size(out.STRATIGRAPHY{1, 1})
            CURRENT = class(out.STRATIGRAPHY{1, time_i}{strat_i, 1});
            
            result.area(strat_i,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.area(1);
            result.class(strat_i,I) = {CURRENT};
            
            % go through all the variables
            for var_i = 1:length(variables)
                if isfield(out.STRATIGRAPHY{1, time_i}{strat_i, 1}.TEMP,(string(variables(var_i))))
                    result.(string(variables(var_i)))(strat_i,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.TEMP.(string(variables(var_i)));
                elseif isfield(out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR,(string(variables(var_i))))
                    result.(string(variables(var_i)))(strat_i,I) = out.STRATIGRAPHY{1, time_i}{strat_i, 1}.STATVAR.(string(variables(var_i)));
                else
                    result.(string(variables(var_i)))(strat_i,I) = NaN;
                end
            end
        end
        
    end
    
end

end