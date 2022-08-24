function result = read_canopy_temperatures(meta)
load([meta.folder_path meta.result_extention meta.run_name '\' meta.run_name '_' meta.run_number '.mat'])

class_name = class(out.STRATIGRAPHY{1, 1}{1, 1});
result.time = out.TIMESTAMP;
result.class = class_name;

if strcmp(class_name(1:10),'VEGETATION')
    result.Tv = nan(size(out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.T,2),size(out.TIMESTAMP,2));
    for i = 1:size(out.STRATIGRAPHY,2)
        result.Tv(:,i) = [out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T];
        
    end
else
    
end

end