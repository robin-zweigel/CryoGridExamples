function result = read_ground_temperatures(meta, depths)
    load([meta.folder_path meta.result_extention meta.run_name '\' meta.run_name '_' meta.run_number '.mat'])
    
    result.time = out.TIMESTAMP;
    
    for i = 1:size(out.STRATIGRAPHY{1, 1},1)
        class_name = class(out.STRATIGRAPHY{1, 1}{i, 1});
        if strcmp(class_name(1:6),'GROUND')
            z_surface = out.STRATIGRAPHY{1, 1}{i, 1}.STATVAR.upperPos;
        end
    end
    result.depths = depths + z_surface;

    for i = 1:size(out.STRATIGRAPHY,2)
        temp.grid = [];
        temp.T = [];
        for j = 1:size(out.STRATIGRAPHY{1,i},1)
            temp.grid = [temp.grid; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos-cumsum(out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick) + out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick./2];
            temp.T = [temp.T; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.T(1); out.STRATIGRAPHY{1,i}{j,1}.STATVAR.T];
        end

        result.T(:,i) = interp1(temp.grid,temp.T,result.depths);

    end
end
