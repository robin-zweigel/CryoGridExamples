function result_out = aggregateResults(sim_location, interpVect)
% Function to aggrate results from a simulation.
% Cas Renette, Léo Martin, Oslo University 2021.

output_list=dir(fullfile(sim_location,'*.mat'));

out2process=load(fullfile(sim_location,output_list(1).name));
if ~isnan(interpVect)
    result_out=usableOUT(out2process.out,interpVect);
else
    result_out=usableOUT(out2process.out);
    result_out.z=result_out.z.TopCell;
end

if length(output_list)>1
    result_out=rmfield(result_out,{'thermCond','snow','lateral','SEB','runoff','excessWater','dimensions'});
    
    h=waitbar(0,'Aggregating and processing CryoGrid results');
    
    frozenGround=cell2mat(permute(struct2cell(result_out.frozenGround.yearly),[3 1 2]));
    
    for i_out=2:length(output_list)
        waitbar(i_out/length(output_list))
        
        out2process=load(fullfile(sim_location,output_list(i_out).name));
        if ~isnan(interpVect)
            out_add=usableOUT(out2process.out,interpVect);
        else
            out_add=usableOUT(out2process.out);
        end
        
        % Aggregate things
        result_out.TIMESTAMP = [result_out.TIMESTAMP out_add.TIMESTAMP];
        result_out.T = [result_out.T out_add.T];
        result_out.waterIce = [result_out.waterIce out_add.waterIce];
        result_out.water = [result_out.water out_add.water];
        result_out.ice = [result_out.ice out_add.ice];
        result_out.air = [result_out.air out_add.air];
        frozenGround = [frozenGround; cell2mat(permute(struct2cell(out_add.frozenGround.yearly),[3 1 2])) ];
    end
    result_out.frozenGround_index=frozenGround(:,1);
    result_out.frozenGround_depth=frozenGround(:,2);
    result_out.frozenGround_z=frozenGround(:,3);
    result_out=rmfield(result_out,'frozenGround');
    close(h);
    
end

end

