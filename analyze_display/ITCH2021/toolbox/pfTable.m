function [ index, index_alongTime, ALT, depth_alongTime, z_alongTime ] = pfTable(Tmat,z_surf,z_midCell)
% Find the permafrost table from a cryoGrid temperature matrix (Tmat, depth
% along rows and time along columns). The function can be used in "light
% mode", where just the temperature matrix is inputed, or using data about
% the model grid to produce serie of ALT and ALz
% L.C.P. Martin, Utrecht University, 2020.

[~,nbt]=size(Tmat);

% Find index
frozen = Tmat<=0;
frozen2 = sum(frozen,2);
frozen3 = frozen2==nbt;
frozen4 = find(frozen3==1);

if isempty(frozen4) % No permafrost, trivial outputs
    index=NaN;
    ALT=NaN;
    index_alongTime=nan(nbt,1);
    depth_alongTime=nan(nbt,1);
    z_alongTime=nan(nbt,1);
    
else % There is permafrost
    
    % find index of frozen ground along time
    index=frozen4(1);
    Tmat=Tmat(1:index,:);
    unfrozen=Tmat>0; % Localize unfrozen ground
    unfrozen2=cumsum(unfrozen,'reverse'); % Sum logical indexing to have ones in the first unfrozen cell
    unfrozen3=unfrozen2==1; % Find the ones
    unfrozen4=find(unfrozen3); % Locate the index of the first unfrozen cell
    [row,col] = ind2sub(size(Tmat),unfrozen4); % convert in row col indices
    index_alongTime=nan(nbt,1); % Keep row indices
    index_alongTime(col)=row+1;
    
    if nargin==1 % Function used in light mode
        ALT=NaN;
        depth_alongTime=nan(nbt,1);
        z_alongTime=nan(nbt,1);
    else % Fuction used with elevation data
        
        % Attribute ALT
        ALT=z_surf - z_midCell(index);
        
        % find depth and z of frozen ground along time
        dealWithNaN=index_alongTime;
        dealWithNaN(isnan(dealWithNaN))=1;
        z_alongTime = z_midCell(dealWithNaN);
        depth_alongTime = z_surf - z_midCell(dealWithNaN);
        z_alongTime = z_alongTime.*(index_alongTime./index_alongTime); % Bring back the NaN where they should be
        depth_alongTime = depth_alongTime.*(index_alongTime./index_alongTime); % Bring back the NaN where they should be
        
    end
    
end

end