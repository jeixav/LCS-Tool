function auxiliaryPosition = auxiliary_position(basePosition,delta)
%AUXILIARY_POSITION Auxiliary grid positions
%   auxiliaryPosition = AUXILIARY_POSITION(position,delta)
%
%   BASE_POSITION gives coordinates in rows. For 2 dimensional problems
%   it must have two columns and one row for each coordinate. For 3
%   dimensional problems it must have three columns.
%
%   AUXILIARY_POSITION_ has the same format as INITIAL_POSITION, but each
%   row gives positions at auxiliary positions. In 2 dimensions, each
%   row has 8 columns with the format:
%      X+DELTA Y X-DELTA Y X Y+DELTA X Y-DELTA

if size(basePosition,2) == 2
    nColAuxiliaryPosition = 8;
elseif size(basePosition,2) == 3
    nColAuxiliaryPosition = 12;
else
    error('Auxiliary position:Column Number Error')
end

auxiliaryPosition = nan(size(basePosition,1),nColAuxiliaryPosition);

auxiliaryPosition(:,1) = basePosition(:,1) + delta;
auxiliaryPosition(:,2) = basePosition(:,2);
auxiliaryPosition(:,3) = basePosition(:,1) - delta;
auxiliaryPosition(:,4) = basePosition(:,2);
auxiliaryPosition(:,5) = basePosition(:,1);
auxiliaryPosition(:,6) = basePosition(:,2) + delta;
auxiliaryPosition(:,7) = basePosition(:,1);
auxiliaryPosition(:,8) = basePosition(:,2) - delta;

end
