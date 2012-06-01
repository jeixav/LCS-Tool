function strainlinePosition = compute_strainline2(flow,strainline,...
    position,scalingFlag)

strainlineIc = initialize_ic_grid(strainline.resolution,flow.domain);

if ~isfield(flow,'cgStrain')
    [flow.cgEigenvector,flow.cgEigenvalue] = arrayfun(@eig_array,...
    cgStrain(:,1),cgStrain(:,2),cgStrain(:,3),...
    'UniformOutput',false);
end

tspan = 0:strainline.timestep:strainline.finalTime;
strainlinePositionFw = xi1_tracing(strainlineIc,tspan,position,...
    flow.resolution,eigenvalue,eigenvector,scalingFlag);
strainlinePositionBw = xi1_tracing(strainlineIc,-tspan,position,...
    flow.resolution,eigenvalue,eigenvector,scalingFlag);

strainlinePosition = cellfun(...
    @(fw,bw)remove_repeated_ic(fw,bw),...
    strainlinePositionFw,strainlinePositionBw,...
    'UniformOutput',false);

% % Discard first line of backward time since it repeats forward time
% if any(x_fw(1,:) ~= x_bw(1,:)) || any(y_fw(1,:) ~= y_bw(1,:))
%     warning('Main:FirstLinesUnequal', ...
%         'Initial positions of forward and backward times unequal')
% else
%     position = [flipud(x_fw); x_bw(2:end,:)];
%     position(:,:,2) = [flipud(y_fw); y_bw(2:end,:)];
% end

end

function position = remove_repeated_ic(positionFw,positionBw)
%REMOVE_REPEATED_IC Remove repeated initial condition.
%   Forward and backward time strainline integrations start from the same
%   initial position. This function concatenates the forward and backward
%   position arrays, but first removes the repeated position.

if ~all(positionFw(1,:) == positionBw(1,:))
    error('Initial positions of forward and backward integrations unequal')
end

position = [flipud(positionFw); positionBw(2:end,:)];

end
