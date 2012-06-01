function strainlinePosition = xi1_tracing(strainlineInitialPosition,...
    tspan,cgPosition,flowResolution,cgEigenvalue,cgEigenvector,scalingFlag)

x0 = strainlineInitialPosition(:,1);
y0 = strainlineInitialPosition(:,2);
xb = x0;
yb = y0;

v1 = cgEigenvector(:,1:2);
clear eigenvector

% Scaling near fixed points
if scalingFlag
    l1 = cgEigenvalue(:,1);
    l2 = cgEigenvalue(:,2);
    v1(:,1) = ((l2 - l1)./(l2 + l1)).^2.*v1(:,1);
    v1(:,2) = ((l2 - l1)./(l2 + l1)).^2.*v1(:,2);
end
clear eigenvalue

% FIXME Cartesian grid is assumed with interp2
position_x_meshgrid = reshape(cgPosition(:,1),fliplr(flowResolution));
position_y_meshgrid = reshape(cgPosition(:,2),fliplr(flowResolution));
v1_x_meshgrid = reshape(v1(:,1),fliplr(flowResolution));
v1_y_meshgrid = reshape(v1(:,2),fliplr(flowResolution));

% Set boundary points to 0 to stop strainlines at domain boundaries
% FIXME Check if this should be done with ODESET's event location
% capabilities
v1_x_meshgrid(1,:) = 0;
v1_y_meshgrid(1,:) = 0;
v1_x_meshgrid(end,:) = 0;
v1_y_meshgrid(end,:) = 0;
v1_x_meshgrid(:,1) = 0;
v1_y_meshgrid(:,1) = 0;
v1_x_meshgrid(:,end) = 0;
v1_y_meshgrid(:,end) = 0;

global XI1_1B XI1_2B
% XI1_1B = interp2(position_x_meshgrid,position_y_meshgrid,v1_x_meshgrid,...
%     xb,yb,'nearest');
% XI1_2B = interp2(position_x_meshgrid,position_y_meshgrid,v1_y_meshgrid,...
%     xb,yb,'nearest');
% save xi1_old.mat xi1_1b xi1_2b
[XI1_1B XI1_2B] = smooth_xi1(xb,yb,position_x_meshgrid,position_y_meshgrid,cat(3,v1_x_meshgrid,v1_y_meshgrid));

v1_meshgrid(:,:,1) = v1_x_meshgrid;
v1_meshgrid(:,:,2) = v1_y_meshgrid;

% FIXME Using MaxStep option is not recommended
ode45_options = odeset('MaxStep',tspan(2) - tspan(1));
solution = ode45(@(t,y)xi1(t,y,v1_meshgrid,position_x_meshgrid,...
    position_y_meshgrid),[tspan(1) tspan(end)],[xb;yb],ode45_options);
Ftmp2 = deval(solution,tspan);
matlinX = transpose(Ftmp2(1:end/2,:));
matlinY = transpose(Ftmp2(end/2+1:end,:));

nInitialPosition = size(strainlineInitialPosition,1);
strainlinePosition = arrayfun(...
    @(col)remove_repeated_points([matlinX(:,col) matlinY(:,col)]),...
    1:nInitialPosition,...
    'UniformOutput',false);

% usePCT = false;
% 
% if usePCT
%     local_rhs = @(t,y)xi1(t,y,v1_meshgrid,position_x_meshgrid,...
%         position_y_meshgrid,interpMethod);
%     spmd
%         ic = codistributed([xb yb],codistributor1d(1));
%         ic = to_coupled(ic);
%         solution = ode45(local_rhs,[tspan(1) tspan(end)],ic,...
%             ode45_options);
%     end
% end
    
end

function position = remove_repeated_points(position)
%REMOVE_REPEATED_POINTS Removes repeated points from strainline.
%   Strainline integration can produce results where the first points of a
%   strainline, or the last points, repeat the same positions. This
%   function removes such repeated points from strainlines.

normPosition = arrayfun(@(row) norm(position(row,:)),1:size(position,1));
diffNormPosition = diff(normPosition);
beginIndex = find(diffNormPosition,1);
endIndex = find(diffNormPosition,1,'last')+1;
position = position(beginIndex:endIndex,:);

end
