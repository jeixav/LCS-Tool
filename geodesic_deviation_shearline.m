function shearline = geodesic_deviation_shearline(flow,shearline,verbose)

% Create interpolants
l1Interpolant = make_interpolant(flow.domain,flow.resolution,...
    flow.cgEigenvalue(:,1));

l2Interpolant = make_interpolant(flow.domain,flow.resolution,...
    flow.cgEigenvalue(:,2));

deltaX = diff(flow.domain(1,:))/(double(flow.resolution(1)) - 1);
deltaY = diff(flow.domain(2,:))/(double(flow.resolution(2)) - 1);

[dl1(:,:,1),dl1(:,:,2)] = gradient(l1Interpolant.Values,deltaX,deltaY);
dl1Interpolant.x = make_interpolant(flow.domain,flow.resolution,dl1(:,:,1));
dl1Interpolant.y = make_interpolant(flow.domain,flow.resolution,dl1(:,:,2));
clear('dl1');

xi1Interpolant.x = make_interpolant(flow.domain,flow.resolution,...
    flow.cgEigenvector(:,1));
xi1Interpolant.y = make_interpolant(flow.domain,flow.resolution,...
    flow.cgEigenvector(:,2));

[dxi1(:,:,1),dxi1(:,:,2)] = gradient(xi1Interpolant.x.Values,deltaX,deltaY);
[dxi1(:,:,3),dxi1(:,:,4)] = gradient(xi1Interpolant.y.Values,deltaX,deltaY);
dxi1Interpolant.xx = make_interpolant(flow.domain,flow.resolution,dxi1(:,:,1));
dxi1Interpolant.xy = make_interpolant(flow.domain,flow.resolution,dxi1(:,:,2));
dxi1Interpolant.yx = make_interpolant(flow.domain,flow.resolution,dxi1(:,:,3));
dxi1Interpolant.yy = make_interpolant(flow.domain,flow.resolution,dxi1(:,:,4));
clear('dxi1')

xi2Interpolant.x = make_interpolant(flow.domain,flow.resolution,...
    flow.cgEigenvector(:,3));
xi2Interpolant.y = make_interpolant(flow.domain,flow.resolution,...
    flow.cgEigenvector(:,4));

[dxi2(:,:,1),dxi2(:,:,2)] = gradient(xi2Interpolant.x.Values,deltaX,deltaY);
[dxi2(:,:,3),dxi2(:,:,4)] = gradient(xi2Interpolant.y.Values,deltaX,deltaY);
dxi2Interpolant.xx = make_interpolant(flow.domain,flow.resolution,dxi2(:,:,1));
dxi2Interpolant.xy = make_interpolant(flow.domain,flow.resolution,dxi2(:,:,2));
dxi2Interpolant.yx = make_interpolant(flow.domain,flow.resolution,dxi2(:,:,3));
dxi2Interpolant.yy = make_interpolant(flow.domain,flow.resolution,dxi2(:,:,4));
clear('dxi2')

etaPosInterpolant.x = make_interpolant(flow.domain,flow.resolution,...
    shearline.etaPos(:,1));
etaPosInterpolant.y = make_interpolant(flow.domain,flow.resolution,...
    shearline.etaPos(:,2));

etaNegInterpolant.x = make_interpolant(flow.domain,flow.resolution,...
    shearline.etaNeg(:,1));
etaNegInterpolant.y = make_interpolant(flow.domain,flow.resolution,...
    shearline.etaNeg(:,2));

alphaInterpolant = make_interpolant(flow.domain,flow.resolution,...
    sqrt(sqrt(l2Interpolant.Values)./(sqrt(l1Interpolant.Values) ...
    + sqrt(l2Interpolant.Values))));

betaInterpolant = make_interpolant(flow.domain,flow.resolution,...
    sqrt(sqrt(l1Interpolant.Values)./(sqrt(l1Interpolant.Values) ...
    + sqrt(l2Interpolant.Values))));

[dAlpha(:,:,1),dAlpha(:,:,2)] = gradient(alphaInterpolant.Values,deltaX,deltaY);
dAlphaInterpolant.x = make_interpolant(flow.domain,flow.resolution,dAlpha(:,:,1));
dAlphaInterpolant.y = make_interpolant(flow.domain,flow.resolution,dAlpha(:,:,2));
clear('dAlpha')

devMode = false;

posOrNeg = 'pos';

if verbose
    progressMax = numel(shearline.positionPos) ...
        + numel(shearline.positionNeg);
    progressBar = ParforProgressStarter2(mfilename,progressMax);
else
    progressBar = false;
end

geodesicDeviationPos = cell(size(shearline.positionPos));
positionPos = shearline.positionPos;
parfor idx = 1:numel(geodesicDeviationPos)
    geodesicDeviationPos{idx} = ...
        geodesic_deviation_individual(positionPos{idx},...
        l1Interpolant,l2Interpolant,dl1Interpolant,alphaInterpolant,...
        betaInterpolant,dAlphaInterpolant,etaPosInterpolant,...
        xi1Interpolant,xi2Interpolant,dxi1Interpolant,dxi2Interpolant,...
        flow,posOrNeg,devMode);
    if verbose
        progressBar.increment(idx) %#ok<PFBNS>
    end
end
shearline.geodesicDeviationPos = geodesicDeviationPos;

posOrNeg = 'neg';
geodesicDeviationNeg = cell(size(shearline.positionNeg));
positionNeg = shearline.positionNeg;
parfor idx = 1:numel(geodesicDeviationNeg)
    geodesicDeviationNeg{idx} = ...
        geodesic_deviation_individual(positionNeg{idx},...
        l1Interpolant,l2Interpolant,dl1Interpolant,alphaInterpolant,...
        betaInterpolant,dAlphaInterpolant,etaPosInterpolant,...
        xi1Interpolant,xi2Interpolant,dxi1Interpolant,dxi2Interpolant,...
        flow,posOrNeg,devMode);
    if verbose
        progressBar.increment(idx) %#ok<PFBNS>
    end
end
shearline.geodesicDeviationNeg = geodesicDeviationNeg;

if verbose
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
end

shearline.averageGeodesicDeviationPos = cellfun(...
    @average_geodesic_deviation_individual,...
    shearline.positionPos,shearline.geodesicDeviationPos);

shearline.averageGeodesicDeviationPos = cellfun(...
    @average_geodesic_deviation_individual,...
    shearline.positionPos,shearline.geodesicDeviationPos);

shearline.averageGeodesicDeviationNeg = cellfun(...
    @average_geodesic_deviation_individual,...
    shearline.positionNeg,shearline.geodesicDeviationNeg);

function interpolant = make_interpolant(domain,resolution,data)

position.x = linspace(domain(1,1),domain(1,2),resolution(1));
position.y = linspace(domain(2,1),domain(2,2),resolution(2));

dataGrid = reshape(data,fliplr(resolution));

interpolant = griddedInterpolant({position.y,position.x},dataGrid);

function [discontIdx,idx1] = ...
    is_element_with_orient_discont(position,domain,resolution,vector)
% Determine if position is between grid points with an orientation
% discontinuity. Includes grid points needed for finite differences.
% Positions:
%     7   6
%
% 8   2   1   5
%       p
% 9   3   4   12
%
%    10  11

vectorX = reshape(vector(:,1),fliplr(resolution));
vectorY = reshape(vector(:,2),fliplr(resolution));

idx1 = position2index(position,domain,resolution);

% Check if boundary element
if idx1(1) == 2 || idx1(1) == resolution(1) || idx1(2) == 2 ...
        || idx1(2) == resolution(2)
    discontIdx = nan;
    return
end

idxX = idx1(1);
idxY = idx1(2);

vector1 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX - 1;
vector2 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxY = idxY - 1;
vector3 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX + 1;
vector4 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX + 1;
idxY = idxY + 1;
vector5 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX - 1;
idxY = idxY + 1;
vector6 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX - 1;
vector7 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX - 1;
idxY = idxY - 1;
vector8 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxY = idxY - 1;
vector9 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX + 1;
idxY = idxY - 1;
vector10 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX + 1;
vector11 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

idxX = idxX + 1;
idxY = idxY + 1;
vector12 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

discontIdx = find([vector2;vector3;vector4;vector5;vector6;vector7;...
    vector8;vector9;vector10;vector11;vector12]*transpose(vector1) < 0);

function [geodesicDeviation,discontPt] = geodesic_deviation_individual(...
    shearlinePosition,l1Interpolant,l2Interpolant,dl1Interpolant,...
    alphaInterpolant,betaInterpolant,dAlphaInterpolant,etaInterpolant,...
    xi1Interpolant,xi2Interpolant,dxi1Interpolant,dxi2Interpolant,flow,...
    posOrNeg,devMode)

position = shearlinePosition;

l1 = l1Interpolant(position(:,2),position(:,1));
l2 = l2Interpolant(position(:,2),position(:,1));

dl1(:,1) = dl1Interpolant.x(position(:,2),position(:,1));
dl1(:,2) = dl1Interpolant.y(position(:,2),position(:,1));

alpha = alphaInterpolant(position(:,2),position(:,1));
beta = betaInterpolant(position(:,2),position(:,1));

dAlpha(:,1) = dAlphaInterpolant.x(position(:,2),position(:,1));
dAlpha(:,2) = dAlphaInterpolant.y(position(:,2),position(:,1));

nPoints = size(position,1);

eta = nan(nPoints,2);
xi1 = nan(nPoints,2);
xi2 = nan(nPoints,2);

if devMode
    discontPt = [];
end

for iPoint = 1:nPoints

    iBE = is_boundary_element(position(iPoint,:),flow.domain,flow.resolution);
    if iBE

        % FIXME Need to check for orientation discontinuity
        eta(iPoint,1) = etaInterpolant.x(position(iPoint,2),...
            position(iPoint,1));
        eta(iPoint,2) = etaInterpolant.y(position(iPoint,2),...
            position(iPoint,1));
        
        xi1(iPoint,1) = xi1Interpolant.x(position(iPoint,2),position(iPoint,1));
        xi1(iPoint,2) = xi1Interpolant.y(position(iPoint,2),position(iPoint,1));
        
        xi2(iPoint,1) = xi2Interpolant.x(position(iPoint,2),position(iPoint,1));
        xi2(iPoint,2) = xi2Interpolant.y(position(iPoint,2),position(iPoint,1));
        
        dxi1.xx(iPoint) = dxi1Interpolant.xx(position(iPoint,2),...
            position(iPoint,1));
        dxi1.xy(iPoint) = dxi1Interpolant.xy(position(iPoint,2),...
            position(iPoint,1));
        dxi1.yx(iPoint) = dxi1Interpolant.yx(position(iPoint,2),...
            position(iPoint,1));
        dxi1.yy(iPoint) = dxi1Interpolant.yy(position(iPoint,2),...
            position(iPoint,1));
        
        dxi2.xx(iPoint) = dxi2Interpolant.xx(position(iPoint,2),...
            position(iPoint,1));
        dxi2.xy(iPoint) = dxi2Interpolant.xy(position(iPoint,2),...
            position(iPoint,1));
        dxi2.yx(iPoint) = dxi2Interpolant.yx(position(iPoint,2),...
            position(iPoint,1));
        dxi2.yy(iPoint) = dxi2Interpolant.yy(position(iPoint,2),...
            position(iPoint,1));
    else
        [discontIdx,idx1] = is_element_with_orient_discont(...
            position(iPoint,:),flow.domain,flow.resolution,...
            flow.cgEigenvector(:,1:2));
    
        if discontIdx

            if devMode
                discontPt = [discontPt iPoint]; %#ok<AGROW>
            end
            
            tmp = ones(11,1);
            tmp(discontIdx) = -1;
            
            % Continuous interpolant for eta
            deltaX = diff(flow.domain(1,:))/(double(flow.resolution(1)) - 1);
            xMin = flow.domain(1,1);
            
            deltaY = diff(flow.domain(2,:))/(double(flow.resolution(2)) - 1);
            yMin = flow.domain(2,1);
            
            position1 = [(idx1(1)-1)*deltaX+xMin (idx1(2)-1)*deltaY+yMin];
            positionX = [position1(1)-deltaX position1(1)];
            positionY = [position1(2)-deltaY position1(2)];
            
            continuousEtaPosInterpolant.x = griddedInterpolant({positionY,...
                positionX},...
                [etaInterpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2) ...
                etaInterpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3); ...
                etaInterpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                etaInterpolant.x.Values(idx1(2),idx1(1))]);
            
            eta(iPoint,1) = continuousEtaPosInterpolant.x(...
                position(iPoint,2),position(iPoint,1));
            
            continuousEtaPosInterpolant.y = griddedInterpolant({positionY,...
                positionX},...
                [etaInterpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2) ...
                etaInterpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3); ...
                etaInterpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                etaInterpolant.y.Values(idx1(2),idx1(1))]);
            
            eta(iPoint,2) = continuousEtaPosInterpolant.y(...
                position(iPoint,2),position(iPoint,1));

            % Continuous interpolant for xi1
            ctsXi1Interpolant.x = griddedInterpolant({positionY,positionX},...
                [xi1Interpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2) ...
                xi1Interpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3);
                xi1Interpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                xi1Interpolant.x.Values(idx1(2),idx1(1))]);
            ctsXi1Interpolant.y = griddedInterpolant({positionY,positionX},...
                [xi1Interpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2) ...
                xi1Interpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3);
                xi1Interpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                xi1Interpolant.y.Values(idx1(2),idx1(1))]);
            xi1(iPoint,1) = ctsXi1Interpolant.x(position(iPoint,2),position(iPoint,1));
            xi1(iPoint,2) = ctsXi1Interpolant.y(position(iPoint,2),position(iPoint,1));
            
            % Continuous interpolant for xi2
            ctsXi2Interpolant.x = griddedInterpolant({positionY,positionX},...
                [xi2Interpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2) ...
                xi2Interpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3);
                xi2Interpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                xi2Interpolant.x.Values(idx1(2),idx1(1))]);
            ctsXi2Interpolant.y = griddedInterpolant({positionY,positionX},...
                [xi2Interpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2) ...
                xi2Interpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3);
                xi2Interpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                xi2Interpolant.y.Values(idx1(2),idx1(1))]);
            xi2(iPoint,1) = ctsXi2Interpolant.x(position(iPoint,2),position(iPoint,1));
            xi2(iPoint,2) = ctsXi2Interpolant.y(position(iPoint,2),position(iPoint,1));
            
            % Continuous interpolant for dxi1
            % Position 1
            dxi1Cts.xx(1) = .5*(xi1Interpolant.x.Values(idx1(2),idx1(1)+1)*tmp(4) ...
                - xi1Interpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1))/deltaX;
            dxi1Cts.xy(1) = .5*(xi1Interpolant.x.Values(idx1(2)+1,idx1(1))*tmp(5) ...
                - xi1Interpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3))/deltaY;
            dxi1Cts.yx(1) = .5*(xi1Interpolant.y.Values(idx1(2),idx1(1)+1)*tmp(4) ...
                - xi1Interpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1))/deltaX;
            dxi1Cts.yy(1) = .5*(xi1Interpolant.y.Values(idx1(2)+1,idx1(1))*tmp(5) ...
                - xi1Interpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3))/deltaY;
            % Position 2
            dxi1Cts.xx(2) = .5*(xi1Interpolant.x.Values(idx1(2),idx1(1)) ...
                - xi1Interpolant.x.Values(idx1(2),idx1(1)-2)*tmp(7))/deltaX;
            dxi1Cts.xy(2) = .5*(xi1Interpolant.x.Values(idx1(2)+1,idx1(1)-1)*tmp(6) ...
                - xi1Interpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaY;
            dxi1Cts.yx(2) = .5*(xi1Interpolant.y.Values(idx1(2),idx1(1)) ...
                - xi1Interpolant.y.Values(idx1(2),idx1(1)-2)*tmp(7))/deltaX;
            dxi1Cts.yy(2) = .5*(xi1Interpolant.y.Values(idx1(2)+1,idx1(1)-1)*tmp(6) ...
                - xi1Interpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaY;
            % Position 3
            dxi1Cts.xx(3) = .5*(xi1Interpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3) ...
                - xi1Interpolant.x.Values(idx1(2)-1,idx1(1)-2)*tmp(8))/deltaX;
            dxi1Cts.xy(3) = .5*(xi1Interpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                - xi1Interpolant.x.Values(idx1(2)-2,idx1(1)-1)*tmp(9))/deltaY;
            dxi1Cts.yx(3) = .5*(xi1Interpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3) ...
                - xi1Interpolant.y.Values(idx1(2)-1,idx1(1)-2)*tmp(8))/deltaX;
            dxi1Cts.yy(3) = .5*(xi1Interpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                - xi1Interpolant.y.Values(idx1(2)-2,idx1(1)-1)*tmp(9))/deltaY;
            % Position 4
            dxi1Cts.xx(4) = .5*(xi1Interpolant.x.Values(idx1(2)-1,idx1(1)+1)*tmp(11) ...
                - xi1Interpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaX;
            dxi1Cts.xy(4) = .5*(xi1Interpolant.x.Values(idx1(2),idx1(1)) ...
                - xi1Interpolant.x.Values(idx1(2)-2,idx1(1))*tmp(10))/deltaY;
            dxi1Cts.yx(4) = .5*(xi1Interpolant.y.Values(idx1(2)-1,idx1(1)+1)*tmp(11) ...
                - xi1Interpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaX;
            dxi1Cts.yy(4) = .5*(xi1Interpolant.y.Values(idx1(2),idx1(1)) ...
                - xi1Interpolant.y.Values(idx1(2)-2,idx1(1))*tmp(10))/deltaY;
            
            ctsDxi1Interpolant.xx = griddedInterpolant({positionY,positionX},...
                [dxi1Cts.xx(3) dxi1Cts.xx(4); dxi1Cts.xx(2) dxi1Cts.xx(1)]);
            dxi1.xx(iPoint) = ctsDxi1Interpolant.xx(position(iPoint,2),position(iPoint,1));
            ctsDxi1Interpolant.xy = griddedInterpolant({positionY,positionX},...
                [dxi1Cts.xy(3) dxi1Cts.xy(4); dxi1Cts.xy(2) dxi1Cts.xy(1)]);
            dxi1.xy(iPoint) = ctsDxi1Interpolant.xy(position(iPoint,2),position(iPoint,1));
            ctsDxi1Interpolant.yx = griddedInterpolant({positionY,positionX},...
                [dxi1Cts.yx(3) dxi1Cts.yx(4); dxi1Cts.yx(2) dxi1Cts.yx(1)]);
            dxi1.yx(iPoint) = ctsDxi1Interpolant.yx(position(iPoint,2),position(iPoint,1));
            ctsDxi1Interpolant.yy = griddedInterpolant({positionY,positionX},...
                [dxi1Cts.yy(3) dxi1Cts.yy(4); dxi1Cts.yy(2) dxi1Cts.yy(1)]);
            dxi1.yy(iPoint) = ctsDxi1Interpolant.yy(position(iPoint,2),position(iPoint,1));
            
            % Continuous interpolant for dxi2
            % Position 1
            dxi2Cts.xx(1) = .5*(xi2Interpolant.x.Values(idx1(2),idx1(1)+1)*tmp(4) ...
                - xi2Interpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1))/deltaX;
            dxi2Cts.xy(1) = .5*(xi2Interpolant.x.Values(idx1(2)+1,idx1(1))*tmp(5) ...
                - xi2Interpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3))/deltaY;
            dxi2Cts.yx(1) = .5*(xi2Interpolant.y.Values(idx1(2),idx1(1)+1)*tmp(4) ...
                - xi2Interpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1))/deltaX;
            dxi2Cts.yy(1) = .5*(xi2Interpolant.y.Values(idx1(2)+1,idx1(1))*tmp(5) ...
                - xi2Interpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3))/deltaY;
            % Position 2
            dxi2Cts.xx(2) = .5*(xi2Interpolant.x.Values(idx1(2),idx1(1)) ...
                - xi2Interpolant.x.Values(idx1(2),idx1(1)-2)*tmp(7))/deltaX;
            dxi2Cts.xy(2) = .5*(xi2Interpolant.x.Values(idx1(2)+1,idx1(1)-1)*tmp(6) ...
                - xi2Interpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaY;
            dxi2Cts.yx(2) = .5*(xi2Interpolant.y.Values(idx1(2),idx1(1)) ...
                - xi2Interpolant.y.Values(idx1(2),idx1(1)-2)*tmp(7))/deltaX;
            dxi2Cts.yy(2) = .5*(xi2Interpolant.y.Values(idx1(2)+1,idx1(1)-1)*tmp(6) ...
                - xi2Interpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaY;
            % Position 3
            dxi2Cts.xx(3) = .5*(xi2Interpolant.x.Values(idx1(2)-1,idx1(1))*tmp(3) ...
                - xi2Interpolant.x.Values(idx1(2)-1,idx1(1)-2)*tmp(8))/deltaX;
            dxi2Cts.xy(3) = .5*(xi2Interpolant.x.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                - xi2Interpolant.x.Values(idx1(2)-2,idx1(1)-1)*tmp(9))/deltaY;
            dxi2Cts.yx(3) = .5*(xi2Interpolant.y.Values(idx1(2)-1,idx1(1))*tmp(3) ...
                - xi2Interpolant.y.Values(idx1(2)-1,idx1(1)-2)*tmp(8))/deltaX;
            dxi2Cts.yy(3) = .5*(xi2Interpolant.y.Values(idx1(2),idx1(1)-1)*tmp(1) ...
                - xi2Interpolant.y.Values(idx1(2)-2,idx1(1)-1)*tmp(9))/deltaY;
            % Position 4
            dxi2Cts.xx(4) = .5*(xi2Interpolant.x.Values(idx1(2)-1,idx1(1)+1)*tmp(11) ...
                - xi2Interpolant.x.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaX;
            dxi2Cts.xy(4) = .5*(xi2Interpolant.x.Values(idx1(2),idx1(1)) ...
                - xi2Interpolant.x.Values(idx1(2)-2,idx1(1))*tmp(10))/deltaY;
            dxi2Cts.yx(4) = .5*(xi2Interpolant.y.Values(idx1(2)-1,idx1(1)+1)*tmp(11) ...
                - xi2Interpolant.y.Values(idx1(2)-1,idx1(1)-1)*tmp(2))/deltaX;
            dxi2Cts.yy(4) = .5*(xi2Interpolant.y.Values(idx1(2),idx1(1)) ...
                - xi2Interpolant.y.Values(idx1(2)-2,idx1(1))*tmp(10))/deltaY;
            
            ctsDxi2Interpolant.xx = griddedInterpolant({positionY,positionX},...
                [dxi2Cts.xx(3) dxi2Cts.xx(4); dxi2Cts.xx(2) dxi2Cts.xx(1)]);
            dxi2.xx(iPoint) = ctsDxi2Interpolant.xx(position(iPoint,2),position(iPoint,1));
            ctsDxi2Interpolant.xy = griddedInterpolant({positionY,positionX},...
                [dxi2Cts.xy(3) dxi2Cts.xy(4); dxi2Cts.xy(2) dxi2Cts.xy(1)]);
            dxi2.xy(iPoint) = ctsDxi2Interpolant.xy(position(iPoint,2),position(iPoint,1));
            ctsDxi2Interpolant.yx = griddedInterpolant({positionY,positionX},...
                [dxi2Cts.yx(3) dxi2Cts.yx(4); dxi2Cts.yx(2) dxi2Cts.yx(1)]);
            dxi2.yx(iPoint) = ctsDxi2Interpolant.yx(position(iPoint,2),position(iPoint,1));
            ctsDxi2Interpolant.yy = griddedInterpolant({positionY,positionX},...
                [dxi2Cts.yy(3) dxi2Cts.yy(4); dxi2Cts.yy(2) dxi2Cts.yy(1)]);
            dxi2.yy(iPoint) = ctsDxi2Interpolant.yy(position(iPoint,2),position(iPoint,1));
            
        else
            
            eta(iPoint,1) = etaInterpolant.x(position(iPoint,2),...
                position(iPoint,1));
            eta(iPoint,2) = etaInterpolant.y(position(iPoint,2),...
                position(iPoint,1));
            
            xi1(iPoint,1) = xi1Interpolant.x(position(iPoint,2),position(iPoint,1));
            xi1(iPoint,2) = xi1Interpolant.y(position(iPoint,2),position(iPoint,1));
            
            xi2(iPoint,1) = xi2Interpolant.x(position(iPoint,2),position(iPoint,1));
            xi2(iPoint,2) = xi2Interpolant.y(position(iPoint,2),position(iPoint,1));
            
            dxi1.xx(iPoint) = dxi1Interpolant.xx(position(iPoint,2),...
                position(iPoint,1));
            dxi1.xy(iPoint) = dxi1Interpolant.xy(position(iPoint,2),...
                position(iPoint,1));
            dxi1.yx(iPoint) = dxi1Interpolant.yx(position(iPoint,2),...
                position(iPoint,1));
            dxi1.yy(iPoint) = dxi1Interpolant.yy(position(iPoint,2),...
                position(iPoint,1));
            
            dxi2.xx(iPoint) = dxi2Interpolant.xx(position(iPoint,2),...
                position(iPoint,1));
            dxi2.xy(iPoint) = dxi2Interpolant.xy(position(iPoint,2),...
                position(iPoint,1));
            dxi2.yx(iPoint) = dxi2Interpolant.yx(position(iPoint,2),...
                position(iPoint,1));
            dxi2.yy(iPoint) = dxi2Interpolant.yy(position(iPoint,2),...
                position(iPoint,1));
            
        end
    end
end

kappa(:,1) = dot([dot(transpose(cat(1,dxi1.xx,dxi1.xy)),xi1,2) ...
    dot(transpose(cat(1,dxi1.yx,dxi1.yy)),xi1,2)],xi2,2);
kappa(:,2) = dot([dot(transpose(cat(1,dxi2.xx,dxi2.xy)),xi2,2) ...
    dot(transpose(cat(1,dxi2.yx,dxi2.yy)),xi2,2)],xi1,2);

switch posOrNeg
    case 'pos'
        geodesicDeviation = abs(1 - alpha) + abs(-1./beta.*dot(dAlpha,eta,2) + ...
            alpha.*kappa(:,1) - beta.*kappa(:,2) + (l1./l2 - 1).*kappa(:,1) - ...
            .5./l2.*dot(dl1,xi2,2));
    case 'neg'
        geodesicDeviation = abs(1 - alpha) + abs(1./beta.*dot(dAlpha,eta,2) + ...
            alpha.*kappa(:,1) + beta.*kappa(:,2) + (l1./l2 - 1).*kappa(:,1) - ...
            .5./l2.*dot(dl1,xi2,2));
    otherwise
        error('posOrNeg invalid')
end

function averageGeodesicDeviation = average_geodesic_deviation_individual(...
    shearlinePosition,geodesicDeviation)

if size(shearlinePosition,1) == 1
    averageGeodesicDeviation = geodesicDeviation;
else
    dPosition = diff(shearlinePosition);
    arcLength = [0; cumsum(hypot(dPosition(:,1),dPosition(:,2)))];
    averageGeodesicDeviation = trapz(arcLength,geodesicDeviation)/...
        arcLength(end);
end

function idx = position2index(position,domain,resolution)

deltaX = diff(domain(1,:))/(double(resolution(1)) - 1);
xMin = domain(1,1);

deltaY = diff(domain(2,:))/(double(resolution(2)) - 1);
yMin = domain(2,1);

% Index of upper-right corner
idxX = ceil((position(1) - xMin)/deltaX) + 1;
idxY = ceil((position(2) - yMin)/deltaY) + 1;
idx = [idxX,idxY];

% iBE is 0 if position is not in a boundary element. Otherwise iBE is set
% as follows:
%
% 2    5    1
%  ┌───────┐
%  │       │
% 6│       │8
%  │       │
%  └───────┘
% 3    7    4
% 
function iBE = is_boundary_element(position,domain,resolution)

idx = position2index(position,domain,resolution);

if idx(1) == resolution(1)
    if idx(2) == resolution(2)
        iBE = 1;
        return
    end
    if idx(2) == 2
        iBE = 4;
        return
    end
    iBE = 8;
    return
end

if idx(1) == 2
    if idx(2) == resolution(2)
        iBE = 2;
        return
    end
    if idx(2) == 2
        iBE = 3;
        return
    end
    iBE = 6;
    return
end

if idx(2) == resolution(2)
    iBE = 5;
    return
end

if idx(2) == 2
    iBE = 7;
    return
end

iBE = 0;