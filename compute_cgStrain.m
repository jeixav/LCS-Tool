function cgStrain = compute_cgStrain(finalPosition,flow,...
    auxiliaryGridRelativeDelta)
%compute_cgStrain   Compute Cauchy-Green strain
%   cgStrain = compute_cgStrain(finalPosition,delta)
%
%   finalPosition gives positions grouped in rows. In 2 dimensions, each
%   row has 8 columns with the format:
%      (x+delta, y) (x-delta, y) (x, y+delta) (x, y-delta)
%
%   cgStrain gives the Cauchy-Green strain tensor. In 2 dimensions it has
%   the format: c_11 c_12 c_22
% On main grid:
% finalPosition = [x1 y1; x2 y2; ... xN yN]
%
% On auxiliary grid:
% finalPosition = [x1+deltaX y1 x1-deltaX y1 x1 y1+deltaY x1 y1-deltaY;
%                  x2+deltaX y2 x2-deltaX y2 x2 y2+deltaY x2 y2-deltaY;
%                                          ...
%                  xN+deltaX yN xN-deltaX yN xN yN+deltaY xN yN-deltaY]

cgStrain = nan(size(finalPosition,1),3);

switch size(finalPosition,2)
    case 2 % Main grid
        
        finalX = reshape(finalPosition(:,1),fliplr(flow.resolution));
        finalY = reshape(finalPosition(:,2),fliplr(flow.resolution));
        
        deltaX = (flow.domain(1,2) - flow.domain(1,1))...
            /double(flow.resolution(1));
        deltaY = (flow.domain(2,2) - flow.domain(2,1))...
            /double(flow.resolution(2));
        
        [gradF11,gradF12] = gradient(finalX,deltaX,deltaY);
        [gradF21,gradF22] = gradient(finalY,deltaX,deltaY);
        
        gradF11 = reshape(gradF11,prod(double(flow.resolution)),1);
        gradF12 = reshape(gradF12,prod(double(flow.resolution)),1);
        gradF21 = reshape(gradF21,prod(double(flow.resolution)),1);
        gradF22 = reshape(gradF22,prod(double(flow.resolution)),1);
        
    case 8 % Auxiliary grid

        finalX = finalPosition(:,1:2:7);
        finalY = finalPosition(:,2:2:8);
        
        deltaX = (flow.domain(1,2) - flow.domain(1,1))...
            /double(flow.resolution(1))*auxiliaryGridRelativeDelta;
        deltaY = (flow.domain(2,2) - flow.domain(2,1))...
            /double(flow.resolution(2))*auxiliaryGridRelativeDelta;
        if deltaX ~= deltaY
            warning('compute_cgStrain:unequalDelta',...
                ['deltaX and deltaY unequal. deltaX = ',...
                num2str(deltaX),', deltaY = ',num2str(deltaY),'.'])
        end
        
        gradF11 = (finalX(:,1) - finalX(:,2))/(2*deltaX);
        gradF12 = (finalX(:,3) - finalX(:,4))/(2*deltaY);
        gradF21 = (finalY(:,1) - finalY(:,2))/(2*deltaX);
        gradF22 = (finalY(:,3) - finalY(:,4))/(2*deltaY);

    otherwise
        error('Number of columns in finalPosition incorrect.')
end

% cgStrain = [c11 c12 c22]
cgStrain(:,1) = gradF11.^2 + gradF21.^2;
cgStrain(:,2) = gradF11.*gradF12 + gradF21.*gradF22;
cgStrain(:,3) = gradF12.^2 + gradF22.^2;

end
