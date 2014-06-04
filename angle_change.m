% angle_change Angular change of vector field between adjacent grid points.
%
% SYNTAX
% [thetaX,thetaY] = angle_change(vector)
% [thetaX,thetaY,thetaMax] = angle_change(vector)

function [thetaX,thetaY,varargout] = angle_change(vector)

% Compare x to x + DeltaX
a = vector(:,:,1:end-1);
b = vector(:,:,2:end);

normA = norm_array(a);
normB = norm_array(b);

% Formula used is: dot(a,b) = norm(a)*norm(b)*cos(theta)
thetaX = acos(squeeze(dot(a,b,1))./normA./normB);

% Compare y to y + DeltaY
a = vector(:,1:end-1,:);
b = vector(:,2:end,:);

normA = norm_array(a);
normB = norm_array(b);

thetaY = acos(squeeze(dot(a,b,1))./normA./normB);

if nargout == 3
    % Calculate the maximum angular change at a point.
    thetaMax = nan(size(thetaY,1),size(thetaX,2));
    for m = 1:size(thetaY,1)
        for n = 1:size(thetaX,2)
            thetaMax(m,n) = max([thetaX(m,n),thetaX(m+1,n),thetaY(m,n),thetaY(m,n+1)]);
        end
    end
    varargout{1} = thetaMax;
end

% Norm of 2-by-m-by-n array.
function normArray = norm_array(array)

normArray = squeeze(sqrt(array(1,:,:).^2 + array(2,:,:).^2));
