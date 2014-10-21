% derivative Ocean data velocity
%
% SYNTAX
% derivative_ = derivative(time,position,VX_interpolant,VY_interpolant)
%
% INPUT ARGUMENTS
% time: scalar
% position: [x1;y1;x2;y2;...;xN;yN]
% VX_interpolant: griddedInterpolant for x-component of velocity
% VY_interpolant: griddedInterpolant for y-component of velocity
%
% OUTPUT ARGUMENT
% derivative_: [xVelocity1;yVelocity1;xVelocity2;yVelocity2;...;xVelocityN;yVelocityN]

function  derivative_ = derivative(time,position,VX_interpolant,VY_interpolant)

nPosition = numel(position)/2;
derivative_ = nan(nPosition*2,1);

% x-positions
derivative_(1:2:end-1) = VX_interpolant(time*ones(nPosition,1),position(2:2:end),position(1:2:end-1));

% y-positions
derivative_(2:2:end) = VY_interpolant(time*ones(nPosition,1),position(2:2:end),position(1:2:end-1));
