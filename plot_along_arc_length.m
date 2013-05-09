% plot_along_arc_length Plot value along curve as a function of curve length.
%
% SYNTAX
% hPlot = plot_along_arc_length(hAxes,position,value);

function hPlot = plot_along_arc_length(hAxes,position,value)

xD = diff(position(:,1));
yD = diff(position(:,2));
arcLength = [0; cumsum(sqrt(xD.^2 + yD.^2))];

hPlot = plot(hAxes,arcLength,value);

