function plot_ftle(axes,flow)

ftle = compute_ftle(flow.cgEigenvalue(:,2),abs(diff(flow.timespan)));
ftle = reshape(ftle,fliplr(flow.resolution));
h = pcolor(axes,...
    linspace(flow.domain(1,1),flow.domain(1,2),flow.resolution(1)),...
    linspace(flow.domain(2,1),flow.domain(2,2),flow.resolution(2)),...
    ftle);
set(h,'tag','ftle')
shading(axes,'interp')
