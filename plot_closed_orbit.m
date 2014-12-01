function h = plot_closed_orbit(hAxes,closedOrbit)

nLambda = size(closedOrbit,1);
nPoincareSection = size(closedOrbit,2);

h = cell(size(closedOrbit));

for iPs = 1:nPoincareSection
    for iL = 1:nLambda
        nOrbit = numel(closedOrbit{iL,iPs});
        h{iL,iPs} = gobjects(1,nOrbit);
        for iO = 1:nOrbit
            x = closedOrbit{iL,iPs}{iO}(:,1);
            y = closedOrbit{iL,iPs}{iO}(:,2);
            % Skip closedOrbit elements with only NaNs
            if ~all(isnan([x,y]))
                h{iL,iPs}(iO) = plot(hAxes,x,y);
            end
        end
    end
end

h = [h{:}];
if verLessThan('matlab','8.4')
    h = h(find(h)); %#ok<FNDSB>
else
    h = h(isgraphics(h));
end
