function [positionPeriodic,varargout] = apply_periodic_bc(position,periodicBc,domain)

nargoutchk(1,2)

positionPeriodic = position;

for m = 1:numel(periodicBc)
    if periodicBc(m)
        positionPeriodic(:,m) = mod(position(:,1) - domain(m,1),diff(domain(m,:))) + domain(m,1);
        if nargout > 1
            if m == 2
                error('Code below not validated for periodicBc(2) == true')
            end
            idx = find(diff(floor(position(:,m)/diff(domain(m,:)))));
            varargout{1} = idx;
        end
    end
end
