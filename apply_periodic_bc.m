function [positionPeriodic,varargout] = apply_periodic_bc(position,periodicBc,domain)

nargoutchk(1,2)

if periodicBc(1)
    positionPeriodic(:,1) = mod(position(:,1),diff(domain(1,:))) + domain(1,1);
    positionPeriodic(:,2) = position(:,2);
    if nargout > 1
        idx = find(diff(floor(position(:,1)/diff(domain(1,:)))));
        varargout{1} = idx;
    end
end

if periodicBc(2)
    warning('Periodic BC in y not programmed')
end
