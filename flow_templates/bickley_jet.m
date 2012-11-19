% Describes a bickley jet.
%
% References: doi:10.1016/j.physd.2012.06.012, doi:10.1063/1.3271342

function bickleyJet = bickley_jet_non_dim(perturbationCase,p)

narginchk(1,2)

defaultP.u = 62.66;
defaultP.c2 = .205*defaultP.u;
defaultP.c3 = .461*defaultP.u;
defaultP.a = 6.371e6;
defaultP.lengthX = pi*defaultP.a;
defaultP.lengthY = 1.77e6;

defaultP.epsilon1 = .075;
defaultP.epsilon2 = .4;
defaultP.epsilon3 = .3;

if nargin < 2
    p = defaultP;
else
    fieldnamesDefaultP = fieldnames(defaultP);
    for i = 1:numel(fieldnamesDefaultP)
        if ~isfield(p,fieldnamesDefaultP{i}) || ...
                isempty(p.(fieldnamesDefaultP{i}))
            p.(fieldnamesDefaultP{i}) = defaultP.(fieldnamesDefaultP{i});
        end
    end
end

k = @(n) 2*n*pi/p.lengthX;

sigma1 = .5*k(2)*(p.c2 - p.c3);
sigma2 = 2*sigma1;

timescale = p.lengthX/p.u;
lengthscaleX = p.lengthX/2/pi;
lengthscaleY = p.lengthY;

switch perturbationCase
    case {1,3}
        p.epsilon1 = p.epsilon1/10;
        p.epsilon2 = p.epsilon2/10;
end

switch perturbationCase
    case 1
        % Time-periodic psi1, case 1 on page 1691 of
        % doi:10.1016/j.physd.2012.06.012.
 
        f1 = @(t)exp(1i*sigma1*t);
        f2 = @(t)exp(1i*sigma2*t);

        u = @(t,x)-p.c3 - p.u*(tanh(x(2)/p.lengthY)^2 - 1) ...
            + 2*p.u*sech(x(2)/p.lengthY)^2 ...
            *tanh(x(2)/p.lengthY) ...
            *(real(p.epsilon1*f1(t)*exp(k(1)*x(1)*1i) ...
            + p.epsilon2*f2(t)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*cos(k(3)*x(1)));
        v = @(t,x)-p.u*p.lengthY*sech(x(2)/p.lengthY)^2 ...
            *(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*x(1)*1i) ...
            + p.epsilon2*f2(t)*k(2)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*x(1)));
        
        lengthscale = [lengthscaleX; lengthscaleY];
        bickleyJet.flow.derivative = @(t,x)...
            [2*u(timescale*t,lengthscale.*x)/lengthscaleX; ...
            v(timescale*t,lengthscale.*x)/lengthscaleY]*timescale;

        dux = @(t,x)-2*p.u ...
            *sech(x(2)/p.lengthY)^2*tanh(x(2)/p.lengthY) ...
            *(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*x(1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*x(1)));
        duy = @(t,x)-2*p.u/p.lengthY*sech(x(2)/p.lengthY)^2 ...
            *((3*tanh(x(2)/p.lengthY)^2 - 1) ...
            *(real(p.epsilon1*f1(t)*exp(k(1)*x(1)*1i)) ...
            + real(p.epsilon2*f2(t)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*cos(k(3)*x(1))) + tanh(x(2)/p.lengthY));
        dvx = @(t,x)-p.u*p.lengthY*sech(x(2)/p.lengthY)^2 ...
            *(real(p.epsilon1*f1(t)*k(1)^2*exp(k(1)*x(1)*1i)) ...
            + real(p.epsilon2*f2(t)*k(2)^2*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)^2*cos(k(3)*x(1)));
        dvy = @(t,x)2*p.u*sech(x(2)/p.lengthY)^2*tanh(x(2)/p.lengthY) ...
            *(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*x(1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*x(1)));
        
        bickleyJet.flow.dDerivative = ...
            @(t,x)[2*dux(timescale*t,lengthscale.*x)/lengthscaleX*lengthscaleX ...
            2*duy(timescale*t,lengthscale.*x)/lengthscaleX*lengthscaleY; ...
            dvx(timescale*t,lengthscale.*x)/lengthscaleY*lengthscaleX ...
            dvy(timescale*t,lengthscale.*x)/lengthscaleY*lengthscaleY]*timescale;
    
    case {2,3}
        % Time-aperiodic psi1, case 2 on page 1691 of
        % doi:10.1016/j.physd.2012.06.012

        % Duffing oscillator
        % FIXME Value copied from Javier's lcsgeo_bickleyduffing_u1
        beronVeraNT = 5;
        beronVeraW = 5;
        phiTimespan = [0 beronVeraNT*beronVeraW];

        phiInitial = [0 0];
        phiSol = ode45(@d_phi,phiTimespan,phiInitial);
        
        timeResolution = 1e3;
        phi2 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),...
            timeResolution),2);
        
        % Computational optimization -- solve forcing function once for
        % entire timespan and use interpolation when integrating forced
        % flow.
        phi2Int = griddedInterpolant(linspace(phiTimespan(1),...
            phiTimespan(2),timeResolution),phi2);
        
        beronVeraT = max(2*pi./abs([sigma1 sigma2]));
        
        % Find maximum value of phi
        % FIXME Sufficently large maxSamples found heuristically.
        % maxSamples = 1e3;
        phi2Max = max(phi2);
        
        amplitudeCorrection = .015;
        beronVeraMagicScaleAmp = 1.75*amplitudeCorrection;
        beronVeraMagicScaleTime = beronVeraT/beronVeraW;
        f1 = @(t)beronVeraMagicScaleAmp*phi2Int(t ...
            /beronVeraMagicScaleTime)*phi2Max;
        f2 = @(t)f1(t);
        
        u = @(t,x)-p.c3 - p.u*(tanh(x(2)/p.lengthY)^2 - 1) ...
            + 2*p.u*sech(x(2)/p.lengthY)^2 ...
            *tanh(x(2)/p.lengthY) ...
            *(real(p.epsilon1*f1(t)*exp(k(1)*x(1)*1i)) ...
            + real(p.epsilon2*f2(t)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*cos(k(3)*x(1)));
        v = @(t,x)-p.u*p.lengthY*sech(x(2)/p.lengthY)^2 ...
            *(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*x(1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*x(1)));
        
        bickleyJet.flow.derivative = @(t,x)...
            [u(t*timescale,[x(1)*lengthscaleX x(2)*lengthscaleY]); ...
            v(t*timescale,[x(1)*lengthscaleX x(2)*lengthscaleY])];

        % FIXME These need not be defined unless using the equation of
        % variation method
        dux = @(t,x)-2*p.u ...
            *sech(x(2)/p.lengthY)^2*tanh(x(2)/p.lengthY) ...
            *(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*x(1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*x(1)));
        duy = @(t,x)-2*p.u/p.lengthY*sech(x(2)/p.lengthY)^2 ...
            *((3*tanh(x(2)/p.lengthY)^2 - 1) ...
            *(real(p.epsilon1*f1(t)*exp(k(1)*x(1)*1i)) ...
            + real(p.epsilon2*f2(t)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*cos(k(3)*x(1))) + tanh(x(2)/p.lengthY));
        dvx = @(t,x)-p.u*p.lengthY*sech(x(2)/p.lengthY)^2 ...
            *(real(p.epsilon1*f1(t)*k(1)^2*exp(k(1)*x(1)*1i)) ...
            + real(p.epsilon2*f2(t)*k(2)^2*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)^2*cos(k(3)*x(1)));
        dvy = @(t,x)2*p.u*sech(x(2)/p.lengthY)^2*tanh(x(2)/p.lengthY) ...
            *(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*x(1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*x(1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*x(1)));
        
        bickleyJet.flow.dDerivative = @(t,x)[dux(t,x) duy(t,x); ...
            dvx(t,x) dvy(t,x)];

    otherwise
        error('Invalid perturbation case selected')
end

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true false];

bickleyJet.flow = set_flow_domain([0 2*pi; [-1 1]*2],bickleyJet.flow);

finalTime = 4;
bickleyJet.flow = set_flow_timespan([0 finalTime],bickleyJet.flow);

bickleyJet.flow = set_flow_resolution(100,bickleyJet.flow);

bickleyJet.strainline = set_strainline_resolution(uint64([2 1]*10));
bickleyJet.strainline = set_strainline_max_length(20,bickleyJet.strainline);
bickleyJet.strainline = set_strainline_geodesic_deviation_tol(inf,...
    bickleyJet.strainline);
bickleyJet.strainline = set_strainline_length_tol(0,bickleyJet.strainline);
bickleyJet.strainline.filteringMethod = 'superminimize';
bickleyJet.strainline.filteringParameters = struct('distance',0,...
    'resolution',uint64([5 2]));

bickleyJet.shearline = set_shearline_resolution(uint64([2 1]*5));
bickleyJet.shearline = set_shearline_max_length(.1,bickleyJet.shearline);
bickleyJet.shearline = set_shearline_average_geodesic_deviation_tol(...
    [inf inf],bickleyJet.shearline);

% Forced-damped Duffing oscillator used with aperiodic forcing
function dPhi = d_phi(tau,phi)
dPhi(2,1) = nan;

dPhi(1) = phi(2);
dPhi(2) = -.1*phi(2) - phi(1)^3 + 11*cos(tau);
