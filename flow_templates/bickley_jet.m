% bickley_jet Defines a bickley jet.
%
% DESCRIPTION
% bickley_jet(p)
% Time periodic forcing case from doi:10.1016/j.physd.2012.06.012.
% 
% REFERENCES
% doi:10.1016/j.physd.2012.06.012
% doi:10.1063/1.3271342

function bickleyJet = bickley_jet_coupled(perturbationCase,p)

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
    case {2,3}
        % Time-aperiodic psi1, cases 2 and 3 on page 1691 of
        % doi:10.1016/j.physd.2012.06.012.
        
        % Duffing oscillator
        % FIXME Value copied from Beron-Vera's lcsgeo_bickleyduffing_u1
        beronVeraNT = 5;
        beronVeraW = 5;
        phiTimespan = [0,beronVeraNT*beronVeraW];
        
        phiInitial = [0,0];
        phiSol = ode45(@d_phi,phiTimespan,phiInitial);
        
        timeResolution = 1e3;
        phi2 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),2);
        
        % Computational optimization -- solve forcing function once for
        % entire timespan and use interpolation when integrating forced
        % flow.
        phi2Int = griddedInterpolant(linspace(phiTimespan(1),phiTimespan(2),timeResolution),phi2);
        
        beronVeraT = max(2*pi./abs([sigma1 sigma2]));
        
        % Find maximum value of phi
        % FIXME Sufficently large maxSamples found heuristically.
        % maxSamples = 1e3;
        phi2Max = max(phi2);
        
        amplitudeCorrection = .015;
        beronVeraMagicScaleAmp = 1.75*amplitudeCorrection;
        beronVeraMagicScaleTime = beronVeraT/beronVeraW;
        
        f1 = @(t)beronVeraMagicScaleAmp*phi2Int(t/beronVeraMagicScaleTime)*phi2Max;
        f2 = @(t)f1(t);
        
    otherwise
        error('Invalid perturbation case selected')
end

    function derivative_ = derivative(t,x,useEoV)
        
        t = timescale*t;
        
        if useEoV
            idx1 = 1:6:size(x,1)-5;
            idx2 = 2:6:size(x,1)-4;
        else
            idx1 = 1:2:size(x,1)-1;
            idx2 = 2:2:size(x,1);
        end
        
        derivative_ = nan(size(x));
        
        % u
        derivative_(idx1) = -p.c3 - p.u*(tanh(x(idx2)).^2 - 1) + 2*p.u*sech(x(idx2)).^2.*tanh(x(idx2)).*(real(p.epsilon1*f1(t)*exp(k(1)*lengthscaleX*x(idx1)*1i) + p.epsilon2*f2(t)*exp(k(2)*lengthscaleX*x(idx1)*1i)) + p.epsilon3*cos(.5*k(3)*lengthscaleX*x(idx1)));
        derivative_(idx1) = derivative_(idx1)/lengthscaleX*timescale;
        
        % v
        derivative_(idx2) = -p.u*p.lengthY*sech(x(idx2)).^2.*(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*lengthscaleX*x(idx1)*1i) + p.epsilon2*f2(t)*k(2)*exp(k(2)*lengthscaleX*x(idx1)*1i)) + .5*p.epsilon3*k(3)*sin(.5*k(3)*lengthscaleX*x(idx1)));
        derivative_(idx2) = derivative_(idx2)/lengthscaleY*timescale;
        
        if useEoV
            idx3 = 3:6:size(x,1)-3;
            idx4 = 4:6:size(x,1)-2;
            idx5 = 5:6:size(x,1)-1;
            idx6 = 6:6:size(x,1);
            
            % dux
            dux = -2*p.u*sech(x(idx2)).^2.*tanh(x(idx2)).*(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*lengthscaleX*x(idx1)*1i)) + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*lengthscaleX*x(idx1)*1i)) + .5*p.epsilon3*k(3)*sin(.5*k(3)*lengthscaleX*x(idx1)));
            dux = dux*timescale;
            
            % duy
            duy = -2*p.u/p.lengthY*sech(x(idx2)).^2.*((3*tanh(x(idx2)).^2 - 1).*(real(p.epsilon1*f1(t)*exp(k(1)*lengthscaleX*x(idx1)*1i)) + real(p.epsilon2*f2(t)*exp(k(2)*lengthscaleX*x(idx1)*1i)) + p.epsilon3*cos(.5*k(3)*lengthscaleX*x(idx1))) + tanh(x(idx2)));
            duy = duy/lengthscaleX*lengthscaleY*timescale;
            
            % dvx
            dvx = -p.u*p.lengthY*sech(x(idx2)).^2.*(real(p.epsilon1*f1(t)*k(1)^2*exp(k(1)*lengthscaleX*x(idx1)*1i)) + real(p.epsilon2*f2(t)*k(2)^2*exp(k(2)*lengthscaleX*x(idx1)*1i)) + .5^2*p.epsilon3*k(3)^2*cos(.5*k(3)*lengthscaleX*x(idx1)));
            dvx = dvx/lengthscaleY*lengthscaleX*timescale;
            
            % dvy
            dvy = 2*p.u*sech(x(idx2)).^2.*tanh(x(idx2)).*(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*lengthscaleX*x(idx1)*1i)) + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*lengthscaleX*x(idx1)*1i)) + .5*p.epsilon3*k(3)*sin(.5*k(3)*lengthscaleX*x(idx1)));
            dvy = dvy*timescale;
            
            % Perform matrix multiplication manually
            derivative_(idx3) = dux.*x(idx3) + duy.*x(idx5);
            derivative_(idx4) = dux.*x(idx4) + duy.*x(idx6);
            derivative_(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
            derivative_(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
        end
    end

bickleyJet.flow.derivative = @(t,x,useEoV)derivative(t,x,useEoV);

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true,false];
bickleyJet.flow.coupledIntegration = true;
bickleyJet.flow = set_flow_domain([0,2*pi;[-1,1]*2.25],bickleyJet.flow);
bickleyJet.flow = set_flow_timespan([0,20],bickleyJet.flow);
bickleyJet.flow = set_flow_resolution([4096,1637],bickleyJet.flow);

bickleyJet.strainline = set_strainline_resolution([2,1]*10);
bickleyJet.strainline = set_strainline_max_length(20,bickleyJet.strainline);
bickleyJet.strainline = set_strainline_geodesic_deviation_tol(1e-5,bickleyJet.strainline);
bickleyJet.strainline = set_strainline_length_tol(4,bickleyJet.strainline);
bickleyJet.strainline.filteringMethod = 'superminimize';
bickleyJet.strainline.filteringParameters = struct('distance',.5,'resolution',uint64([5,2]));

bickleyJet.shearline = set_shearline_resolution([2,1]*10);
bickleyJet.shearline = set_shearline_max_length(10,bickleyJet.shearline);
bickleyJet.shearline = set_shearline_average_geodesic_deviation_tol([inf,inf],bickleyJet.shearline);

end

% Forced-damped Duffing oscillator used with aperiodic forcing
function dPhi = d_phi(tau,phi)
dPhi(2,1) = nan;

dPhi(1) = phi(2);
dPhi(2) = -.1*phi(2) - phi(1)^3 + 11*cos(tau);
end
