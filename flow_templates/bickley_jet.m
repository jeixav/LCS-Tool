% bickley_jet Defines a bickley jet.
%
% DESCRIPTION
% bickley_jet(p)
% Time periodic forcing case from doi:10.1016/j.physd.2012.06.012.
% 
% REFERENCES
% doi:10.1016/j.physd.2012.06.012
% doi:10.1063/1.3271342

function bickleyJet = bickley_jet(perturbationCase,p)

narginchk(1,2)

defaultP.u = 62.66;
defaultP.c2 = .205*defaultP.u;
defaultP.c3 = .461*defaultP.u;

defaultP.a = earthRadius;

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

k = @(n)2*n*pi/p.lengthX;

sigma1 = .5*k(2)*(p.c2 - p.c3);
sigma2 = 2*sigma1;

% timescale = p.lengthX/p.u;
% lengthscaleX = p.lengthX/2/pi;
% lengthscaleY = p.lengthY;

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
        
        timeResolution = 1e5;
        phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
                
        % Computational optimization -- solve forcing function once for
        % entire timespan and use interpolation when integrating forced
        % flow.
        phi1Int = griddedInterpolant(linspace(phiTimespan(1),phiTimespan(2),timeResolution),phi1);
        
        beronVeraT = max(2*pi./abs([sigma1,sigma2]));
        
        % Find maximum value of phi
        % FIXME Sufficently large maxSamples found heuristically.
        % maxSamples = 1e3;
        phi1Max = max(phi1);
        
        amplitudeCorrection = .015;
        beronVeraMagicScaleAmp = 1.75*amplitudeCorrection;
        beronVeraMagicScaleTime = beronVeraT/beronVeraW;
        
        f1 = @(t)beronVeraMagicScaleAmp*phi1Int(t/beronVeraMagicScaleTime)*phi1Max;
        f2 = f1;
    otherwise
        error('Invalid perturbation case selected')
end

    function derivative_ = derivative(t,x,useEoV)
        
        if useEoV
            idx1 = 1:6:size(x,1)-5;
            idx2 = 2:6:size(x,1)-4;
        else
            idx1 = 1:2:size(x,1)-1;
            idx2 = 2:2:size(x,1);
        end
        
        derivative_ = nan(size(x));
        
        % u
        derivative_(idx1) = (cosh(x(idx2)./p.lengthY).*(p.u + p.u.*p.epsilon3.*cos(k(3).*x(idx1)).*sinh(x(idx2)./p.lengthY)) + 2.*p.u.*real(p.epsilon1.*f1(t).*exp(k(1).*x(idx1).*1i)).*sinh(x(idx2)./p.lengthY) + 2.*p.u.*real(p.epsilon2.*f2(t).*exp(k(2).*x(idx1).*1i)).*sinh(x(idx2)./p.lengthY))./cosh(x(idx2)./p.lengthY).^3 - p.c3;
        
        % v
        derivative_(idx2) = - (p.lengthY.*p.u.*(imag(p.epsilon1.*f1(t).*k(1).*exp(k(1).*x(idx1).*1i)) + imag(p.epsilon2.*f2(t).*k(2).*exp(k(2).*x(idx1).*1i))))./cosh(x(idx2)./p.lengthY).^2 - (p.lengthY.*p.u.*p.epsilon3.*k(3).*sin(k(3).*x(idx1)))./cosh(x(idx2)./p.lengthY);
        
        if useEoV
            idx3 = 3:6:size(x,1)-3;
            idx4 = 4:6:size(x,1)-2;
            idx5 = 5:6:size(x,1)-1;
            idx6 = 6:6:size(x,1);
            
            dux = -(2.*p.u.*imag(p.epsilon1.*f1(t).*k(1).*exp(k(1).*x(idx1).*1i)).*sinh(x(idx2)./p.lengthY) + 2.*p.u.*imag(p.epsilon2.*f2(t).*k(2).*exp(k(2).*x(idx1).*1i)).*sinh(x(idx2)./p.lengthY) + p.u.*p.epsilon3.*k(3).*sin(k(3).*x(idx1)).*cosh(x(idx2)./p.lengthY).*sinh(x(idx2)./p.lengthY))./cosh(x(idx2)./p.lengthY).^3;
            duy = ((sinh(x(idx2)./p.lengthY).*(p.u + p.u.*p.epsilon3.*cos(k(3).*x(idx1)).*sinh(x(idx2)./p.lengthY)))./p.lengthY + (2.*p.u.*real(p.epsilon1.*f1(t).*exp(k(1).*x(idx1).*1i)).*cosh(x(idx2)./p.lengthY))./p.lengthY + (2.*p.u.*real(p.epsilon2.*f2(t).*exp(k(2).*x(idx1).*1i)).*cosh(x(idx2)./p.lengthY))./p.lengthY + (p.u.*p.epsilon3.*cos(k(3).*x(idx1)).*cosh(x(idx2)./p.lengthY).^2)./p.lengthY)./cosh(x(idx2)./p.lengthY).^3 - (3.*sinh(x(idx2)./p.lengthY).*(cosh(x(idx2)./p.lengthY).*(p.u + p.u.*p.epsilon3.*cos(k(3).*x(idx1)).*sinh(x(idx2)./p.lengthY)) + 2.*p.u.*real(p.epsilon1.*f1(t).*exp(k(1).*x(idx1).*1i)).*sinh(x(idx2)./p.lengthY) + 2.*p.u.*real(p.epsilon2.*f2(t).*exp(k(2).*x(idx1).*1i)).*sinh(x(idx2)./p.lengthY)))./(p.lengthY.*cosh(x(idx2)./p.lengthY).^4);
            dvx = - (p.lengthY.*p.u.*(real(p.epsilon1.*f1(t).*k(1).^2.*exp(k(1).*x(idx1).*1i)) + real(p.epsilon2.*f2(t).*k(2).^2.*exp(k(2).*x(idx1).*1i))))./cosh(x(idx2)./p.lengthY).^2 - (p.lengthY.*p.u.*p.epsilon3.*k(3).^2.*cos(k(3).*x(idx1)))./cosh(x(idx2)./p.lengthY);
            dvy = (2.*p.u.*sinh(x(idx2)./p.lengthY).*(imag(p.epsilon1.*f1(t).*k(1).*exp(k(1).*x(idx1).*1i)) + imag(p.epsilon2.*f2(t).*k(2).*exp(k(2).*x(idx1).*1i))))./cosh(x(idx2)./p.lengthY).^3 + (p.u.*p.epsilon3.*k(3).*sin(k(3).*x(idx1)).*sinh(x(idx2)./p.lengthY))./cosh(x(idx2)./p.lengthY).^2;
            
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
bickleyJet.flow.coupledIntegration = 1e5;
bickleyJet.flow = set_flow_domain([0,p.lengthX;[-1,1]*2.2599*p.lengthY],bickleyJet.flow);
bickleyJet.flow = set_flow_timespan([0,20*p.lengthX/p.u],bickleyJet.flow);
bickleyJet.flow = set_flow_resolution(50,bickleyJet.flow);
bickleyJet.flow.cgStrainMethod.name = 'finiteDifference';
bickleyJet.flow.cgStrainMethod.eigenvalueFromMainGrid = false;

bickleyJet.strainline = set_strainline_resolution([2,1]*10);
bickleyJet.strainline = set_strainline_max_length(1e6,bickleyJet.strainline);
bickleyJet.strainline = set_strainline_geodesic_deviation_tol(1e-5,bickleyJet.strainline);
bickleyJet.strainline = set_strainline_length_tol(4,bickleyJet.strainline);
bickleyJet.strainline.filteringMethod = 'superminimize';
bickleyJet.strainline.filteringParameters = struct('distance',.5,'resolution',uint64([5,2]));

bickleyJet.shearline = set_shearline_resolution([2,1]*10);
bickleyJet.shearline = set_shearline_max_length(1e6,bickleyJet.shearline);
bickleyJet.shearline = set_shearline_average_geodesic_deviation_tol([inf,inf],bickleyJet.shearline);

end

% Forced-damped Duffing oscillator used with aperiodic forcing
function dPhi = d_phi(tau,phi)
dPhi(2,1) = nan;

dPhi(1) = phi(2);
dPhi(2) = -.1*phi(2) - phi(1)^3 + 11*cos(tau);
end
