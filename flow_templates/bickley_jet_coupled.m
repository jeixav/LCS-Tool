% bickley_jet Defines a bickley jet.
%
% DESCRIPTION
% bickley_jet(p)
% Time periodic forcing case from doi:10.1016/j.physd.2012.06.012.
% 
% REFERENCES
% doi:10.1016/j.physd.2012.06.012
% doi:10.1063/1.3271342

function bickleyJet = bickley_jet_coupled(p)

narginchk(0,1)

defaultP.u = 62.66;
defaultP.c2 = .205*defaultP.u;
defaultP.c3 = .461*defaultP.u;
defaultP.a = 6.371e6;
defaultP.lengthX = pi*defaultP.a;
defaultP.lengthY = 1.77e6;

defaultP.epsilon1 = .075;
defaultP.epsilon2 = .4;
defaultP.epsilon3 = .3;

if nargin < 1
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

p.epsilon1 = p.epsilon1/10;
p.epsilon2 = p.epsilon2/10;

% Time-periodic psi1, case 1 on page 1691 of
% doi:10.1016/j.physd.2012.06.012.
f1 = @(t)exp(1i*sigma1*t);
f2 = @(t)exp(1i*sigma2*t);

    function derivative_ = derivative(t,x)
        
        t = timescale*t;
        
        idx1 = 1:6:size(x,1)-5;
        idx2 = 2:6:size(x,1)-4;
        
        derivative_ = nan(size(x));
        
        % u
        derivative_(idx1) = -p.c3 - p.u*(tanh(x(idx2)).^2 - 1) ...
            + 2*p.u*sech(x(idx2)).^2.*tanh(x(idx2)) ...
            .*(real(p.epsilon1*f1(t)*exp(k(1)*lengthscaleX*x(idx1)*1i) ...
            + p.epsilon2*f2(t)*exp(k(2)*lengthscaleX*x(idx1)*1i)) ...
            + p.epsilon3*cos(k(3)*lengthscaleX*x(idx1)));
        derivative_(idx1) = 2*derivative_(idx1)/lengthscaleX*timescale;
        
        % v
        derivative_(idx2) = -p.u*p.lengthY*sech(x(idx2)).^2 ...
            .*(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*lengthscaleX*x(idx1)*1i) ...
            + p.epsilon2*f2(t)*k(2)*exp(k(2)*lengthscaleX*x(idx1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*lengthscaleX*x(idx1)));
        derivative_(idx2) = derivative_(idx2)/lengthscaleY*timescale;
        
        idx3 = 3:6:size(x,1)-3;
        idx4 = 4:6:size(x,1)-2;
        idx5 = 5:6:size(x,1)-1;
        idx6 = 6:6:size(x,1);
        
        % dux
        dux = -2*p.u*sech(x(idx2)).^2.*tanh(x(idx2)) ...
            .*(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*lengthscaleX*x(idx1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*lengthscaleX*x(idx1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*lengthscaleX*x(idx1)));
        dux = 2*dux*timescale;
        
        % duy
        duy = -2*p.u/p.lengthY*sech(x(idx2)).^2 ...
            .*((3*tanh(x(idx2)).^2 - 1) ...
            .*(real(p.epsilon1*f1(t)*exp(k(1)*lengthscaleX*x(idx1)*1i)) ...
            + real(p.epsilon2*f2(t)*exp(k(2)*lengthscaleX*x(idx1)*1i)) ...
            + p.epsilon3*cos(k(3)*lengthscaleX*x(idx1))) + tanh(x(idx2)));
        duy = 2*duy/lengthscaleX*lengthscaleY*timescale;
        
        % dvx
        dvx = -p.u*p.lengthY*sech(x(idx2)).^2 ...
            .*(real(p.epsilon1*f1(t)*k(1)^2*exp(k(1)*lengthscaleX*x(idx1)*1i)) ...
            + real(p.epsilon2*f2(t)*k(2)^2*exp(k(2)*lengthscaleX*x(idx1)*1i)) ...
            + p.epsilon3*k(3)^2*cos(k(3)*lengthscaleX*x(idx1)));
        dvx = dvx/lengthscaleY*lengthscaleX*timescale;
        
        % dvy
        dvy = 2*p.u*sech(x(idx2)).^2.*tanh(x(idx2)) ...
            .*(imag(p.epsilon1*f1(t)*k(1)*exp(k(1)*lengthscaleX*x(idx1)*1i)) ...
            + imag(p.epsilon2*f2(t)*k(2)*exp(k(2)*lengthscaleX*x(idx1)*1i)) ...
            + p.epsilon3*k(3)*sin(k(3)*lengthscaleX*x(idx1)));
        dvy = dvy*timescale;

        % Perform matrix multiplication manually
        derivative_(idx3) = dux.*x(idx3) + duy.*x(idx5);
        derivative_(idx4) = dux.*x(idx4) + duy.*x(idx6);
        derivative_(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
        derivative_(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
        
    end

bickleyJet.flow.derivative = @(t,x)derivative(t,x);

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true false];
bickleyJet.flow = set_flow_domain([0 2*pi; [-1 1]*2],bickleyJet.flow);
bickleyJet.flow = set_flow_timespan([0 4],bickleyJet.flow);
bickleyJet.flow = set_flow_resolution(200,bickleyJet.flow);

bickleyJet.strainline = set_strainline_resolution([2 1]*5);
bickleyJet.strainline = set_strainline_max_length(20,bickleyJet.strainline);
bickleyJet.strainline = set_strainline_geodesic_deviation_tol(1e-5,...
    bickleyJet.strainline);
bickleyJet.strainline = set_strainline_length_tol(4,bickleyJet.strainline);
bickleyJet.strainline.filteringMethod = 'superminimize';
bickleyJet.strainline.filteringParameters = struct('distance',.5,...
    'resolution',uint64([5 2]));

bickleyJet.shearline = set_shearline_resolution([2 1]*10);
bickleyJet.shearline = set_shearline_max_length(20,bickleyJet.shearline);
bickleyJet.shearline = set_shearline_average_geodesic_deviation_tol(...
    [inf inf],bickleyJet.shearline);

end
