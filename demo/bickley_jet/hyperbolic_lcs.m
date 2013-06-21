%% Input parameters
u = 62.66;
earthRadius = 6.371e6;

lengthX = pi*earthRadius;
lengthY = 1.77e6;

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true,false];
perturbationCase = 3;
bickleyJet.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,lengthX,lengthY,perturbationCase),bickleyJet.flow);
bickleyJet.flow = set_flow_resolution([500,200],bickleyJet.flow);
% magicNumber gives the domain an aspect ratio similar to that used in
% doi:10.1016/j.physd.2012.06.012 and ensures grid spacing is equal in the
% x and y directions. In doi:10.1016/j.physd.2012.06.012, 
% magicNumber = 2.2599.
magicNumber = .5*pi*earthRadius/lengthY*double(bickleyJet.flow.resolution(2)-1)/double(bickleyJet.flow.resolution(1)-1);
bickleyJet.flow = set_flow_domain([0,lengthX;[-1,1]*magicNumber*lengthY],bickleyJet.flow);
bickleyJet.flow = set_flow_timespan([0,4*lengthX/u],bickleyJet.flow);

bickleyJet.strainline = set_strainline_max_length(1e8);
gridSpace = diff(bickleyJet.flow.domain(1,:))/(double(bickleyJet.flow.resolution(1))-1);
localMaxDistance = 8*gridSpace;

%% Compute λ₂ and ξ₁
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(bickleyJet.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(bickleyJet.flow.resolution),2]);

%% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue2,diff(bickleyJet.flow.timespan));
hAxes = setup_figure(bickleyJet.flow.domain);
hImagesc = imagesc(bickleyJet.flow.domain(1,:),bickleyJet.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow

%% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda2_max(localMaxDistance,bickleyJet.strainline.maxLength,cgEigenvalue2,cgEigenvector1,bickleyJet.flow.domain);

%% Plot strainlines
hStrainline = cellfun(@(position)plot(position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','w')
set(hStrainline,'lineWidth',2)
hStrainlineInitialPosition = arrayfun(@(idx)plot(strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:numel(strainlinePosition));
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','k')
set(hStrainlineInitialPosition,'MarkerFaceColor','k')
