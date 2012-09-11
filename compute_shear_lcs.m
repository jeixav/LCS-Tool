function [flow,shearline] = compute_shear_lcs(flow,shearline)

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))

    % FIXME Code block copy-pasted from compute_strain_lcs -- make it a
    % shared function
    if isfield(flow,'symDerivative') && ~isfield(flow,'dDerivative')
        symJacDy = symJacDerivative(flow.symDerivative);
        
        jacDyScalar11 = matlabFunction(symJacDy{1,1},'vars',{'t','x','y'});
        jacDyScalar12 = matlabFunction(symJacDy{1,2},'vars',{'t','x','y'});
        jacDyScalar21 = matlabFunction(symJacDy{2,1},'vars',{'t','x','y'});
        jacDyScalar22 = matlabFunction(symJacDy{2,2},'vars',{'t','x','y'});
        
        flow.dDerivative = @(t,y)[jacDyScalar11(t,y(1),y(2)) ...
            jacDyScalar12(t,y(1),y(2)); jacDyScalar21(t,y(1),y(2)) ...
            jacDyScalar22(t,y(1),y(2))];
    end

    verbose.progress = true;
    verbose.stats = false;
%     eig_cgStrainMethod.name = 'fd';
%     eig_cgStrainMethod.params = struct('eigenvalueFromMainGrid',true);
    eig_cgStrainMethod.name = 'eov';
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        eig_cgStrainMethod,verbose);
end

if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    verbose = true;
    shearline = compute_shearline(flow,shearline,verbose);
end

if ~all(isfield(shearline,{'geodesicDeviationPos','geodesicDeviationNeg',...
        'averageGeodesicDeviationPos','averageGeodesicDeviationNeg'}))
    shearline = geodesic_deviation_shearline(flow,shearline);
end

if ~all(isfield(shearline,{'filteredIndexPos','filteredIndexNeg'}))
    shearline.filteredIndexPos = filter_shearline(...
        shearline.averageGeodesicDeviationPos,...
        shearline.averageGeodesicDeviationPosTol);
    shearline.filteredIndexNeg = filter_shearline(...
        shearline.averageGeodesicDeviationNeg,...
        shearline.averageGeodesicDeviationNegTol);
end

function index = filter_shearline(averageGeodesicDeviation,tol)

% FIXME Workaround while boundary element geodesic is not calculated and
% yields NaNs
if isinf(tol)
    warning('averageGeodesicDeviation = inf workaround used')
    index = true(size(averageGeodesicDeviation));
else
    index = averageGeodesicDeviation <= tol;
end
