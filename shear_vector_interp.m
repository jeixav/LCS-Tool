function shearVectorInterp = shear_vector_interp(position,cgPosition,...
    cgResolution,cgEigenvector,cgEigenvalue)

meshgridCgPositionX = reshape(cgPosition(:,1),fliplr(cgResolution));
meshgridCgPositionY = reshape(cgPosition(:,2),fliplr(cgResolution));

cgEigenvectorInterp(:,1) = interp2(meshgridCgPositionX,meshgridCgPositionY,...
    reshape(cgEigenvector(:,1),fliplr(cgResolution)),...
    position(:,1),position(:,2));
cgEigenvectorInterp(:,1) = interp2(meshgridCgPositionX,meshgridCgPositionY,...
    reshape(cgEigenvector(:,1),fliplr(cgResolution)),...
    position(:,1),position(:,2));

% cgEigenvectorInterp = interp2(cgPosition,cgEigenvector,position);
cgEigenvalueInterp = interp2(cgPosition,cgEigenvalue,position);

shearVectorInterp = shear_vector(cgEigenvectorInterp,cgEigenvalueInterp,...
    plusOrMinus);

end

function shearVector = shear_vector(cgEigenvector,cgEigenvalue,plusOrMinus)
%SHEAR_VECTOR Normalized Lagrangian shear vector.
%   Equation 18 from Physica D submitted article.

denominator = sqrt(cgEigenvalue(1)) + sqrt(cgEigenvalue(2));

term1 = sqrt(sqrt(lambda2)/denominator)*cgEigenvector(:,1);
term2 = sqrt(sqrt(lambda1)/denominator)*cgEigenvector(:,2);

switch plusOrMinus
    case 'plus'
        shearVector = term1 + term2;
    case 'minus'
        shearVector = term1 - term2;
    otherwise
        error('Invalid value for plusOrMinus')
end

end
