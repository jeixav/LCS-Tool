function o = cgStrain_stats(cgStrain,cgStrainEigenvalue,verbose)

prodCgStrainD = prod(cgStrainEigenvalue,2);

o(1) = 1 - max(cgStrainEigenvalue(:,1));
detCgStrain = arrayfun(@(idx)det(cgStrain(:,:,idx)),1:size(cgStrain,3));
o(2) = mean(abs(detCgStrain-1));

if verbose
    fprintf('lambda_1*lambda2:\n')
    fprintf('\tmin = %g\n',min(prodCgStrainD))
    fprintf('\tmax = %g\n',max(prodCgStrainD))
    fprintf('\tmean = %g\n',mean(prodCgStrainD))
    fprintf('\tmedian = %g\n',median(prodCgStrainD))
    fprintf('1 - max(lambda_1) = %g\n',o(1))
    fprintf('mean(abs(detCgStrain-1)) = %g\n',o(2))
end
