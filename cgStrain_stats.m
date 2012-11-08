function o = cgStrain_stats(cgStrain,cgStrainEigenvalue,verbose)

prodCgStrainD = prod(cgStrainEigenvalue,2);

o(1) = 1 - max(cgStrainEigenvalue(:,1));
detCgStrain = arrayfun(@(idx)det(cgStrain(:,:,idx)),1:size(cgStrain,3));
o(2) = mean(abs(detCgStrain-1));

if verbose
    fprintf('lambda_1*lambda2:\n')
    fprintf('min = %g\n',min(prodCgStrainD))
    fprintf('max = %g\n',max(prodCgStrainD))
    fprintf('mean = %g\n',mean(prodCgStrainD))
    fprintf('median = %g\n',median(prodCgStrainD))
    
    fprintf('\n')
    fprintf('1 - max(lambda_1) = %g\n',o(1))
    fprintf('mean(abs(detCgStrain-1)) = %g\n',o(2))
end
