function cgStrain_stats(cgStrainEigenvalue)

prodCgStrainD = prod(cgStrainEigenvalue,2);

fprintf('lambda_1*lambda2:\n')
fprintf('min = %g\n',min(prodCgStrainD))
fprintf('max = %g\n',max(prodCgStrainD))
fprintf('mean = %g\n',mean(prodCgStrainD))
fprintf('median = %g\n',median(prodCgStrainD))

% fprintf('\n')
% fprintf('max lambda_1 = %g\n',max(cgStrainD(:,1)))
% detCgStrain = cellfun(@det,cgStrain);
% fprintf('mean(abs(detCgStrain-1)) = %g\n',mean(abs(detCgStrain-1)))
