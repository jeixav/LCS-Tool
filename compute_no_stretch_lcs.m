function [flow,noStretchLine] = compute_no_stretch_lcs(flow,noStretchLine)

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    verbose.progress = true;
    verbose.stats = false;
    cgStrainMethod.name = 'eov';
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        cgStrainMethod,verbose);
end

if ~all(isfield(noStretchLine,{'chiPos','chiNeg','positionPos',...
        'positionNeg'}))
    verbose = true;
    noStretchLine = compute_no_stretch_line(flow,noStretchLine,verbose);
end
