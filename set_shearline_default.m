function shearline = set_shearline_default(shearline)
% Add default values to shearline structure.

% if ~isfield(shearline,'odeSolverOptions')
%     defaultTimestep = .025;
%     % Kludge: Store fixed timestep in odeset's InitialStep parameter.
%     shearline.odeSolverOptions = odeset('InitialStep',defaultTimestep);
% end

