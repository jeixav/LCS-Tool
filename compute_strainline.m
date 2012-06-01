function strainlinePosition = compute_strainline(flow,strainline,...
    position,eigenvalue,eigenvector,scalingFlag)

strainlineIc = initialize_ic_grid(strainline.resolution,flow.domain);

tspan = 0:strainline.timestep:strainline.finalTime;

if ~isempty(odeget(strainline.odeSolverOptions,'OutputFcn')) && ...
        strcmp(func2str(odeget(strainline.odeSolverOptions,'OutputFcn')),...
        'ode_progress_bar')
    progressBar = ConsoleProgressBar;
    progressBar.setText(mfilename)
    progressBar.setTextPosition('left')
    progressBar.setElapsedTimeVisible(1);
    progressBar.setRemainedTimeVisible(1);
    progressBar.setLength(20)
    progressBar.start
else
    progressBar = false;
end

strainlinePositionFw = xi1_tracing(strainlineIc,tspan,position,...
    flow.resolution,eigenvalue,eigenvector,scalingFlag);

if isa(progressBar,'ConsoleProgressBar')
        progressBar.setValue(progressBar.maximum/2)
end

strainlinePositionBw = xi1_tracing(strainlineIc,-tspan,position,...
    flow.resolution,eigenvalue,eigenvector,scalingFlag);

if isa(progressBar,'ConsoleProgressBar')
    progressBar.setValue(progressBar.maximum)
    progressBar.stop
    fprintf('\n')
end

strainlinePosition = cellfun(...
    @(fw,bw)remove_repeated_ic(fw,bw),...
    strainlinePositionFw,strainlinePositionBw,...
    'UniformOutput',false);

is_position_outside_domain(strainlinePosition,flow.domain)

end

function position = remove_repeated_ic(positionFw,positionBw)
%REMOVE_REPEATED_IC Remove repeated initial condition.
%   Forward and backward time strainline integrations start from the same
%   initial position. This function concatenates the forward and backward
%   position arrays, but first removes the repeated position.

if ~all(positionFw(1,:) == positionBw(1,:))
    error('Initial positions of forward and backward integrations unequal')
end

position = [flipud(positionFw); positionBw(2:end,:)];

end
