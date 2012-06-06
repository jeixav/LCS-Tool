function noStretchLine = set_no_stretch_line_final_time(finalTime,...
    noStretchLine)
% Set the no stretch line initial final time and delete fields that depend
% on it from the no stretch line structure.

validateattributes(finalTime,{'double'},{'scalar','positive'})

noStretchLine.finalTime = finalTime;

if isfield(noStretchLine,'positionPos')
    noStretchLine = rmfield(noStretchLine,'positionPos');
end

if isfield(noStretchLine,'positionNeg')
    noStretchLine = rmfield(noStretchLine,'positionNeg');
end
