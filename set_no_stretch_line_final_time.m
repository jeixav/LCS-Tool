function noStretchLine = set_no_stretch_line_final_time(finalTime,...
    noStretchLine)
% Set the no stretch line initial final time and delete fields that depend
% on it from the no stretch line structure.

validateattributes(finalTime,{'double'},{'scalar','positive'})

noStretchLine.finalTime = finalTime;

fieldsToDelete = {'positionNeg','positionPos'};

for iField = 1:length(fieldsToDelete)
    if isfield(noStretchLine,fieldsToDelete{iField})
        noStretchLine = rmfield(noStretchLine,fieldsToDelete{iField});
    end
end
