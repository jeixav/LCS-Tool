function superminIndex = superminimize2(position,value,superminDistance)

superminIndex = nan(size(position));
counterIndex = 0;

while ~all(isnan(value))
    counterIndex = counterIndex + 1;
    [~,minIndex] = min(value);
    value(minIndex) = nan;
    minInterval = position(minIndex) + .5*superminDistance*[-1 1];

    % CleanSuperminIndex is the same as superminIndex, but with NaN values
    % removed.
    cleanSuperminIndex = superminIndex(~isnan(superminIndex));
    if ~any((position(cleanSuperminIndex) > minInterval(1)) & ...
            (position(cleanSuperminIndex) <= minInterval(2)))
        superminIndex(counterIndex) = minIndex;
    end
end

superminIndex = superminIndex(~isnan(superminIndex));