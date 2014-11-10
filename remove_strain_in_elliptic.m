% remove_strain_in_elliptic Remove strainlines inside elliptic LCSs
function strainlinePosition = remove_strain_in_elliptic(strainlinePosition,ellipticLcs)

strainlinePositionNew = cell(size(strainlinePosition));
strainlinePositionStart = cell(size(strainlinePosition));
strainlinePositionEnd = cell(size(strainlinePosition));

for idx = 1:numel(strainlinePosition)
    in = inpolygon(strainlinePosition{idx}(:,1),strainlinePosition{idx}(:,2),ellipticLcs(:,1),ellipticLcs(:,2));
    if any(in)
        % FIXME Algorithm below assumes a strainline segment crosses an
        % elliptic region once at most. If this is not the case, strainline
        % segments _outside_ elliptic LCSs will be discarded.
        inStart = find(in,1,'first');
        if ~isempty(inStart)
            strainlinePositionStart{idx} = strainlinePosition{idx}(1:inStart-1,:);
        end
        inEnd = find(in,1,'last');
        if ~isempty(inEnd)
            strainlinePositionEnd{idx} = strainlinePosition{idx}(inEnd+1:end,:);
        end
    else
        strainlinePositionNew{idx} = strainlinePosition{idx};
    end
end

idx = cellfun(@(input)isempty(input),strainlinePositionNew);
strainlinePosition = strainlinePosition(~idx);
idx = cellfun(@(input)isempty(input),strainlinePositionStart);
strainlinePosition = [strainlinePosition,strainlinePositionStart(~idx)];
idx = cellfun(@(input)isempty(input),strainlinePositionEnd);
strainlinePosition = [strainlinePosition,strainlinePositionEnd(~idx)];
