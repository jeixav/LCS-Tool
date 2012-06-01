function side = side_point_line(point,line)
%SIDE_POINT_LINE Determine on which side of a line a point lies

cellPoint = mat2cell(point,ones(size(point,1),1));

% Formula based on determinant
% Reference: http://stackoverflow.com/questions/1560492/
side = sign(cellfun(@(iPoint) det([line(2,:)-line(1,:); ...
    iPoint-line(1,:)]),cellPoint));

end
