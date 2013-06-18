% Enlarge a circular mask array to an array of size sizeArray centered at
% centre. This involves "padding" circleArray with zeros.
function centredMask = centred_mask(centre,sizeArray,circleMask)

% validateattributes(sizeArray,uint,{'size',[1,2],'>=',1},mfilename,'sizeArray',2)
% validateattributes(centre,uint,{'size',[1,2],'>=',1},mfilename,'centre',1)
% validateattributes(centre(1),uint,{'scalar','<=',sizeArray(1)},mfilename,'centre(1)',1)
% validateattributes(centre(2),uint,{'scalar','<=',sizeArray(2)},mfilename,'centre(2)',1)
% validateattributes(circleMask,{'logical'},{'square'})

if ~rem(size(circleMask,1),2)
    error('circleMask must have odd dimensions')
else
    radius = (size(circleMask,1) - 1)/2 + 1;
end

% Enlarged array notation
% ┌─┬─┬─┐
% │ │1│ │
% │ ├─┤ │
% │4│K│2│
% │ ├─┤ │
% │ │3│ │
% └─┴─┴─┘

% Enlarge array with false entries
array = true(sizeArray);

enlarge1nRows = radius - centre(1);
if enlarge1nRows < 0
    enlarge1nRows = 0;
end
enlarge1nCols = sizeArray(2);
enlarge1 = false(enlarge1nRows,enlarge1nCols);

enlarge3nRows = centre(1) + radius - 1 - sizeArray(1);
if enlarge3nRows < 0
    enlarge3nRows = 0;
end
enlarge3nCols = sizeArray(2);
enlarge3 = false(enlarge3nRows,enlarge3nCols);

enlarge2nRows = sizeArray(1) + enlarge1nRows + enlarge3nRows;
enlarge2nCols = centre(2) + radius - 1 - sizeArray(2);
enlarge2 = false(enlarge2nRows,enlarge2nCols);

enlarge4nRows = sizeArray(1) + enlarge1nRows + enlarge3nRows;
enlarge4nCols = radius - centre(2);
enlarge4 = false(enlarge4nRows,enlarge4nCols);

enlargedArray = [enlarge4,[enlarge1;array;enlarge3],enlarge2];

% Enlarge circle mask with false entries
enlarge1nRows = centre(1) - radius;
if enlarge1nRows < 0
    enlarge1nRows = 0;
end
enlarge1nCols = size(circleMask,1);
enlarge1 = false(enlarge1nRows,enlarge1nCols);

enlarge3nRows = sizeArray(1) - (radius - 1) - centre(1);
if enlarge3nRows < 0
    enlarge3nRows = 0;
end
enlarge3nCols = size(circleMask,1);
enlarge3 = false(enlarge3nRows,enlarge3nCols);

enlarge2nRows = size(circleMask,1) + enlarge1nRows + enlarge3nRows;
enlarge2nCols = sizeArray(2) - radius - centre(2) + 1;
enlarge2 = false(enlarge2nRows,enlarge2nCols);

enlarge4nRows = size(circleMask,1) + enlarge1nRows + enlarge3nRows;
enlarge4nCols = centre(2) - radius;
enlarge4 = false(enlarge4nRows,enlarge4nCols);

enlargedCircle = [enlarge4,[enlarge1;circleMask;enlarge3],enlarge2];

centredMask = enlargedArray & enlargedCircle;
corner1 = [radius-centre(1)+1,radius-centre(2)+1];
corner1(corner1 < 1) = 1;
corner2 = corner1 + sizeArray - [1,1];
centredMask = centredMask(corner1(1):corner2(1),corner1(2):corner2(2));
