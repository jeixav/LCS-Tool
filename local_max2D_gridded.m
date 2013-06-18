% Locally largest element in 2D array.
%
% SYNTAX
% [c,i] = local_max(array,radius)
%
% EXAMPLE
% C = peaks;
% [c,i] = local_max2D_gridded(C,uint8(5));
% imagesc(C)
% hold on
% plot(i(:,2),i(:,1),'ko')

function [c,i] = local_max2D_gridded(array,radius)

validateattributes(array,{'numeric'},{'2d'})
validateattributes(radius,{'uint8','uint16','uint32','uint64'},{'scalar','>',0})

circle = circle_mask(radius);

c = nan(size(array));
i = false(size(array));

for m = 1:size(array,1)
    for n = 1:size(array,2)
        centredMask = centred_mask([m,n],size(array),circle);
        if ~any(array(m,n) < array(centredMask))
            c(m,n) = array(m,n);
            i(m,n) = true;
        end
    end
end

c = c(~isnan(c));
[row,col] = find(i);
i = [row,col];
