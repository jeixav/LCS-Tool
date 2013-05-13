function d = drectangle(p,x1,x2,y1,y2,periodicBc)
% Compute signed distance function for rectangle with corners (x1,y1), 
% (x2,y1), (x1,y2), (x2,y2).
%
% Copied from DistMesh, http://persson.berkeley.edu/distmesh/, (GPL license.)
% Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

d = -min(min(min(-y1+p(2),y2-p(2)),-x1+p(1)),x2-p(1));

if periodicBc(1)
    d = -min(-y1+p(2),y2-p(2));
end

if periodicBc(2)
    warning('Periodic BC in y not programmed')
end