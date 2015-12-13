function [index, pInterior] = essentialBoundaryOnCircleFilter(p,xr,yr,r)
[n, ~] = size(p);
allIndex = 1:n;
distanceSquare = @(p) sum(bsxfun(@minus,p,[xr yr]).^2,2);
boolean = r^2 -distanceSquare(p) < eps;
pInterior = p(~boolean);
index = allIndex(boolean);
end
