function [index, boolean] = essentialBoundaryOnCircleFilter(p,t,problem)

tol = 10000*eps;

[n, ~] = size(p);
e = boundedges(p,t);
e = unique(e)';

switch problem
    case 1
          index = e; 
    case 2
          boundaryPoints = p(e,:);
          filtro = ( abs(boundaryPoints(:,1)) < tol & ~(abs(boundaryPoints(:,2)-1) < tol) ) | ( abs(boundaryPoints(:,2)) < tol & ~(abs(boundaryPoints(:,1)-1) < tol));
          index = e(~filtro);
end
b = zeros(1,n);
b(index) = 1;
boolean = logical(b);
end
