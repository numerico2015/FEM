function [A, b] = getStiffnessMatrixAndLVector(p,t,beta,f) 

[circleBoundary,~] = essentialBoundaryOnCircleFilter(p,0,0,1);

[m,~] = size(p);

[~,l] = size(circleBoundary);

n = m - l;

A = zeros(n,n);

b = zeros(n,1);

tmpi = 0;

for i = 1:m
    if ~any(circleBoundary == i)
        tmpi=tmpi+1;
        tis = t(t(:,1) == i | t(:,2) == i | t(:,3) == i,:,:);
    
        b(tmpi) = getLValue(i,tis,p,f);
        
        tmpj = 0;
        
        for j = 1:m
            if ~any(circleBoundary == j)
                tmpj=tmpj+1;
                A(tmpi,tmpj) = getOperatorValue(j,i,beta,tis,p);
            end
        end
    end
end

end