function bi = getLValue(i,tis,p,f)
[n,~] = size(tis);
bi = 0;

for m = 1:n
    tri = tis(m,:);
    
    [at,phit,~]=getLocalBase(p,tri);
    
    tmploc = 1:3;
    
    tmpi = tmploc(tri == i);
    
    phi = phit{tmpi};
    
    product = @(p) f(p)* phi(p);
    
    a12 = (at(1,:)+at(2,:))./2;
    
    a23 = (at(2,:)+at(3,:))./2;
    
    a13 = (at(1,:)+at(3,:))./2;
    
    sumaNodos = product(a12) + product(a23) + product(a13);
    
    bi = bi + (abs(det(bsxfun(@minus,p(tri(2:end),:),p(tri(1),:))))/2)*(1/3)*sumaNodos;
end

% bi = (1/6)*f(p(i,:))*bi;

end