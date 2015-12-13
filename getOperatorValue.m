function aji = getOperatorValue(j,i,beta,tis,p)

aji = 0;

tijs = tis(tis(:,1) == j | tis(:,2) == j | tis(:,3) == j,:,:);

[n,~] = size(tijs);

for m = 1:n;
    
    triangle = tijs(m,:);
    
    tmploc = 1:3;
    
    tmpi = tmploc(triangle == i);
    
    tmpj = tmploc(triangle == j);
    
    [at,~,gradphit]=getLocalBase(p,triangle);
        
    T_area = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2; 
    
    aji = aji + dot(gradphit(tmpi,:) + beta * (1/6),gradphit(tmpj,:)) * T_area;  
    
end

end