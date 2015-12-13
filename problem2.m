function passed = problem2(h,beta, tol)
    f = @(x) 2 - x * beta';
    
    [p,t] = pacman(h);
    
    [A, b] = getStiffnessMatrixAndLVector(p,t,beta,f);
    
    [~,pInterior] = essentialBoundaryOnCircleFilter(p,0,0,1);
    
    uh = A\b;
    
    u = (1 - sum(pInterior.^2,2))./2;
    
    passed = norm(u-uh,2)<tol;
    
    figure(1)
    plot(uh,'r');
    hold on;
    plot(u,'g');
    legend('reconstruction','solution')
end