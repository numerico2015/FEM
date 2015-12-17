function [passed,p,t,u,uhtotal] = problem1(h,beta, tol)
    f = @(x) 2 - x * beta';
    
    [p,t] = circleMesh(h);
    
    [A, b] = getStiffnessMatrixAndLVector(p,t,beta,f,1);
    
    [~,boolean] = essentialBoundaryOnCircleFilter(p,t,1);
    
    uh = A\b;
    
    [n,~] = size(p);
    
    uhtotal = zeros(n,1);
    
    uhtotal(~boolean) = uh;
    
    u = (1 - sum(p.^2,2))./2;
    
    passed = norm(u-uhtotal,2)<tol;
    
    figure(1)
    plot(uhtotal,'r');
    hold on;
    plot(u,'g');
    legend('reconstruction','solution')
end