function e = errores(problem,beta,norma,h)

     % norma = 1 es norma infinito
     % norma = 2
     
      f = @(x) 2 - x*beta';

    if norma == 2
    
        e = norm2Error(h,beta,problem);
    
    end
   
    if norma == 1
        
        if problem == 1

            [Puntos,Triangulos] = circleMesh(h);

        end

        if problem == 2
    
            [Puntos,Triangulos] = pacman(h);
    
        end
        
        [~,boolean] = essentialBoundaryOnCircleFilter(Puntos,Triangulos,problem);
        
        [A, b] = getStiffnessMatrixAndLVector(Puntos,Triangulos,beta,f,problem);
        
        uh = A\b;
    
        [n,~] = size(Puntos);
    
        uhtotal = zeros(n,1);
    
        uhtotal(~boolean) = uh;

        u = (1-sum(Puntos.^2,2))./2;
        
        e = max(abs(uhtotal-u));
        
    end
end