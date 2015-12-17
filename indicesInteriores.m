function IndexInt =  indicesInteriores(p,t,problem)

    [circleBoundary,~] = essentialBoundaryOnCircleFilter(p,t,problem);
    
    [m,~] = size(p);

    IndexInt = zeros(m,1);
    
    tmpi = 0;
    
    for i = 1:m
        
        if ~any(circleBoundary == i)
            
            tmpi=tmpi+1;
        
            IndexInt(i) = tmpi;
            
        end
        
    end
    