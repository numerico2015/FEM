function e = norm2Error(h,beta,problem)

f = @(x) 2 - x*beta';
if problem == 1

    [Puntos,Triangulos] = circleMesh(h);

end

if problem == 2
    
    [Puntos,Triangulos] = pacman(h);
    
end
[A, b] = getStiffnessMatrixAndLVector(Puntos,Triangulos,beta,f,problem);

uh = A\b;

u =@(p) (1-p(1)^2 - p(2)^2)/2;

[cantidadDeTriangulos,~] = size(Triangulos);

IndexInt =  indicesInteriores(Puntos,Triangulos,problem);

e = 0;

for t=1:cantidadDeTriangulos
    
   %indices de los puntos del triangulo
    
   triangulo = Triangulos(t,:);
   
   indice_1 = triangulo(1);
    
   indice_2 = triangulo(2);
   
   indice_3 = triangulo(3);
   
   
    
    
    %---------- Triangulo con 2 puntos en la frontera ------------%
    
    if ((0 == IndexInt(indice_1)) && (0 == IndexInt(indice_2))) && ~(0 == IndexInt(indice_3))
        
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        base = phit{3};
        
        phi = @(p) uh(IndexInt(indice_3))*base(p);
        
        %---puntos intermedios----%
        
        p13 = (at(1,:) + at(2,:))/2;
        
        p23 = (at(3,:) + at(2,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p13)-phi(p13))^2 + (u(p23)-phi(p23))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
        
        
    end
    
    if (0 == IndexInt(indice_1)) && (0 == IndexInt(indice_3)) && ~(0 == IndexInt(indice_2))
    
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        phi = @(p) uh(IndexInt(indice_2))*phit{2}(p);
        
        %---puntos intermedios----%
        
        p12 = (at(1,:) + at(2,:))/2;
        
        p23 = (at(3,:) + at(2,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p12)-phi(p12))^2 + (u(p23)-phi(p23))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
        
    end
    
    if ((0 == IndexInt(indice_2)) && (0 == IndexInt(indice_3))) && ~(0 == IndexInt(indice_1))
    
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        phi = @(p) uh(IndexInt(indice_1))*phit{1}(p);
        
        %---puntos intermedios----%
        
        p12 = (at(1,:) + at(2,:))/2;
        
        p13 = (at(3,:) + at(1,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p12)-phi(p12))^2 + (u(p13)-phi(p13))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
    end
    
     %---------- Triangulo con 1 punto en la frontera ------------%
    
     if (0 == IndexInt(indice_1)) && ~(0 == IndexInt(indice_2)) && ~(0 == IndexInt(indice_3))
    
         
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        phi = @(p) uh(IndexInt(indice_2))*phit{2}(p) + uh(IndexInt(indice_3))*phit{3}(p);
        
        %---puntos intermedios----%
        
        p12 = (at(1,:) + at(2,:))/2;
        
        p23 = (at(3,:) + at(2,:))/2;
        
        p13 = (at(3,:) + at(1,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p12)-phi(p12))^2 +(u(p23)-phi(p23))^2+ (u(p13)-phi(p13))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
         
         
    end
    
    if ~(0 == IndexInt(indice_1)) && (0 == IndexInt(indice_2)) && ~(0 == IndexInt(indice_3))
    
         
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        phi = @(p) uh(IndexInt(indice_1))*phit{1}(p) + uh(IndexInt(indice_3))*phit{3}(p);
        
        %---puntos intermedios----%
        
        p12 = (at(1,:) + at(2,:))/2;
        
        p23 = (at(3,:) + at(2,:))/2;
        
        p13 = (at(3,:) + at(1,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p12)-phi(p12))^2 +(u(p23)-phi(p23))^2+ (u(p13)-phi(p13))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
    
    end
    
    if ~(0 == IndexInt(indice_1)) && ~(0 == IndexInt(indice_2)) && (0 == IndexInt(indice_3))
    
         
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        phi = @(p) uh(IndexInt(indice_1))*phit{1}(p) + uh(IndexInt(indice_2))*phit{2}(p);
        
        %---puntos intermedios----%
        
        p12 = (at(1,:) + at(2,:))/2;
        
        p23 = (at(3,:) + at(2,:))/2;
        
        p13 = (at(3,:) + at(1,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p12)-phi(p12))^2 +(u(p23)-phi(p23))^2+ (u(p13)-phi(p13))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
    
    
    end
    
    
    %---------- Triangulo en el interior ------------%
    
    
    if ~(0 == IndexInt(indice_1)) && ~(0 == IndexInt(indice_2)) && ~(0 == IndexInt(indice_3))
    
         
        [at,phit,~]=getLocalBase(Puntos,triangulo);
        
        phi = @(p) uh(IndexInt(indice_1))*phit{1}(p) + uh(IndexInt(indice_2))*phit{2}(p) + uh(IndexInt(indice_3))*phit{3}(p);
        
        %---puntos intermedios----%
        
        p12 = (at(1,:) + at(2,:))/2;
        
        p23 = (at(3,:) + at(2,:))/2;
        
        p13 = (at(3,:) + at(1,:))/2;
        
        %-- Area del trinagulo----%
        
        Area_T = abs(det(bsxfun(@minus,at((2:end),:),at(1,:))))/2;
        
        %--- Integro usando cuadratura---%
        
        error_T =  ((u(p12)-phi(p12))^2 +(u(p23)-phi(p23))^2+ (u(p13)-phi(p13))^2)*Area_T/3;
        
        % sumo al error
        
        e = e +error_T;
    end
    
end
    
e = sqrt(e);

            
            
            
            
            
            
            
            
        
        
        
        
    

end