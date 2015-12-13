function [at,phit,gradphit]=getLocalBase(P,Triangle)

%Dado el triangulo r, arma una matriz at cuyas columnas son las coordenadas
%de los vertices de dicho triangulo

at = P(Triangle,:);

a1=at(1,:);
a2=at(2,:);
a3=at(3,:);

%Armado de las phi correspondientes a cada uno de los vertices

v1=[-(a3(2)-a2(2)) a3(1)-a2(1)];
v2=[-(a3(2)-a1(2)) a3(1)-a1(1)];
v3=[-(a2(2)-a1(2)) a2(1)-a1(1)];

phit=cell(3,1);

phit{1}=@(p) (((p-a2)*v1')/((a1-a2)*v1'));
phit{2}=@(p) (((p-a1)*v2')/((a2-a1)*v2'));
phit{3}=@(p) (((p-a1)*v3')/((a3-a1)*v3'));

%Armado de los gradientes de dichas phi

gradphit=zeros(3,2);
gradphit(1,:) =  v1./((a1-a2)*v1');   %[v1(1)/((a1-a2)*v1') v1(2)/((a1'-a2')*v1')];
gradphit(2,:) =  v2./((a2-a1)*v2');   %[v2(1)/((a2'-a1')*v2') v2(2)/((a2'-a1')*v2')];
gradphit(3,:) =  v3./((a3-a1)*v3');   %[v3(1)/((a3'-a1')*v3') v3(2)/((a3'-a1')*v3')];