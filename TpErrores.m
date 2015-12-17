h = [0.05,0.1,0.2];

n = length(h);

beta = [0,0];

einf1 = zeros(1, n);
einf2 = zeros(1, n);
e21 = zeros(1, n);
e22 = zeros(1, n);

for i=1:n
    %%%%%%%%%%%Problema 1 norma infinito%%%%%%%%%%%%%%%%%%%
    einf1(i) = errores(1, beta, 1, h(i));
    %%%%%%%%%%%Problema 2 norma infinito%%%%%%%%%%%%%%%%%%%
    einf2(i) = errores(2, beta, 1, h(i));
    %%%%%%%%%%%Problema 1 norma 2%%%%%%%%%%%%%%%%%%%
    e21(i) = errores(1, beta, 2, h(i));
    %%%%%%%%%%%Problema 2 norma 2%%%%%%%%%%%%%%%%%%%
    e22(i) = errores(2, beta, 2, h(i));
end

figure(1)
plot(h, einf1);
title('Error en norma infinito para el problema 1');

figure(2)
plot(h, einf2);
title('Error en norma infinito para el problema 2');

figure(3)
plot(h, e21);
title('Error en norma 2 para el problema 1');

figure(4)
plot(h, e22);
title('Error en norma 2 para el problema 2');