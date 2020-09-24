function [F] = fFromA(X,Y)
mu = 398600;
r = 6678;
A = [0 1 0 0;
     -mu*(X^2+Y^2)^(-3/2)+3*mu*X^2*(X^2+Y^2)^(-5/2) 0 mu*X*3*Y*(X^2+Y^2)^(-5/2) 0;
     0 0 0 1;
     mu*X*3*Y*(X^2+Y^2)^(-5/2) 0 -mu*(X^2+Y^2)^(-3/2)+mu*Y^2*3*(X^2+Y^2)^(-5/2) 0;];
% A = [0 1 0 0;
%     -mu/r^3 0 0 0;
%     0 0 0 1;
%     0 0 -mu/r^3 0];

F = expm(A*10);


end