function dydt = odefunL(t,y,deltau)

K = 1000;
g = 50;
G = 6.673e-11;
M = 5.98e24;
R = 6.37e6;
x3_0 = 1000;
% deltau = 10;

u0 = (G*M*x3_0*exp((-G*M*t)/(R^2*K)))/(R^2*K);
u = u0 + deltau*abs(cos(t));

A = [0 1 0;(-2*G*M/(R^3)) (-g/x3_0) ((-K*u)/x3_0^2); 0 0 0];

B = [0;(K/x3_0);-1];

dydt = A*y + B*u;