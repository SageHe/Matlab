function dydt = odefunNL(t,y,deltau)
K = 1000;
g = 50;
G = 6.673e-11;
M = 5.98e24;
R = 6.37e6;
x3_0 = 1000;
% deltau = 10;

u0 = (G*M*x3_0*exp((-G*M*t)/(R^2*K)))/(R^2*K);
u = u0 + deltau*abs(cos(t));

dydt = zeros(3,1);

dydt(1) = y(2);
dydt(2) = ((K*u - g*y(2))/y(3)) - ((G*M)/(R + y(1))^2);
dydt(3) = -u;
end