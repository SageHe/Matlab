function [dydt] = odeFunc_States(t,y)

mu = 398600;
r = 6678;
Re = 6378;
we = 2*pi/86400;
X = y(1);
Xdot = y(2);
Y = y(3);
Ydot = y(4);

r = sqrt(X^2+Y^2);

dydt(1) = y(2);
dydt(2) = -(mu*X)/r^3;
dydt(3) = y(4);
dydt(4) = -(mu*Y)/r^3;

dydt = dydt';
end