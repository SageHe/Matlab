function dydt = odefun(t,y0,params)
dydt = zeros(6,1);

x = y0(1);
y = y0(2);
z = y0(3);
u = y0(4);
v = y0(5);
w = y0(6);

r = sqrt(x^2 + y^2 + z^2);

J2 = params(1);
R0 = params(2);
mu = params(3);
coeff = -(3/2)*mu*R0^2*J2/(r^7);

uj2_1 = coeff*(x^2 + y^2 - 4*z^2)*x;
uj2_2 = coeff*(x^2 + y^2 - 4*z^2)*y;
uj2_3 = coeff*(3*(x^2 + y^2) - 2*z^2)*z;

dydt(1) = u;
dydt(2) = v;
dydt(3) = w;
dydt(4) = (-mu*x/(r^3)) + uj2_1;
dydt(5) = (-mu*y/(r^3)) + uj2_2;
dydt(6) = (-mu*z/(r^3)) + uj2_3;

end