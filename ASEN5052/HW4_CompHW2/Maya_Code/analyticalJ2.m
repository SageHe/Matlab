function dIC = analyticalJ2(t,IC,var)
x = IC(1);
y = IC(2);
z = IC(3);
u = IC(4);
v = IC(5);
w = IC(6);

mu = var(1);
RE = var(2);
J2 = var(3);

r = sqrt(x^2 + y^2 + z^2);

coeff = (-3/2)*mu*RE^2*J2/2/r^7;
u1 = x*(x^2 + y^2 - 4*z^2);
u2 = y*(x^2 + y^2 - 4*z^2);
u3 = z*(3*(x^2 + y^2) - 2*z^2);

dx = u;
dy = v;
dz = w;
du = -mu*x/(r^3) + u1*coeff;
dv = -mu*y/(r^3) + u2*coeff;
dw = -mu*z/(r^3) + u3*coeff;


dIC = [dx dy dz du dv dw]';

end