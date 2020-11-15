function dIC = analyticalg(t,IC,var)
x = IC(1);
y = IC(2);
z = IC(3);
u = IC(4);
v = IC(5);
w = IC(6);

mu = var{1};
RE = var{2};
U = var{3};

r = sqrt(x^2 + y^2 + z^2);


dx = u;
dy = v;
dz = w;
du = -mu*x/(r^3) + U(1);
dv = -mu*y/(r^3) + U(2);
dw = -mu*z/(r^3) + U(3);


dIC = [dx dy dz du dv dw]';

end