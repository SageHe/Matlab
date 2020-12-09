function dydt = odefunp4(t,y0,params)
dydt = zeros(6,1);

x = y0(1);
y = y0(2);
z = y0(3);
u = y0(4);
v = y0(5);
w = y0(6);

mu = params{1};
Re = params{2};
U = params{3};

r = sqrt(x^2 + y^2 + z^2);

dydt(1) = u;
dydt(2) = v;
dydt(3) = w;
dydt(4) = -mu*x/(r^3) + U(1);
dydt(5) = -mu*y/(r^3) + U(2);
dydt(6) = -mu*z/(r^3) + U(3);
end