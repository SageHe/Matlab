function zd = keplerJ2CD_wPhi_ODE(t,Z)
Re = 6378.1363; %Earth radius in km

x = Z(1);
y = Z(2);
z = Z(3);

xdot = Z(4);
ydot = Z(5);
zdot = Z(6);

mu = Z(7);
J2 = Z(8);
CD = Z(9);

xs1 = Z(10);
xy1 = Z(11);
xz1 = Z(12);

xs2 = Z(13);
ys2 = Z(14);
zs2 = Z(15);

xs3 = Z(16);
ys3 = Z(17);
zs3 = Z(18);

r = sqrt(x^2 + y^2 + z^2);
end