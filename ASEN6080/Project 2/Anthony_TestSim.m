clear all
n = 7;
syms x y z xdot ydot zdot mu R AU P_Phi r Cr ri_sc xs ys zs Ri_c Am Phi_SRP uhat mu_s w_E Xs Ys Zs T tau ti timinus sigma Dt real 
sym_X = [x y z xdot ydot zdot Cr];
r = [x;y;z]; 
Rs_c = [xs;ys;zs];
r_norm = (x^2+y^2+z^2)^(1/2);
Rs_c_norm = (xs^2 + ys^2 + zs^2)^(1/2);
rs_sc = Rs_c - r;
rs_sc_norm = ((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(1/2);
uhat = rs_sc/rs_sc_norm;
Asrp = -Cr*(P_Phi/rs_sc_norm^2)*(Am)*AU^2*uhat;
A3bp = mu_s*(rs_sc/(rs_sc_norm^3) - Rs_c/(Rs_c_norm^3));
U = (mu/r_norm);
% U_J2J3 = (mu/r)*(1-J2*R^2*(3*z^2-r^2)/(2*r^4)-J3*R^3*z*(5*z^2-3*r^2)/(2*r^6));
ax = diff(U,x);
ay = diff(U,y);
az = diff(U,z);
a = [ax;ay;az] + A3bp + Asrp
dUdZ = [xdot ydot zdot a' 0];
A = jacobian(dUdZ,sym_X)