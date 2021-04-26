clearvars -except A ax ay az Asrp A3bp a

syms x y z xs ys zs xdot ydot zdot mu mu_s Am_mat Cr AU ri_sc a_srp a_3BG P_phi real

R_sc = [x y z]';
R_sc_norm = (x^2 + y^2 + z^2)^(1/2);

V_sc = [xdot ydot zdot]';

Rs_e = [xs ys zs]';
Rs_e_norm = (xs^2 + ys^2 + zs^2)^(1/2);

ri_sc = Rs_e - R_sc;
ri_sc_norm = ((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(1/2);
uhat = ri_sc/ri_sc_norm;

U = (mu/R_sc_norm);

a_kep = simplify(jacobian(U,[x,y,z]))';
% ax = diff(U,x);ay = diff(U,y);az = diff(U,z);
% a_kep = [ax ay az]';

a_srp = -Cr*(P_phi/ri_sc_norm^2)*Am_mat*AU^2*uhat;

a_3BG = mu_s*((ri_sc/ri_sc_norm^3) - (Rs_e/Rs_e_norm^3));

a_tot = a_kep + a_srp + a_3BG;

f = [xdot ydot zdot a_tot' 0];

state = [x y z xdot ydot zdot Cr];

Amat = jacobian(f,state);
