function zd = keplerJ2_wPhi_ODE(t,Z)
mu = 3.986004415e5;
x = Z(1);
y = Z(2);
z = Z(3);
xdot = Z(4);
ydot = Z(5);
zdot = Z(6);
J2 = Z(7);
Phi_flat = Z(8:end);
Phi = reshape(Phi_flat,7,7);
rdot_vec = [xdot;ydot;zdot];

r = sqrt(x^2 + y^2 + z^2);
R = 6378.0;

xdd = (mu*((J2*R^2*x)/(x^2 + y^2 + z^2)^2 - (2*J2*R^2*x*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*x*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + 1))/(x^2 + y^2 + z^2)^(3/2);
ydd = (mu*((J2*R^2*y)/(x^2 + y^2 + z^2)^2 - (2*J2*R^2*y*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*y*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + 1))/(x^2 + y^2 + z^2)^(3/2);
zdd = - (mu*((2*J2*R^2*z)/(x^2 + y^2 + z^2)^2 + (2*J2*R^2*z*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*z*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + 1))/(x^2 + y^2 + z^2)^(3/2);

vdot_vec = [xdd;ydd;zdd];

xdot_vec = [rdot_vec;vdot_vec];

Amat(1,:) = [0, 0, 0, 1, 0, 0, 0];
Amat(2,:) = [0, 0, 0, 0, 1, 0, 0];
Amat(3,:) = [0, 0, 0, 0, 0, 1, 0];
Amat(4,:) = [(mu*(12*J2*R^2*x^4 + 9*J2*R^2*x^2*y^2 - 81*J2*R^2*x^2*z^2 - 3*J2*R^2*y^4 + 9*J2*R^2*y^2*z^2 + 12*J2*R^2*z^4 + 4*x^6 + 6*x^4*y^2 + 6*x^4*z^2 - 2*y^6 - 6*y^4*z^2 - 6*y^2*z^4 - 2*z^6))/(2*(x^2 + y^2 + z^2)^(9/2)), (3*mu*x*y*(5*J2*R^2*x^2 + 5*J2*R^2*y^2 - 30*J2*R^2*z^2 + 2*x^4 + 4*x^2*y^2 + 4*x^2*z^2 + 2*y^4 + 4*y^2*z^2 + 2*z^4))/(2*(x^2 + y^2 + z^2)^(9/2)), (3*mu*x*z*(15*J2*R^2*x^2 + 15*J2*R^2*y^2 - 20*J2*R^2*z^2 + 2*x^4 + 4*x^2*y^2 + 4*x^2*z^2 + 2*y^4 + 4*y^2*z^2 + 2*z^4))/(2*(x^2 + y^2 + z^2)^(9/2)), 0, 0, 0, -(3*R^2*mu*x*(x^2 + y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(7/2))];
Amat(5,:) = [(3*mu*x*y*(5*J2*R^2*x^2 + 5*J2*R^2*y^2 - 30*J2*R^2*z^2 + 2*x^4 + 4*x^2*y^2 + 4*x^2*z^2 + 2*y^4 + 4*y^2*z^2 + 2*z^4))/(2*(x^2 + y^2 + z^2)^(9/2)), (mu*(- 3*J2*R^2*x^4 + 9*J2*R^2*x^2*y^2 + 9*J2*R^2*x^2*z^2 + 12*J2*R^2*y^4 - 81*J2*R^2*y^2*z^2 + 12*J2*R^2*z^4 - 2*x^6 - 6*x^4*z^2 + 6*x^2*y^4 - 6*x^2*z^4 + 4*y^6 + 6*y^4*z^2 - 2*z^6))/(2*(x^2 + y^2 + z^2)^(9/2)), (3*mu*y*z*(15*J2*R^2*x^2 + 15*J2*R^2*y^2 - 20*J2*R^2*z^2 + 2*x^4 + 4*x^2*y^2 + 4*x^2*z^2 + 2*y^4 + 4*y^2*z^2 + 2*z^4))/(2*(x^2 + y^2 + z^2)^(9/2)), 0, 0, 0, -(3*R^2*mu*y*(x^2 + y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(7/2))];
Amat(6,:) = [(3*mu*x*z*(15*J2*R^2*x^2 + 15*J2*R^2*y^2 - 20*J2*R^2*z^2 + 2*x^4 + 4*x^2*y^2 + 4*x^2*z^2 + 2*y^4 + 4*y^2*z^2 + 2*z^4))/(2*(x^2 + y^2 + z^2)^(9/2)), (3*mu*y*z*(15*J2*R^2*x^2 + 15*J2*R^2*y^2 - 20*J2*R^2*z^2 + 2*x^4 + 4*x^2*y^2 + 4*x^2*z^2 + 2*y^4 + 4*y^2*z^2 + 2*z^4))/(2*(x^2 + y^2 + z^2)^(9/2)), -(mu*(9*J2*R^2*x^4 + 18*J2*R^2*x^2*y^2 - 72*J2*R^2*x^2*z^2 + 9*J2*R^2*y^4 - 72*J2*R^2*y^2*z^2 + 24*J2*R^2*z^4 + 2*x^6 + 6*x^4*y^2 + 6*x^2*y^4 - 6*x^2*z^4 + 2*y^6 - 6*y^2*z^4 - 4*z^6))/(2*(x^2 + y^2 + z^2)^(9/2)), 0, 0, 0, -(3*R^2*mu*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(7/2))];
Amat(7,:) = [0, 0, 0, 0, 0, 0, 0];

Phi_dot = Amat*Phi;

Phi_dot_flat = reshape(Phi_dot,49,1);
zd = [xdot_vec;0;Phi_dot_flat];
end