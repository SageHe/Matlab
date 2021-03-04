function zd = keplerJ2J3_wPhi_ODE(t,Z)
mu = 3.986004415e5;
x = Z(1);
y = Z(2);
z = Z(3);
xdot = Z(4);
ydot = Z(5);
zdot = Z(6);
J2 = Z(7);
J3 = Z(8);
Phi_flat = Z(9:end);
Phi = reshape(Phi_flat,8,8);
rdot_vec = [xdot;ydot;zdot];

r = sqrt(x^2 + y^2 + z^2);
R = 6378.0;

xdd = (mu*((J2*R^2*x)/(x^2 + y^2 + z^2)^2 - (2*J2*R^2*x*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3 + (3*J3*R^3*x*z)/(x^2 + y^2 + z^2)^3 - (3*J3*R^3*x*z*(3*x^2 + 3*y^2 - 2*z^2))/(x^2 + y^2 + z^2)^4))/(x^2 + y^2 + z^2)^(1/2) - (mu*x*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + (J3*R^3*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^3) + 1))/(x^2 + y^2 + z^2)^(3/2);
ydd = (mu*((J2*R^2*y)/(x^2 + y^2 + z^2)^2 - (2*J2*R^2*y*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3 + (3*J3*R^3*y*z)/(x^2 + y^2 + z^2)^3 - (3*J3*R^3*y*z*(3*x^2 + 3*y^2 - 2*z^2))/(x^2 + y^2 + z^2)^4))/(x^2 + y^2 + z^2)^(1/2) - (mu*y*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + (J3*R^3*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^3) + 1))/(x^2 + y^2 + z^2)^(3/2);
zdd = - (mu*((2*J2*R^2*z)/(x^2 + y^2 + z^2)^2 + (2*J3*R^3*z^2)/(x^2 + y^2 + z^2)^3 - (J3*R^3*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^3) + (2*J2*R^2*z*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3 + (3*J3*R^3*z^2*(3*x^2 + 3*y^2 - 2*z^2))/(x^2 + y^2 + z^2)^4))/(x^2 + y^2 + z^2)^(1/2) - (mu*z*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + (J3*R^3*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^3) + 1))/(x^2 + y^2 + z^2)^(3/2);

vdot_vec = [xdd;ydd;zdd];

xdot_vec = [rdot_vec;vdot_vec];

Amat(1,:) = [0, 0, 0, 1, 0, 0, 0, 0];
Amat(2,:) = [0, 0, 0, 0, 1, 0, 0, 0];
Amat(3,:) = [0, 0, 0, 0, 0, 1, 0, 0];
Amat(4,:) = [(6*mu*x^2*(x^2 + y^2 + z^2)^3 - 2*mu*(x^2 + y^2 + z^2)^4 + 2*J2*R^2*mu*(x^2 + y^2 + z^2)^3 + 189*J3*R^3*mu*x^4*z - 126*J3*R^3*mu*x^2*z^3 + 35*J2*R^2*mu*x^4*(x^2 + y^2 + z^2) + 6*J3*R^3*mu*z*(x^2 + y^2 + z^2)^2 + 14*J3*R^3*mu*z^3*(x^2 + y^2 + z^2) - 25*J2*R^2*mu*x^2*(x^2 + y^2 + z^2)^2 - 5*J2*R^2*mu*y^2*(x^2 + y^2 + z^2)^2 + 10*J2*R^2*mu*z^2*(x^2 + y^2 + z^2)^2 - 105*J3*R^3*mu*x^2*z*(x^2 + y^2 + z^2) - 21*J3*R^3*mu*y^2*z*(x^2 + y^2 + z^2) + 35*J2*R^2*mu*x^2*y^2*(x^2 + y^2 + z^2) - 70*J2*R^2*mu*x^2*z^2*(x^2 + y^2 + z^2) + 189*J3*R^3*mu*x^2*y^2*z)/(2*(x^2 + y^2 + z^2)^(11/2)), (3*mu*x*y*(35*J3*R^3*x^2*z + 35*J3*R^3*y^2*z - 70*J3*R^3*z^3 + 5*J2*R^2*x^4 + 10*J2*R^2*x^2*y^2 - 25*J2*R^2*x^2*z^2 + 5*J2*R^2*y^4 - 25*J2*R^2*y^2*z^2 - 30*J2*R^2*z^4 + 2*x^6 + 6*x^4*y^2 + 6*x^4*z^2 + 6*x^2*y^4 + 12*x^2*y^2*z^2 + 6*x^2*z^4 + 2*y^6 + 6*y^4*z^2 + 6*y^2*z^4 + 2*z^6))/(2*(x^2 + y^2 + z^2)^(11/2)), (6*mu*x*z*(x^2 + y^2 + z^2)^3 - 126*J3*R^3*mu*x*z^4 + 189*J3*R^3*mu*x^3*z^2 + 6*J3*R^3*mu*x*(x^2 + y^2 + z^2)^2 - 21*J3*R^3*mu*x^3*(x^2 + y^2 + z^2) - 21*J3*R^3*mu*x*y^2*(x^2 + y^2 + z^2) + 10*J2*R^2*mu*x*z*(x^2 + y^2 + z^2)^2 - 70*J2*R^2*mu*x*z^3*(x^2 + y^2 + z^2) + 35*J2*R^2*mu*x^3*z*(x^2 + y^2 + z^2) + 189*J3*R^3*mu*x*y^2*z^2 + 35*J2*R^2*mu*x*y^2*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2)), 0, 0, 0, -(3*R^2*mu*x*(x^2 + y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)), -(5*R^3*mu*x*z*(3*x^2 + 3*y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(9/2))];
Amat(5,:) = [(3*mu*x*y*(35*J3*R^3*x^2*z + 35*J3*R^3*y^2*z - 70*J3*R^3*z^3 + 5*J2*R^2*x^4 + 10*J2*R^2*x^2*y^2 - 25*J2*R^2*x^2*z^2 + 5*J2*R^2*y^4 - 25*J2*R^2*y^2*z^2 - 30*J2*R^2*z^4 + 2*x^6 + 6*x^4*y^2 + 6*x^4*z^2 + 6*x^2*y^4 + 12*x^2*y^2*z^2 + 6*x^2*z^4 + 2*y^6 + 6*y^4*z^2 + 6*y^2*z^4 + 2*z^6))/(2*(x^2 + y^2 + z^2)^(11/2)), (6*mu*y^2*(x^2 + y^2 + z^2)^3 - 2*mu*(x^2 + y^2 + z^2)^4 + 2*J2*R^2*mu*(x^2 + y^2 + z^2)^3 + 189*J3*R^3*mu*y^4*z - 126*J3*R^3*mu*y^2*z^3 + 35*J2*R^2*mu*y^4*(x^2 + y^2 + z^2) + 6*J3*R^3*mu*z*(x^2 + y^2 + z^2)^2 + 14*J3*R^3*mu*z^3*(x^2 + y^2 + z^2) - 5*J2*R^2*mu*x^2*(x^2 + y^2 + z^2)^2 - 25*J2*R^2*mu*y^2*(x^2 + y^2 + z^2)^2 + 10*J2*R^2*mu*z^2*(x^2 + y^2 + z^2)^2 - 21*J3*R^3*mu*x^2*z*(x^2 + y^2 + z^2) - 105*J3*R^3*mu*y^2*z*(x^2 + y^2 + z^2) + 35*J2*R^2*mu*x^2*y^2*(x^2 + y^2 + z^2) - 70*J2*R^2*mu*y^2*z^2*(x^2 + y^2 + z^2) + 189*J3*R^3*mu*x^2*y^2*z)/(2*(x^2 + y^2 + z^2)^(11/2)), (6*mu*y*z*(x^2 + y^2 + z^2)^3 - 126*J3*R^3*mu*y*z^4 + 189*J3*R^3*mu*y^3*z^2 + 6*J3*R^3*mu*y*(x^2 + y^2 + z^2)^2 - 21*J3*R^3*mu*y^3*(x^2 + y^2 + z^2) - 21*J3*R^3*mu*x^2*y*(x^2 + y^2 + z^2) + 10*J2*R^2*mu*y*z*(x^2 + y^2 + z^2)^2 - 70*J2*R^2*mu*y*z^3*(x^2 + y^2 + z^2) + 35*J2*R^2*mu*y^3*z*(x^2 + y^2 + z^2) + 189*J3*R^3*mu*x^2*y*z^2 + 35*J2*R^2*mu*x^2*y*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2)), 0, 0, 0, -(3*R^2*mu*y*(x^2 + y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)), -(5*R^3*mu*y*z*(3*x^2 + 3*y^2 - 4*z^2))/(2*(x^2 + y^2 + z^2)^(9/2))];
Amat(6,:) = [(6*mu*x*z*(x^2 + y^2 + z^2)^3 - 126*J3*R^3*mu*x*z^4 + 189*J3*R^3*mu*x^3*z^2 + 6*J3*R^3*mu*x*(x^2 + y^2 + z^2)^2 - 21*J3*R^3*mu*x^3*(x^2 + y^2 + z^2) - 21*J3*R^3*mu*x*y^2*(x^2 + y^2 + z^2) + 10*J2*R^2*mu*x*z*(x^2 + y^2 + z^2)^2 - 70*J2*R^2*mu*x*z^3*(x^2 + y^2 + z^2) + 35*J2*R^2*mu*x^3*z*(x^2 + y^2 + z^2) + 189*J3*R^3*mu*x*y^2*z^2 + 35*J2*R^2*mu*x*y^2*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2)), (6*mu*y*z*(x^2 + y^2 + z^2)^3 - 126*J3*R^3*mu*y*z^4 + 189*J3*R^3*mu*y^3*z^2 + 6*J3*R^3*mu*y*(x^2 + y^2 + z^2)^2 - 21*J3*R^3*mu*y^3*(x^2 + y^2 + z^2) - 21*J3*R^3*mu*x^2*y*(x^2 + y^2 + z^2) + 10*J2*R^2*mu*y*z*(x^2 + y^2 + z^2)^2 - 70*J2*R^2*mu*y*z^3*(x^2 + y^2 + z^2) + 35*J2*R^2*mu*y^3*z*(x^2 + y^2 + z^2) + 189*J3*R^3*mu*x^2*y*z^2 + 35*J2*R^2*mu*x^2*y*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2)), -(2*mu*(x^2 + y^2 + z^2)^4 - 6*mu*z^2*(x^2 + y^2 + z^2)^3 + 4*J2*R^2*mu*(x^2 + y^2 + z^2)^3 + 126*J3*R^3*mu*z^5 - 189*J3*R^3*mu*x^2*z^3 - 189*J3*R^3*mu*y^2*z^3 + 70*J2*R^2*mu*z^4*(x^2 + y^2 + z^2) + 12*J3*R^3*mu*z*(x^2 + y^2 + z^2)^2 - 98*J3*R^3*mu*z^3*(x^2 + y^2 + z^2) + 5*J2*R^2*mu*x^2*(x^2 + y^2 + z^2)^2 + 5*J2*R^2*mu*y^2*(x^2 + y^2 + z^2)^2 - 50*J2*R^2*mu*z^2*(x^2 + y^2 + z^2)^2 + 63*J3*R^3*mu*x^2*z*(x^2 + y^2 + z^2) + 63*J3*R^3*mu*y^2*z*(x^2 + y^2 + z^2) - 35*J2*R^2*mu*x^2*z^2*(x^2 + y^2 + z^2) - 35*J2*R^2*mu*y^2*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2)), 0, 0, 0, -(3*R^2*mu*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)), (R^3*mu*(3*x^4 + 6*x^2*y^2 - 24*x^2*z^2 + 3*y^4 - 24*y^2*z^2 + 8*z^4))/(2*(x^2 + y^2 + z^2)^(9/2))];
Amat(7,:) = [0, 0, 0, 0, 0, 0, 0, 0];
Amat(8,:) = [0, 0, 0, 0, 0, 0, 0, 0];
Phi_dot = Amat*Phi;

Phi_dot_flat = reshape(Phi_dot,64,1);
zd = [xdot_vec;0;0;Phi_dot_flat];
end