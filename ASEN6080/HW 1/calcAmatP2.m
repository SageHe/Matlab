clear all;close all;

syms x y z xdot ydot zdot mu J2 r R
%Define r and gravitational potential function for J2
r = sqrt(x^2 + y^2 + z^2);
U = (mu/r)*(1 - ((J2*R^2)/(2*r^4))*(3*z^2 - r^2));
%Define state vector
state = [x y z xdot ydot zdot J2];
%Calculate acclerations
rdd = jacobian(U,[x,y,z]);

f = [xdot ydot zdot rdd 0];
%Calculate A matrix
Amat = simplify(jacobian(f,state));