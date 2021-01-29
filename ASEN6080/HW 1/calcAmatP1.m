clear all;close all; clc
% Taking gravitational potential function from Roy P. 506, 
%U = (Gm/r)*[1 - (sum. n=2 to inf.)*(R/r)^n*J_n*P_n(sin(phi))]
%Define symbolic variables
syms x y z xdot ydot zdot mu J2 J3 r R
%define potential function, r, and state
r = sqrt(x^2 + y^2 + z^2);
U = (mu/r)*(1 - ((J2*R^2)/(2*r^4))*(3*z^2 - r^2) - ((J3*R^3*z)/(2*r^6))*(5*z^2 - 3*r^2));
state = [x y z xdot ydot zdot mu J2 J3];
%calculate accelerations
rdd = jacobian(U,[x,y,z]);

f = [xdot ydot zdot rdd 0 0 0];

Amat = simplify(jacobian(f,state));



