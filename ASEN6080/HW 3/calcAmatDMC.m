clear all;close all;clc

syms x y z xdot ydot zdot mu J2 r R Bw tau
%Define r and gravitational potential function for J2
% Amat = zeros(9,9);

r = sqrt(x^2 + y^2 + z^2);
U = (mu/r)*(1 - ((J2*R^2)/(2*r^4))*(3*z^2 - r^2));
%Define state vector
state = [x y z xdot ydot zdot];
%Calculate acclerations
rdd = jacobian(U,[x,y,z]);

f = [xdot ydot zdot rdd];
Bw = (1/tau)*eye(3);
%Calculate A matrix
Amat(1:6,1:6) = simplify(jacobian(f,state));
Amat(7:9,1:3) = zeros(3);
Amat(1:3,7:9) = zeros(3);
Amat(4:6,7:9) = eye(3);
Amat(7:9,7:9) = -Bw;
