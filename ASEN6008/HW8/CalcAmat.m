clear all;close all;clc

syms x y z xdot ydot zdot mu

R1 = sqrt((x + mu)^2 + y^2 + z^2);
R2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

xdd = 2*ydot + x - (1 - mu)*((x+mu)/(R1^3)) - mu*((x - 1 + mu)/(R2^3));
ydd = -2*xdot + y - (1 - mu)*(y/(R1^3)) - mu*(y/(R2^3));
zdd = -(1 - mu)*(z/(R1^3)) - mu*(z/(R2^3));

f = [xdot ydot zdot xdd ydd zdd];

state = [x y z xdot ydot zdot];

Amat = simplify(jacobian(f,state));
