%% Problem 4
clear all;close all;clc

k = 398600;
r0 = 6678;
w0 = sqrt(k/(r0^3));

A = [0 1 0 0;w0^2+(2*k/(r0^3)) 0 0 2*r0*w0;0 0 0 1;0 (-2*w0/r0) 0 0];

tspan = [0 10];
y0 = eye(4);
[t,y] = ode45(@(t,y)odefunp4(t,y,A),tspan,y0);

Phi10 = reshape(y(end,:),size(A));

tspan = [0 100];
y0 = eye(4);
[t,y] = ode45(@(t,y)odefunp4(t,y,A),tspan,y0);

Phi100 = reshape(y(end,:),size(A));
