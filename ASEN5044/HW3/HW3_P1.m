%% Problem 1
clear all;close all;clc
k = 398600;
r0 = 6678;
w0 = sqrt(k/(r0^3));

A = [0 1 0 0;w0^2+(2*k/(r0^3)) 0 0 2*r0*w0;0 0 0 1;0 (-2*w0/r0) 0 0];
B = [0 0;1 0;0 0;0 1/r0];
C = [1 0 0 0;0 0 1 0];
D = [0 0;0 0];

Sysc = ss(A,B,C,D);

Sysd = c2d(Sysc,10,'zoh');

[F,G,H,M] = ssdata(Sysd);