clear all;close all;clc 
%% Problem 1 part C
% part i
% Define constants
 J2 = 1e-3;
 mu = 4e5;
 Re = 6400;
 rho = (2*pi)/(365*24*3600);
 temp = ((2*rho)/(3*J2*Re^2*sqrt(mu)))^(2/7);
 a_M = 1/temp;
 % pot curve of relation showing inc. for SSO as a func. of orbit radius
 a = linspace(6400,a_M,100);
 i = acos(-(a./a_M).^(7/2));
 
 figure
 plot(a,i)
 xlabel('Orbit Radius (Km)')
 ylabel('Inclination (rad)')
 title('Orbit Radius VS Inclination')
 grid on
 grid minor
 
 