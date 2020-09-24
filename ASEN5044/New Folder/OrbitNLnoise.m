%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theodore Trozinski 103308050
% ASEN 5044 Estimation 
% Project, Report 1
% Created: April 14, 2020
% Purpose: ODE45 fcn for problem 2D Non-Linear Orbit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUT] = OrbitNLnoise(t,x,mu,noise)
% Input
x = x';
% Prelim Calculashuns
r = sqrt(x(1)^2 + x(3)^2);
% Diff Eq
dx = [x(2);
      -mu*x(1)/r^3;
      x(4);
      -mu*x(3)/r^3];
dx(2) = dx(2) + noise(1);
dx(4) = dx(4) + noise(2);
% Set Ouptut
OUT = dx;
end