%% Solution to written question
clear all;close all;clc
%define initial starting points of each function
y = -0.5;
x = 0.5;
xp1 = [];
yp1 = [];
%define given eqns to perfomr 2D NR on
for i = 1:15
    f1 = x + tan(y);
    f2 = y + 2*cos(x) - exp(x);
    %define Jacobian with respective partials
    J = det([1 sec(y)^2;-2*cos(x) - exp(x) 1]);
    %defind h and k matrices used to solve system
    h = (1/J)*det([-f1 sec(y)^2;-f2 1]);
    k = (1/J)*det([1 -f1;-2*sin(x) - exp(x) -f2]);
    %determine next iteration of X values
    xp1 = [xp1 x + h;];
    yp1 = [yp1 y + k];
    x = xp1(end);
    y = yp1(end);
end
xp1 = xp1';
yp1 = yp1';

fprintf