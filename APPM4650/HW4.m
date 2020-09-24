%Sage Herrin, Numerical Analysis 1, Prof. Coffey
clear all;close all;clc
%% Solutions to question 3.5, 1
%First solution
A = [1 0 0 0 0 0 0 0;1 1 1 1 0 0 0 0;0 0 0 0 1 0 0 0;0 0 0 0 1 1 1 1;...
    0 1 2 3 0 -1 0 0;0 0 2 6 0 0 -2 0;0 0 2 0 0 0 0 0;0 0 0 0 0 0 2 6];

b = [0 1 1 2 0 0 0 0]';

x = A\b

%second solution
syms a0 b0 c0 d0 a1 b1 c1 d1
eqn1 = a0 == 0;
eqn2 = a0 + b0 + c0 + d0 == 1;
eqn3 = a1 == 1;
eqn4 = a1 + b1 + c1 + d1 == 2;
eqn5 = b0 + 2*c0 + 3*d0 - b1 == 0;
eqn6 = 2*c0 + 6*d0 - 2*c1 == 0;
eqn7 = 2*c0 == 0;
eqn8 = 2*c1 + 6*d1 == 0;

[s1 s2 s3 s4 s5 s6 s7 s8] = solve(eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8,a0, b0, c0, d0, a1, b1, c1, d1);

S = [s1 s2 s3 s4 s5 s6 s7 s8]

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
xp1(end)

yp1(end)

