%% Problem 2a, solve given IVP using Euler's method and compare to exact soln
clear all;clc
t = 0:.1:1;
h = .1;
y(1) = 2;
f = @(t,y) -5*y + 6*exp(t);
for i = 1:length(t) - 1
    y(i + 1) = y(i) + h*f(t(i),y(i));
end
% compare estimated values to actual values
y
exact = exp(-5.*t) + exp(t)
error = abs(exact - y)
%% Problem 2b
clear all;clc
t = 0:0.1:1
h = 0.1
y(1) = exp(1);
f = @(t,y) -10*y + 10*t + 1;
for i = 1:length(t) - 1
    y(i+1) = y(i) + h*f(t(i),y(i));
end
y
exact = exp(-10.*t + 1) + t
error = abs(exact - y)
%% Chapter 6, problem 6a
clear all;clc
A = [0 1 -2;1 -1 1;1 0 -1];
determ = det(A)
