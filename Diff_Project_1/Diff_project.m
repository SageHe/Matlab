%Question 1
clear all;close all;clc
T = input('Amount of time in years interest is being compounded over \n');
t = 0:.1:T;
A = compound(4,t,750000);
B = compound(12,t,750000);
C = compound_cont(.03,t,750000);

hold on 

plot(t,A)
plot(t,B)
plot(t,C)

%Question 2
%Set yprime equal to zero and solve for y to determine equilibium
%solutions => y = (12p/r) 

%This equilibrium equation is unstable since all solutions near the equilibrium equation deviate away from the equilibrium solution 
