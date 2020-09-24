%Housekeeping
clear all;close all;clc

%Sage Herrin, ASEN 2012 Homework 7, Due 12/4/17
%% Problem 1
%Determine poles of transfer function that represents navigational control system for aircraft to three sig figs
syms S
A = vpasolve(0 == 3*S^6 + 7*S^4 - 27*S^3 + 13*S^2 - S + 4);
%print values to outuput file and round to three sig figs
fileID = fopen('Homework_7.txt','w');
fprintf(fileID,'Sage Herrin, SID #106071909, Homework #7, due 12/4/17 \n');
fprintf(fileID,'Problem 1 pole values of transfer function: \n');
fprintf(fileID,'%.3f \n',A(1));
fprintf(fileID,'%.3f \n',A(2));
fprintf(fileID,'%.3f \n',A(3));
fprintf(fileID,'%.3f \n',A(4));
fprintf(fileID,'%.3f \n',A(5));
fprintf(fileID,'%.3f \n',A(6));

%% Problem 2
%Determine average surface temperature for given case to within 4 sig figs
%of accuracy

heat_xfer = 750; %desired rate of heat flux removed from the hot metal
h = 1;
T_surround = 110; %Degrees F
epsilon = 0.93; %emmisivity
sigma = 1.174*10^-9; %Stefan-Boltzman constant
syms T_s
avg_surface_temp = vpasolve(heat_xfer == h*(T_s - T_surround) + epsilon*sigma*(T_s^4 - T_surround^4));
%print avg_surface_temp to text file to 4 sig figs
fprintf(fileID,'Problem 2 Average surface temperature values in Fahrenheit: \n');
fprintf(fileID,'%.4f \n',avg_surface_temp(1));
fprintf(fileID,'%.4f \n',avg_surface_temp(2));
fprintf(fileID,'%.4f \n',avg_surface_temp(3));
fprintf(fileID,'%.4f \n',avg_surface_temp(4));

%% Problem 3
%solve system of linear eqns using Gaussian elimination without pivoting.
%Check answers using \ operator and with inverse operator (^-1) for A mat
%and time results and compare.
%Linear eqns
%A mat
A = [7 3 17;-4 0 2;4 3 -9];
b = [13;-2;-5];
tic
%Determine size of A and b mat
[m,m] = size(A);
[j,k] = size(b);
z = zeros(j,k);
for i = 1:m - 1
    n = -A(i+1:m,i)/A(i,i);
    A(i+1:m,:) = A(i+1:m,:) + n*A(i,:);
    b(i+1:m,:) = b(i+1:m,:) + n*b(i,:);
end
%Back sub
x(m,:) = b(m,:)/A(m,m);
for j = m-1:-1:1
    x(j,:) = (b(j,:) - A(j,j+1:m)*x(j+1:m,:))/A(j,j);
end
t1 = toc;
fprintf(fileID,'Problem 3 X,Y, and Z coefficients, respectively, solved for using Gaussian substitution without pivoting: \n');
fprintf(fileID,'%.4f \n',x(1));
fprintf(fileID,'%.4f \n',x(2));
fprintf(fileID,'%.4f \n',x(3));

%Compare to backslash method and time
tic
x1 = A\b;
t2 = toc;

fprintf(fileID,'X,Y, and Z coefficients, respectively, solved for using backslash method with time elapsed: \n');
fprintf(fileID,'%.4f \n',x1(1));
fprintf(fileID,'%.4f \n',x1(2));
fprintf(fileID,'%.4f \n',x1(3));
fprintf(fileID,'Time ellasped to perform calculation in seconds: \n');
fprintf(fileID,'%.6f \n',t2);

%Compare to inverse method
tic
x2 = A^-1*b;
t3 = toc;
    
fprintf(fileID,'X,Y, and Z coefficients, respectively, solved for using inverse method with time elapsed: \n');
fprintf(fileID,'%.4f \n',x2(1));
fprintf(fileID,'%.4f \n',x2(2));
fprintf(fileID,'%.4f \n',x2(3));
fprintf(fileID,'Time ellasped to perform calculation in seconds: \n');
fprintf(fileID,'%.6f \n',t3);
    
    
    
