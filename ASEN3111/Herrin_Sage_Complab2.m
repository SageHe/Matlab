%Author:Sage Herrin
%Created: 9/21/18
%SID:106071909
%Comp Lab 2
%ASEN 3111, Aerodynamics
%Housekeeping
clear all;close all;clc
%Flow parameters
c = 2;
alpha = 12;
V_inf = 68;
P_inf = 101.3e3;
rho_inf = 1.225;
N = 10000;
%%Write matlab function that plots stream lines, equipotential lines, and
%%pressure contours for flow about a thin symmetric airfoil
[V,P] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N);
%Error in velocity and pressure with variation in N;
N = [2:100];
val_errors = zeros(length(N),2);
for i = 1:length(N)
[V_meas,P_meas] = Plot_Airfoil_Flow_Error(c,alpha,V_inf,P_inf,rho_inf,N(i),V,P);
val_errors(i,1) = (norm(V - V_meas,'fro'))/(norm(V,'fro'));
val_errors(i,2) = (norm(P - P_meas,'fro'))/(norm(P,'fro'));
end
val_errors = val_errors*100;
figure(4)
plot(N,val_errors);
xlabel('Number of Vortices N')
ylabel('Percent Error')
title('Error in Pressure and Velocity VS N')
legend('Velocity','Pressure')
