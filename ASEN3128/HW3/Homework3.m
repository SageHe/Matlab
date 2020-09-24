%Assignment 3, linearized model of quadcopter dynamics and implementation
%of control laws
clear all;clc
global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4
%Constants for every scenario
g = 9.81; %acceleration due to gravity, m/s^2
m = .068; %kg
rad = .06; %m
I_x = 6.8e-5; %Components of moment of inertia matrix
I_y = 9.2e-5;
I_z = 1.35e-4;
k = .0024;
alpha = 2e-6;
eta = 1e-3;
%linearized model for steady hover trim state
f = -m*g;
f1 = f/4;
f2 = f/4;
f3 = f/4;
f4 = f/4;
%Treat initial conditions as disturbances in "top 8" equations
% d_Lc = 0; %L_c moment pert.
% d_Mc = 0; %M_c moment pert.
% d_Nc = 0; %N_c moment pert.
% d_p =0; % y_axis angular velocity
% d_q = 0; %z_axis angular velocity
% d_phi = deg2rad(5); % bank angle pert.
% d_theta = 0; %elevation angle pert.
% d_Zc = 0; %aero_control force pert. in Z-axis direction
dp = 0;
dq = 0;
dr = 0;
dphi = deg2rad(5);
dtheta = 0;
du = 0;
dv = 0;
dw = 0;
% inish_condish = [d_Lc d_Mc d_Nc d_p d_q d_phi d_theta d_Zc];
inish_condish = [dp dq dr dphi dtheta du dv dw];
tspan = [0 5];
[t y] = ode45('hw3ode',tspan,inish_condish);
figure
hold on
plot(t,y(:,6))
plot(t,y(:,7))
plot(t,y(:,8))
