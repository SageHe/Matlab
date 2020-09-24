clear all;close all; clc
global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4 X Y Z K1_roll K2_roll K1_pitch K2_pitch K_y feedback
%Constants for every scenario
g = 9.81; %acceleration due to gravity, m/s^2
m = .068; %kg
rad = .06; %m
I_x = 6.8e-5; %Components of moment of inertia matrix
I_y = 9.2e-5;
I_z = 1.35e-4;
k = .0024;
alpha = 2e-6;
eta = 2.5e-1;
K_y = .004;
K1_roll = .0017;
K2_roll = .0106;
K1_pitch = .0023;
K2_pitch = .0144;
feedback = 1;
%Rotor forces
f = -m*g;
f1 = f/4;
f2 = f/4;
f3 = f/4;
f4 = f/4;


clear t y
u_E = 0;
v_E = 0; %m/s
w_E = 0;
X = 0 + -eta*u_E^2;
Y = 0 + -eta*v_E^2;
Z = -m*g + -eta*w_E^2;
p = 0;
q = 0;
r = 0;
phi = 0;
theta = 0;
psi = 0;
N = 0;
E = 0;
D = -1;
vel = [u_E v_E w_E];
omega = [p q r];
euler = [psi theta phi];
pos = [N E D];
phi_d = 0;
theta_d = deg2rad(-5);
inish_condish = [vel omega euler pos phi_d theta_d];
tspan = [0 5];
[t,y] = ode45('HW5ode',tspan,inish_condish);
figure
% subplot(1,3,1)
hold on
grid on
grid minor
plot(t,y(:,1))
plot(t,y(:,2))
plot(t,y(:,3))
title('Velocity Components Vs Time for \theta = -5^{\circ}, Nonlinear Model')
xlabel('Time (s)')
ylabel('m/s')
legend('u','v','w')
figure
hold on
grid on
grid minor
plot(t,y(:,10))
plot(t,y(:,11))
plot(t,y(:,12))
title('Position Vs Time')
xlabel('Time (s)')
ylabel('Meters')
legend('N','E','D')
figure
hold on
grid on
grid minor
plot(t,y(:,4))
plot(t,y(:,5))
plot(t,y(:,6))
title('Roll, Pitch, and Yaw Rate Vs Time')
xlabel('Time (s)')
ylabel('Rads/s')
legend('p','q','r')
figure
hold on
grid on
grid minor
plot(t,y(:,7))
plot(t,y(:,8))
plot(t,y(:,9))
title('Bank Elevation and Azimuth Vs Time')
xlabel('Time (s)')
ylabel('Radians')
legend('\psi','\theta','\phi')

data1 = load('A5_M_1125.mat');

% subplot(1,2,1)
figure
grid on
grid minor
hold on
plot(data1.rt_estim.time(:),-data1.rt_estim.signals.values(:,1))
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,2))
plot(data1.rt_estim.time(:),-data1.rt_estim.signals.values(:,3))
title('Spider Copter Positoin Vs Time')
xlabel('time (s)')
ylabel('Meters')
legend('N','E','D')

