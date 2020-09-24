% %Assignment 3, linearized model of quadcopter dynamics and implementation
% %of control laws
clear all;clc
global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4 K1_roll K2_roll K1_pitch K2_pitch K_y
% %Constants for every scenario
g = 9.81; %acceleration due to gravity, m/s^2
m = .068; %kg
rad = .06; %m
I_x = 6.8e-5; %Components of moment of inertia matrix
I_y = 9.2e-5;
I_z = 1.35e-4;
k = .0024;
alpha = 2e-6;
eta = 1e-3;
K1_roll = .0017;
K2_roll = .0106;
K1_pitch = .0023;
K2_pitch = .0144;
K_y = .004;
% %linearized model for steady hover trim state
f = -m*g;
f1 = f/4;
f2 = f/4;
f3 = f/4;
f4 = f/4;
% %Treat initial conditions as disturbances in "top 8" equations
% d_Lc = 0; %L_c moment pert.
% d_Mc = 0; %M_c moment pert.
% d_Nc = 0; %N_c moment pert.
% d_p =0; % y_axis angular velocity
% d_q = 0; %z_axis angular velocity
% d_phi = deg2rad(5); % bank angle pert.
% d_theta = 0; %elevation angle pert.
% d_Zc = 0; %aero_control force pert. in Z-axis direction
dp = 0;
dq = 0.1;
dr = 0;
dphi = 0;
dtheta = 0;
du = 0;
dv = 0; 
dw = 0;

inish_condish = [dp dq dr dphi dtheta du dv dw];
tspan = [0 5];
[t y] = ode45('hw4ode',tspan,inish_condish);
figure
subplot(1,3,1)
hold on
grid on
grid minor
plot(t,y(:,1))
plot(t,y(:,2))
plot(t,y(:,3))
title('Roll Pitch and Yaw Rates Vs Time for dq = 0.1 rad/s, Linear Model')
xlabel('Time (s)')
ylabel('Rad/s')
legend('Roll Rate','Pitch Rate','Yaw Rate')
% figure 
subplot(1,3,2)
hold on
grid on
grid minor
plot(t,y(:,4))
plot(t,y(:,5))
title('Bank and Pitch Angle Vs Time for dq = 0.1 rad/s, Linear Model')
xlabel('Time (s)')
ylabel('Radians')
legend('\phi','\theta')
subplot(1,3,3)
% figure
hold on
grid on
grid minor
plot(t,y(:,6))
plot(t,y(:,7))
plot(t,y(:,8))
title('Velocity Components Vs Time for dq = 0.1 rad/s, Linear Model')
xlabel('Time(s)')
ylabel('m/s')
legend('u','v','w')
%Data plotting from lab
data1 = load('W_1138_A4');

figure
grid on
grid minor
plot3(data1.rt_estim.signals.values(:,1),data1.rt_estim.signals.values(:,2),-data1.rt_estim.signals.values(:,3))
title('Spiderbot Copter Position For Implemented Controls ')
xlabel('N (m)')
ylabel('E (m)')
zlabel('D (m)')
legend('Spider Copter Position')
% subplot(1,2,1)
figure
grid on
grid minor
hold on
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,10))
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,11))
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,12))
title('Spider Copter Roll Pitch and Yaw Rates Vs Time, Implemented Controls')
xlabel('time (s)')
ylabel('rad/s')
legend('p','q','r')
figure
grid on
grid minor
hold on
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,4));
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,5));
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,6));
title('Spider copter Euler Angles Vs Time, Implemented Controls')
xlabel('Time (s)')
ylabel('Radians')
legend('\psi','\theta','\phi')
