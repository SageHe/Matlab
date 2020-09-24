clear all;close all;clc
%define constants
global ss theta_0 u_0 A
g = 9.81;
S = 510.9667;%m^2
b = 59.6433;%m
c_bar = 8.3241;%m
theta_0 = 0;
W = 2.8317e6;%Newtons
m = W/9.81;
V = 157.886;%m/s
rho = .6532;%kg/m^3
u_0 = V;%m/s
xi = deg2rad(-6.8);
C_w0 = W/(.5*rho*(u_0^2)*S);%weight coeff
I_x = 24675586.69;%kgm^2
I_y = 44877574.145;%kgm^2
I_z = 67384152.115;%kgm^2
I_zx = 1315143.4;%kgm^2
%non-dimensional stability matrix
ndstab = [-0.1080 -0.1060 0.1043;0.2193 -4.920 -1.023;0 -5.921 -23.92;0 5.896 -6.314];
%Create dimensional stability matrix for body frame using above constants 
bs = [rho*u_0*S*C_w0*sin(theta_0) + .5*rho*u_0*S*ndstab(1,1) -rho*u_0*S*C_w0*cos(theta_0) + .5*rho*u_0*S*ndstab(1,2) .5*rho*u_0*c_bar*S*ndstab(1,3);...
        .5*rho*u_0*S*ndstab(2,1) .5*rho*u_0*S*ndstab(2,2) .5*rho*u_0*c_bar*S*ndstab(2,3);...
        .25*rho*u_0*c_bar*S*ndstab(3,1) .25*rho*u_0*c_bar*S*ndstab(3,2) .25*rho*u_0*(c_bar^2)*S*ndstab(3,3);...
        .25*rho*c_bar*S*ndstab(4,1) .25*rho*c_bar*S*ndstab(4,2) .25*rho*(c_bar^2)*S*ndstab(4,3)];
    
%change above matrix from body frame to stability frame using following
%conversion, taken from Appendix B12 of textbook
ss = zeros(4,3);
ss(1,1) = bs(1,1)*cos(xi)^2 - (bs(2,1) + bs(1,2))*sin(xi)*cos(xi) + bs(2,2)*sin(xi)^2;
ss(2,1) = bs(2,1)*cos(xi)^2 + (bs(1,1) - bs(2,2))*sin(xi)*cos(xi) - bs(1,2)*sin(xi)^2;
ss(3,1) = bs(3,1)*cos(xi) - bs(3,2)*sin(xi);
ss(4,1) = -bs(4,2)*sin(xi)*cos(xi);
ss(1,2) = bs(1,2)*cos(xi)^2 - (bs(2,2) - bs(1,1))*sin(xi)*cos(xi) - bs(2,1)*sin(xi)^2;
ss(2,2) = bs(2,2)*cos(xi)^2 + (bs(1,2) + bs(2,1))*sin(xi)*cos(xi) + bs(1,1)*sin(xi)^2;
ss(3,2) = bs(3,2)*cos(xi) + bs(3,1)*sin(xi);
ss(4,2) = bs(4,2)*cos(xi)^2;
ss(1,3) = bs(1,3)*cos(xi) - bs(2,3)*sin(xi);
ss(2,3) = bs(2,3)*cos(xi) + bs(1,3)*sin(xi);
ss(3,3) = bs(3,3);
ss(4,3) = bs(4,3)*cos(xi);
%the above stability matrix converted to the stability frame is used to
%construct the A matrix of the linearized longitudinal model
A = [(ss(1,1)/m) (ss(2,1)/m) 0 -g*cos(theta_0);...
    (ss(1,2)/(m - ss(4,2))) (ss(2,2)/(m - ss(4,2))) ((ss(3,2) + m*u_0)/(m - ss(4,2))) ((-m*g*sin(theta_0))/(m - ss(4,2)));...
    ((1/I_y)*(ss(1,3) + ((ss(4,3)*ss(1,2))/(m - ss(4,2))))) (1/I_y)*(ss(2,3) + ((ss(4,3)*ss(2,2))/(m - ss(4,2))))...
    ((1/I_y)*(ss(3,3) + ((ss(4,3)*(ss(3,2) + m*u_0)/(m - ss(4,2)))))) (-(ss(4,3)*m*g*sin(theta_0))/(I_y*(m - ss(4,2))));...
    0 0 1 0];
    
%determine the eigen values and eigenvectors of A using built in eig
%function
[V,D] = eig(A);

%Calculate the natural frequencies and damping ratios using built in
%function damp
[Wn,Z] = damp(D);

%compare the eigenvalues determined above to the approx given in class and
%compare oscillation period of phugoid mode above to Lanchester approx. in
%textbook
%given approx of short period mode evalues
eval_SP_1 = (ss(3,3)/(2*I_y)) + sqrt((ss(3,3))^2 + (4*I_y*u_0*(ss(4,3))))*(1/(2*I_y));
eval_SP_2 = (ss(3,3)/(2*I_y)) - sqrt((ss(3,3))^2 + (4*I_y*u_0*(ss(4,3))))*(1/(2*I_y));

%calculated phugoid period from previous values
T_ph_1 = (2*pi)/0.0938;
%phuoid oscillation period from textbook
T_ph_2 = 0.138*u_0;

%simulate linearized longitudinal dynamics to verify trim state is an
%equilibrium and perturb initial states one at a time
delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 500];
[t,y] = ode45('hw7ode',tspan,inish_condish);

%Verify that trim state is equilibrium, plot output of 6 diff. eqs to show
%that that all are 0 with no perturbations => equilibrium
figure(1)
grid on
grid minor
subplot(3,2,1)
plot(t,y(:,1))
title('\Deltau Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta w Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta q Vs Time')
xlabel('Time (s)')
ylabel('rad/s')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \theta Vs Time')
xlabel('Time (s)')
ylabel('rads')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta X Vs Time')
xlabel('Time (s)')
ylabel('m')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Z Vs Time')
xlabel('Time (s)')
ylabel('m')
sgtitle('No Perturbation')
%Disturbance of delta_u = 10m/s
delta_u = 10; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0; %rads
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 500];
[t,y] = ode45('hw7ode',tspan,inish_condish);

figure(2)
grid on
grid minor
subplot(3,2,1)
plot(t,y(:,1))
title('\Deltau Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta w Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta q Vs Time')
xlabel('Time (s)')
ylabel('rad/s')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \theta Vs Time')
xlabel('Time (s)')
ylabel('rads')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta X Vs Time')
xlabel('Time (s)')
ylabel('m')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Z Vs Time')
xlabel('Time (s)')
ylabel('m')
sgtitle('\Deltau=10 m/s')
%disturbance of delta_w = 10m/s
delta_u = 0; %m/s
delta_w = 10; %m/s
delta_q = 0; %rad/s
delta_theta = 0; %rads
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 500];
[t,y] = ode45('hw7ode',tspan,inish_condish);

figure(3)
grid on
grid minor
subplot(3,2,1)
plot(t,y(:,1))
title('\Deltau Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta w Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta q Vs Time')
xlabel('Time (s)')
ylabel('rad/s')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \theta Vs Time')
xlabel('Time (s)')
ylabel('rads')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta X Vs Time')
xlabel('Time (s)')
ylabel('m')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Z Vs Time')
xlabel('Time (s)')
ylabel('m')
sgtitle('\Deltaw=10 m/s')
%disturbance of delta_q=0.1 rad/s
delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0.1; %rad/s
delta_theta = 0; %rads
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 500];
[t,y] = ode45('hw7ode',tspan,inish_condish);

figure(4)
grid on
grid minor
subplot(3,2,1)
plot(t,y(:,1))
title('\Deltau Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta w Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta q Vs Time')
xlabel('Time (s)')
ylabel('rad/s')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \theta Vs Time')
xlabel('Time (s)')
ylabel('rads')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta X Vs Time')
xlabel('Time (s)')
ylabel('m')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Z Vs Time')
xlabel('Time (s)')
ylabel('m')
sgtitle('\Deltaq=0.1 rads/s')
%disturbance of delta_theta=0.1 rads
delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0.1; %rads
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 500];
[t,y] = ode45('hw7ode',tspan,inish_condish);

figure(5)
grid on
grid minor
subplot(3,2,1)
plot(t,y(:,1))
title('\Deltau Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta w Vs Time')
xlabel('Time (s)')
ylabel('m/s')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta q Vs Time')
xlabel('Time (s)')
ylabel('rad/s')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \theta Vs Time')
xlabel('Time (s)')
ylabel('rads')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta X Vs Time')
xlabel('Time (s)')
ylabel('m')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Z Vs Time')
xlabel('Time (s)')
ylabel('m')
sgtitle('\Delta\theta=0.1 radians')
