%% HW 7
%% Part B
clear;close all;clc
%Define given constants
dt = 0.5;
Omega = 0.045;
%Define given F matrix
F = [1 (sin(Omega*dt)/Omega) 0 -((1-cos(Omega*dt))/Omega);
    0 cos(Omega*dt) 0 -sin(Omega*dt);
    0 (1-cos(Omega*dt))/Omega 1 sin(Omega*dt)/Omega
    0 sin(Omega*dt) 0 cos(Omega*dt)];
%Define ICs
mu0 = [0 85*cos(pi/4) 0 -85*sin(pi/4)]';
P0 = diag([10 2 10 2]);
mu = [];
sigma = [];
%Calculate mu and P for the specified number of steps
for i = 1:300
   mu = [mu F^(i)*mu0];
   P = F^(i)*P0*(F')^(i);
   sigma = [sigma;[P(1,1) P(2,2) P(3,3) P(4,4)]];
end
t = linspace(0,150,size(mu,2));
sigma = 2*sqrt(sigma);
sigma = sigma';
figure(1)
subplot(4,1,1)
hold on
plot(t,mu(1,:))
plot(t,mu(1,:)+sigma(1,:),'--')
plot(t,mu(1,:)-sigma(1,:),'--')
grid on
grid minor
xlabel('Time (s)')
ylabel('East Position (m)')
title('East Position VS Time')
subplot(4,1,2)
hold on
plot(t,mu(2,:))
plot(t,mu(2,:)+sigma(2,:),'--')
plot(t,mu(2,:)-sigma(2,:),'--')
grid on
grid minor
xlabel('Time (s)')
ylabel('East Velocity (m/s)')
title('East Velocity VS Time')
subplot(4,1,3)
hold on
plot(t,mu(3,:))
plot(t,mu(3,:)+sigma(3,:),'--')
plot(t,mu(3,:)-sigma(3,:),'--')
grid on
grid minor
xlabel('Time (s)')
ylabel('North Position (m)')
title('North Position VS Time')
subplot(4,1,4)
hold on
plot(t,mu(4,:))
plot(t,mu(4,:)+sigma(4,:),'--')
plot(t,mu(4,:)-sigma(4,:),'--')
grid on
grid minor
xlabel('Time (s)')
ylabel('North Velocity (m/s)')
title('North Velocity VS Time')
%Plot only 2sigma vals 
figure
subplot(4,1,1)
plot(t,sigma(1,:))
grid on
grid minor
xlabel('Time (s)')
ylabel('2\sigma (m)')
title('2\sigma VS Time')
subplot(4,1,2)
plot(t,sigma(2,:))
grid on
grid minor
ylim([2.8 2.84])
xlabel('Time (s)')
ylabel('2\sigma (m/s)')
title('2\sigma VS Time')
subplot(4,1,3)
plot(t,sigma(3,:))
grid on
grid minor
xlabel('Time (s)')
ylabel('2\sigma (m)')
title('2\sigma VS Time')
subplot(4,1,4)
plot(t,sigma(4,:))
grid on
grid minor
ylim([2.8 2.84])
xlabel('Time (s)')
ylabel('2\sigma (m/s)')
title('2\sigma VS Time')
%% Part C 
close all;clc
%Define given Omega values for a and b
Omega_a = .045;
Omega_b = -.045;
%Define given initial conditions
mua0 = [0 85*cos(pi/4) 0 -85*sin(pi/4)]';
Pa0 = diag([10 4 10 4]);
mub0 = [3200 85*cos(pi/4) 3200 -85*sin(pi/4)]';
Pb0 = diag([11 3.5 11 3.5]);
%Defind Fa and Fb matrices
Fa = [1 (sin(Omega_a*dt)/Omega_a) 0 -((1-cos(Omega_a*dt))/Omega_a);
    0 cos(Omega_a*dt) 0 -sin(Omega_a*dt);
    0 (1-cos(Omega_a*dt))/Omega_a 1 sin(Omega_a*dt)/Omega_a
    0 sin(Omega_a*dt) 0 cos(Omega_a*dt)];

Fb = [1 (sin(Omega_b*dt)/Omega_b) 0 -((1-cos(Omega_b*dt))/Omega_b);
    0 cos(Omega_b*dt) 0 -sin(Omega_b*dt);
    0 (1-cos(Omega_b*dt))/Omega_b 1 sin(Omega_b*dt)/Omega_b
    0 sin(Omega_b*dt) 0 cos(Omega_b*dt)];
%Define given constant xi_R and eta_R values
xi_R = 100;
eta_R = 100;
%Run simulation for 150 seconds
S = [1 0 0 0;0 0 1 0];
clear mu
clear P
P = [];
prob = [];
for i = 1:300
    mu = S*Fa^(i)*mua0 - S*Fb^(i)*mub0;
    mu = mu';
    P = Fa^(i)*Pa0*(Fa^(i))' + Fb^(i)*Pb0*(Fb^(i))';
    P = [P(1,1) P(3,3)];
    prob = [prob mvncdf([-100 -100],[100 100],mu,P)];
end
figure
plot(t,prob)
grid on
grid minor
xlabel('Time (s)')
ylabel('Probability')
title('Probablility of Collision VS Time')