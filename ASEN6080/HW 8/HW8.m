%% HW 8 -- Particle filter. Build linear KF and EKF (very basic) to compare to performance of the Particle filter
% Problem 1 -- linear dynamcis, linear measurement -- Implement and run KF
% with SNC, generate requested plots
clear all;close all;clc
rng(1); %fix random number seed for repeatability
t = linspace(0,2*pi,50);
k = 1;
A = [0 1;-1 0];
x0 = [0;1];
x(:,1) = x0;
dt = t(2) - t(1);

P0 = diag([0.2^2 0.2^2]);
P_plus(:,:,1) = P0;
eta = 0.1^2*randn(1,numel(t));
dx_minus(:,1) = zeros(2,1);
dx_plus(:,1) = zeros(2,1);

for i = 2:numel(t)
    Phi = expm(A*dt);
    x(:,i) = Phi*x(:,i-1);
    
    y(i) = x(1,i) + eta(i);
    %time update
%     dx_minus(:,i) = Phi*dx_plus(:,i-1);
%     P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi';
%     %process observations
%     r(i) = y(i) - x(1,i);
end