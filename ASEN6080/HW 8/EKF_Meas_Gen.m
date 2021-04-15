clear all;close all;clc
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
rng('default'); %fix random number seed for repeatability
t = linspace(0,2*pi,50);
k = 1;
x0 = [0;1];
x(:,1) = x0;
dt = t(2) - t(1);
n = 2;
eta = 0.1^2*randn(1,numel(t));

%generate truth state and measurements
tspan = linspace(0,2*pi,50);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];
[~,x] = ode45(@(t,Z) EKF_ode(t,Z,n),tspan,Z,opts);
y = x(:,1) + eta';
