clear all;close all;clc

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('Project2_Prob2_truth_traj_50days.mat')

tspan = Tt_50;
state = Xt_50(1,1:7);
n = 7;

Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);

Z = [state';Phi_flat];

AU = 149597870.700; %AU in km
Phi = 1357;
c = 299792.458;
const.n = 7;
const.mu = 3.98600433e5;
const.mu_s = 132712440017.987;
const.Am = 0.01*1e-6;
const.P_Phi = Phi/c;
const.JD_o = 2456296.25;
const.AU = AU;

[time,X] = ode45(@(t,Z) kepler_wPhi_mu_SRP_3BP_ODE(t,Z,const),tspan,Z,opts);
