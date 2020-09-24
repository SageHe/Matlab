%% House Keeping
clear all;
close all;
clc

%% Define constants

%% Enter Equation for rigid arm
% Assign coeffcient 
Kd = 0;
Kp = 2;
WD = 2.5;

% Natural Frequency

% Zeta
% Set up and build transfer function
num = Kp;
den = [ 0 1 (2+Kp) ];

sysTF = tf(num,den);

[x t] = step(sysTF);

x = 2*WD*x;

figure(1); clf;
hold on
% y = ones(1,length(t))*(2*WD);
% plot(t,y);
plot(t,x);

u = -Kp.*(x - WD);