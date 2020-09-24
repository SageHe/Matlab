%% House Keeping
clear all;
close all;
clc

%% Define constants

%% Enter Equation for rigid arm
% Assign coeffcient 
Kd = 0;
Kp = 200;
WD = 2.5;

% Natural Frequency

% Zeta
% Set up and build transfer function
num = Kp;
den = [ 0 1 (2+Kp) ];

sysTF = tf(num,den);

[x t] = step(sysTF);

x = WD*x;

figure(1); clf;
hold on
% y = ones(1,length(t))*(2*WD);
% plot(t,y);
plot(t,x);

u_p = -Kp.*(x - WD);
%%PI control
Kp = 200;
Ki = 200;
WD = 2.5;

num_int = [Kp Ki];
den_int = [1 (2+Kp) Ki];

int_sysTF = tf(num_int,den_int);
[x_int t_int] = step(int_sysTF);

x_int = WD*x_int;
figure(2); clf;
plot(t_int,x_int)
W_ss = ((Kp*WD)/(2+Kp));
u_i = -Kp.*(x - WD) - Ki.*(W_ss - WD).*t_int;

figure(3)
hold on
plot(t,u_p)
plot(t_int,u_i)
title('P Vs PI Control Law')
xlabel('Time')
ylabel('Control')
legend('Proportional Control','Proportional Integral Control')
