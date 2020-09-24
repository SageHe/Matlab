%Adam Boylston, Sage Herrin
%% Housekeeping
close all
clear all
clc

%% 4

syms alpha beta gamma N0

% a

% equilibrium points: (0,0), (0,N0), (N0,0)
equib = [0 0; 0 N0; N0 0]; %first equib of 0 0 may not be an equilibrium point

% b
%partials of dS_dt and dZ_dt with respect to S first, then Z
%J = [(-bet*Z), (-bet*S); (bet*Z - gamma - alpha*Z), (bet*S - gamma - alpha*S)];
for i = 1:length(equib)
    J{1,i} = [(-beta*equib(i,2)), (-beta*equib(i,1)); (beta*equib(i,2) - gamma - alpha*equib(i,2)), (beta*(equib(i,1) - alpha) - gamma)];
end
    

% eig
%Code below was given in project3 pdf
% Assign parameter values
alpha = 0.00002; beta = 0.00003; gamma = 0.000025; N0 = 60000;
% Set length of simulation
tspan = [0 35];
% Set initial conditions
y0 = [59999; 1];
% Solve system
[t,y] = ode45(@(t,y) szr(t,y,alpha,beta,gamma,N0), tspan, y0);

S_t = y(:,1);
Z_t = y(:,2);
R_t = N0 - S_t - Z_t;

figure(1)
plot(t,S_t)
hold on
plot(t,Z_t)
hold on
plot(t,R_t)
title('S, Z, R vs. Time')
xlabel('Time (days)')
ylabel('Humanoids')
legend('Susceptibles','Zombies','Removed')


%% 5.1
% a

rho = 0.2;
%tspan2 = [0 100];
tspan2 = [0 35];
[t2,y2] = ode45(@(t,y) szr_the_sequel(t,y,alpha,beta,gamma,N0,rho), tspan2, y0);
S_t2 = y2(:,1);
Z_t2 = y2(:,2);
R_t2 = N0 - S_t2 - Z_t2;

figure(2)
plot(t2,S_t2)
hold on
plot(t2,Z_t2)
hold on
plot(t2,R_t2)
title('S, Z, R vs. Time with Antidote')
xlabel('Time (days)')
ylabel('Humanoids')
legend('Susceptibles','Zombies','Removed')
% b

%% 5.2

tspan3 = [0 35];
[t3,y3] = ode45(@(t,y) revenge_of_the_szr(t,y,alpha,beta,gamma,N0), tspan3, y0);
S_t3 = y3(:,1);
Z_t3 = y3(:,2);
R_t3 = N0 - S_t3 - Z_t3;

figure(3)
plot(t3,S_t3)
hold on
plot(t3,Z_t3)
hold on
plot(t3,R_t3)
title('S, Z, R vs. Time with Cruise Missile')
xlabel('Time (days)')
ylabel('Humanoids')
legend('Susceptibles','Zombies','Removed')





