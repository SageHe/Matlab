% Script utilies ode45 function to calculate trajectory and thrust of
% bottle rocket along two dimensions (x and z).
% 
% Inputs:
% none, though users should set initial conditions in both this script and
% in rocketODE() as needed
%
% Outputs:
% Graphs of rocket trajectory and thrust profile
% Max height and distance traveled
%
% Author: 212
% Created: 12/2/17
% Modified: 12/2/17

% Housekeeping
clear; clc; close all;
global counter
counter = 1;
%% Constants 
% Environmental
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('rocketODE', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = rocketODE(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
%hold on
%plot(y(y(:,4) == max(y(:,4)),3), max(y(:,4)),'ko','MarkerFaceColor',[0,0,0]);
axis([0 70 0 20])
title('Bottle Rocket Trajectory');
xlabel('Distance [m]');
ylabel('Height [m]');

% Determine max distance
pre_dist = sum(y(:,4) >= 0);
post_dist = pre_dist + 1;
max_dist = linterp([y(post_dist,4) y(pre_dist,4)], [y(post_dist,3) y(pre_dist,3)], 0);

% Label Trajectory plot
maxH = sprintf('Max Height = %.2fm', max(y(:,4)));
maxD = sprintf('Max Distance = %.2fm', max_dist);
text(20,19, maxH);
text(20,18, maxD);

% Thrust plot
subplot(1,2,2);
plot(t,F);
axis([0 0.45 0 200])
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');
