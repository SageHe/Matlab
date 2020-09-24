%Project number 142
%Date created: 11/21/17
%Date modified: 12/8/17

% In this script, ODE45 is used to carry out numerous 4th order 5th degree
% accurate Runge Kutta approximations to compute solutions to diffeential
% equations in order to track the thrust and trajectory of bottle rocket
% filled with water
% 
% Inputs:
% Varying initial conditions including water amaount, initial pressure,
% launch angle, and several others.
%
% Outputs:
% Plot of rocket distance vs height as well as a thrust profile

% Housekeeping
clear; clc; close all;
%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 566059;% Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4)); 
title('Bottle Rocket Trajectory, 70 PSI Gage Pressure');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 60 PSI Gage Pressure');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 359216.9; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 40 PSI Gage Pressure');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, Test Case');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.00075; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, .75 Kg of Water');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values
v_air_w = 0.00125; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 1.25 Kg of Water');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.00135; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 1.35 Kg of Water');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2_3_C_D', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 0.3 C_D');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2_4_C_D', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 0.4 C_D');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2_6_C_D', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 0.6 C_D');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all
figure

%% Given Values
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2_30_Deg', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 30 Degree Launch Angle');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2_60_Deg', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 60 Degree Launch Angle');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all;
figure

%% Given Values 
v_air_w = 0.001; % initial water volume [m^3]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2_65_Deg', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 65 Degree Launch angle');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');

clear all
figure

%% Given Values 
v_air_w = 0.00135; % initial water volume [m^3]
P_air_i = 566059; % Initial air pressure, gage + ambient [Pa]
T_air_i = 300; % air temperature [k]

R = 287; % Specific gas constant [J/(kg*K)]
v_air_i = .002 - v_air_w; % initial air volume [m^3]
mass_w = 1; % mass water [kg]
mass_air_i = (P_air_i/(R*T_air_i))*v_air_i; % initial mass of air [kg]

%% Ode45 calculation
% y = [V_x, V_z, dist, height, V_air, mass_water, mass_air]
[t,y] = ode45('odefun2', [0 5], [0 0 0 .01 v_air_i mass_w mass_air_i]);

% Solve for thrust values at each time point
for i = 1:length(t)
    [~,F(i)] = odefun2(t(i), y(i,:));
end

%% Plots
% Trajectory plot
subplot(1,2,1);
plot(y(:,3), y(:,4));
title('Bottle Rocket Trajectory, 75 Meter distance case');
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
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');