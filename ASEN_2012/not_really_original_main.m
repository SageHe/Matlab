clear all;close all;clc
%% Water Exhaust phase
% global g C_d rho_air_amb Vol_bottle p_amb gamma rho_water D_throat D_bottle A_bottle R M_Bottle C_D p_gage Vol_water_i T_air_i v0 theta x0 y0 l_s tspan init_air_v init_water_mass init_air_mass A_throat bottle_air_mass counter
% 
%  
% 
% 
g = 9.81; % m/s^2 ... acceleration due to gravity
C_d = 0.8; % discharge coefficient
rho_air_amb = 0.961; % kg/m^3 ... ambient air density
Vol_bottle = 0.002; % m^3 ... volume of empty bottle
p_amb = 83426.56; %pa ... atmospheric pressure
gamma = 1.4; % ratio of specific heats for air
rho_water = 1000; % kg/m^3 ... density of water
D_throat = .021; % m ... diameter of throat
D_bottle = .105; % m ... diameter of bottle
A_bottle = pi*(D_bottle/2)^2;
R = 287; % J/kgK ... gas constant of air
M_Bottle = 0.15; % kg ... mass of empty 2-liter bottle with cone and fins
C_D = 0.5; % drag coefficient
p_gage = 344738; % pa ... initial gage pressure of air in bottle
Vol_water_i = 0.001; % m^3 ... initial volume of water inside bottle
T_air_i = 300; % K ... initial temperature of ai
v0 = 0.0; % m/s ... initial velocity of rocket
theta = 45; % ° ... initial angle of rocket
x0 = 0.0; % m ... initial horizontal distance
y0 = 0.01; % m ... initial vertical height
l_s= 0.5; % m ... length of test stand
tspan = [0 5]; %s ... integration time
init_air_v = Vol_bottle - Vol_water_i;
init_water_mass = Vol_water_i*rho_water;
init_air_mass = ((p_gage + p_amb)*init_air_v)/(R*T_air_i);
A_throat = pi*(D_throat/2)^2;
bottle_air_mass = Vol_bottle*rho_air_amb;

%Specify initial conditions
%y_0 = [Vx, Vz, Hx, Hz, init_air_v, water_mass, init_air_mass] 
y_0 = [0 0 0 0.01 init_air_v init_water_mass init_air_mass];

%Options = odeset('Events', @myEvent);
[t,y] = ode45(@odefun2,tspan,y_0);

for i = 1:length(t)
    [~,thrust(i)] = odefun2(t(i), y(i,:));
end

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
plot(t,thrust);
axis([0 0.45 0 200])
title('Thrust Profile')
xlabel('Time [s]');
ylabel('Thrust [N]');