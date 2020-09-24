%ASEN 3128 Homework 1
%Author:Sage Herrin, SID 106071909
%Created:1/25/19
%Modified:1/27/19
%Use ode45 to model flight of a golf ball with initial velocity and
%position and takes in to account drag
clear all;
close all;
clc
%Declare variables to be used and make global
global g cd m D N_0 E_0 D_0 vN_0 vE_0 vD_0 rho area q_0 V_E0 R_E W 
g = 9.81; %m/s
cd = 0.6; %drag coeff
m = .03; %mass in kg
D = 3e-2; %diameter of ball
V_E0 = [0 20 -20];
% vN_0 = 0; %initial N component of velocity
% vE_0 = 20; %initial E component of velocity in m/s
% vD_0 = -20; %initial D component of velocity in m/s
R_E = [0 0 0];
% N_0 = 0; %initial N position
% E_0 = 0; %initial E position
% D_0 = 0; %initial D position
W = [0 0 0];
rho = 1.225; %air density
area = pi*(D/2)^2;
KE = .5*m*(norm(V_E0)^2); % Compute initial kinetic energy to reamain constant throughout simulation 
m = .01;
distance = []; %initialize distance vector
mass = []; %initialize mass vector 
for i = 1:10
    V = (KE*2)/m; %equal to v^2
    V = sqrt(V/2); %velocity corresponding to current mass to keep kinetic energy constant 
    V_E0 = [0 V -V]; %Inertial velocity vector
    inish_condish = [V_E0 R_E]; %initial velocity and position of golf ball
    tspan = [0 5];
    options = odeset('MaxStep', 10^-1);
    [t, y] = ode45('golfball_fun',tspan,inish_condish,options); %call ode45 to numerically integrate differential equations

    ground = find(-y(:,6)<0); %find index of first vertical position less than 0
    ground = ground(1) - 1; %define ground as last vertical position before it becomes negative, equal to when ball hits ground


    vn = y(1:ground,1); %N velocity component of ball
    ve = y(1:ground,2); %E velocity component of ball
    vd = -y(1:ground,3); %D velocity component of ball
    n = y(1:ground,4); %N position of ball
    e = y(1:ground,5); %E position of ball
    d = -y(1:ground,6); %D position of ball


    figure(1) %Plot the flight path of the ball with variation in mass
    hold on
    plot3(n,e,d) 
%     axis equal
    xlabel('N [M]')
    ylabel('E [M]')
    zlabel('-D [M]')
    title('Variation of Flight Path Due to Mass With Zero Wind')
    legend('0.01 Kg','0.02 Kg','0.03 Kg','0304 Kg','0.05 Kg','0.06 Kg','0.07 Kg','0.08 Kg','0.09 Kg','0.1 Kg')
    mass = [mass m];
    m = m + .01; %iteratively increases mass each time ode45 is called 
    distance = [distance e(end)];
end
m = .03; 
figure(3) %Plot the maximum distance traveled by the ball vs mass of the ball
grid on
grid minor
plot(mass,distance)
title('Mass of Ball VS Distance Traveled in E direction')
xlabel('Mass [Kg]')
ylabel('Distance in E direction [M]')
deflection = []; %initialize deflection and wind vectors  
wind_vec = [];
V_E0 = [0 20 -20]; %reinitialize velocity vector
for i = 1:10
    inish_condish = [V_E0 R_E]; %give initial conditions 
    tspan = [0 5];
    options = odeset('MaxStep', 10^-1);
    [t, y] = ode45('golfball_fun',tspan,inish_condish,options);

    ground = find(-y(:,6)<0); %define when ball hits the ground 
    ground = ground(1) - 1;


    vn = y(1:ground,1); %separate velocity and position into components 
    ve = y(1:ground,2);
    vd = -y(1:ground,3);
    n = y(1:ground,4);
    e = y(1:ground,5);
    d = -y(1:ground,6);


    figure(2) %plot flight path of ball with variation in wind
    hold on
    plot3(n,e,d) 
%     axis equal
    xlabel('N [M]')
    ylabel('E [M]')
    zlabel('-D [M]')
    title('Variation of Flight Path Due to Change in Wind')
    legend('0 m/s','2 m/s','4 m/s','6 m/s','8 m/s','10 m/s','12 m/s','14 m/s','16 m/s','18 m/s')
    wind_vec = [wind_vec W(1)];
    deflection = [deflection n(end)];
    W(1) = W(1) + 1;
end
chenge = gradient(deflection,wind_vec);
figure(4)
plot(wind_vec,chenge)
title('Deflection Rate in M per M/s of wind')
xlabel('Wind [M/s]')
ylabel('Deflection Rate [M per M/s]')
