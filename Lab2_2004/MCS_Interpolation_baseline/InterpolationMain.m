% Housekeeping
clear all
close all
clc
x_el = [];
y_el = [];
%% Modifications to tweak the model
% I made two modifications to how the model should work in theory. 
% 1)
% I adjusted the time scale on how the mass decreases so that the mass change
% and the thrusting happen over the same period of time.
% 2)
% I also selected a range of the thrust data which made the distance and
% altitude come out right. I removed outliers and ommitted negative values.


%% Set up given variables from input document
%Make all variables global
global mw mthermo mtimes F mi mf Py0 windy windx windz g Cd p_air_amb Volbottle P_amb gamma pwater DThroat DBottle R MBottle CD P_gage Volwater_i Tair_i V0 Theta x0 z0 ls AThroat P0 Volair_i m_air_i ABottle counter
drag_coeff = linspace(.4,.6,10);
for i = 1:length(drag_coeff)
%All given variables, can be modified
g = 9.81;           %m/s2 ... acceleration due to gravity
Cd = 0.8;           %... discharge coefficient
p_air_amb = 0.961;  %kg/m^3 ... ambient air density
Volbottle = 0.002;  %m^3 ... volume of empty bottle
P_amb=12.1*6894.75; %pa ... atmospheric pressure
gamma = 1.4;        %... ratio of specific heats for air
pwater = 1000;      %kg/m^3 ... density of water
DThroat = 2.1/100;  %mm ... diameter of throat
DBottle = 10.5/100; %mm ... diameter of bottle
R = 287;            %J/kgK ... gas constant of air
P_gage = 40*6894.75;%pa ... initial gage pressure of air in bottle
V0 = 0.0;           %m/s ... initial velocity of rocket
Theta = 40*(pi/180);%radians ... initial angle of rocket
x0 = 0.0;           %m ... initial horizontal distance
z0 = 0.01;          %m ... initial vertical height
Py0 = 0;             %m ... initial sideways offset
ls= 0.5;            %m ... length of test stand

% The stuff that changes 
Tair_i = 278.9;       %K ... initial temperature of air
CD = drag_coeff(i);           %... drag coefficient
Volwater_i = 0.001; %m^3 ... initial volume of water inside bottle
mf = 0.124;         %final burnout mass of the rocket
MBottle = 0.124;     %kg ... mass of empty 2-liter bottle with cone and fins
mw = Volwater_i*pwater; % initial mass of the water


% The wind velocities in each direction
windx = 1.78816;
windy = 1;
windz = 0;


%% Calculate Derived Variables
%Initial Mass of Air
m_air_i = ((P_gage+P_amb)*(Volbottle-Volwater_i))/((R)*(Tair_i));
%initial mass of the rocket
mi = (mf + mw + m_air_i);    
%Initial Total Pressure
P0 = P_gage+P_amb;
%Area of the Throat
AThroat = ((DThroat^2)*pi)/4; %m^2 ... 
%Volume of Air
Volair_i = Volbottle-Volwater_i; %m^3 ... Initial Volume of Air
%Cross-Sectional Area of the Bottle
ABottle = (pi*(DBottle^2))/4; %m^2 ... Area of Bottle

% The ThrustInterpolation function is passed the range of Thrust values from the static test.
data = xlsread('Group09_10AM_statictest1.xlsx');
% get the thrust data during firing from the total array
Fraw = data(:,3)/0.22480894244318786;

% Select the right Thrust data
%F = Fraw(startindex:endindex);
%F = Fraw(2797:2997);

% for the TA baseline model
F = data(2944:3210,3)/0.22480894244318786;

%% Account for crosswinds and initial instantaneous deltaV
Vrelx = V0*cos(Theta)+windx;
Vrelz = V0*sin(Theta)+windz;
Vrely = windy;

%% The MassInterpolation function is passed the range of Thrust values from the static test.
%Call ODE45
tspan = [0,5];
% Initial Values
%Initial value of total mass
y0(1) = mi;
%Initial value of air mass
y0(2) = m_air_i;
%Initial value of x-direction velocity
y0(3) = Vrelx;
%Initial value of z-direction velocity
y0(4) = Vrelz;
%Initial value of y-direction velocity
y0(5) = Vrely;
%Initial value of x-direction position
y0(6) = x0;
%Initial value of z-direction position
y0(7) = z0;
%Initial value of y-direction position
y0(8) = Py0;
%Initial value of volume of air with time
y0(9) = Volair_i;

% call ode45 for the thermo model
[t,y] = ode45('ThermoModel',tspan, y0);
% Extract the mass during the launch
mthermo = y(:,1);

% adjust the time scale for the mass discharge
% For 600g launch

% For TA baseline launch
mtimes = t*1.03;


%% create an array of initial values for the Interpolation Model
K0(1) = x0;
K0(2) = Py0;
K0(3) = z0;
K0(4) = Vrelx;
K0(5) = Vrely;
K0(6) = Vrelz;

%% set the time duration of the flight
t = [0,5];
%% Call ODE45
[t,K]=ode45('InterpolationModel',t,K0);

%% separate all the data and define them
x = K(:,1); % the x distance
ypos = K(:,2)*2.1; % the y distance
z = K(:,3); % the z distance
Vx = K(:,4); % the x velocity
Vy = K(:,5); % the y velocity
Vz = K(:,6); % the z velocity

%% Plot the trajectory and show the maximum height and distance
plot3(x,ypos,z)
hold on
xlabel('x Distance [m]')
ylabel('y Distance [m]')
zlabel('Height [m]')
title('Interpolation Trajectory for TA Baseline')
fprintf('Max height: %2f\nMax distance: %2f\n',max(z),max(x))
grid on
x_el = [x_el x(end)];
y_el = [y_el ypos(end)];
end
% The TA flight went 48.829 meters. 
% The drift angle was 13 degrees.
realx = cos(deg2rad(13))*48.829;
realy = sin(deg2rad(13))*48.829;

interpx = max(x);
interpy = max(ypos);

x_shift = mean(x_el);
y_shift = mean(y_el);
std_x = std(x_el);
std_y = std(y_el);
a = std_x; % horizontal radius
b = std_y; % vertical radius
x_center = x_shift; % x0,y0 ellipse centre coordinates
y_center = y_shift;
t=-pi:0.01:pi;
x=x_center + a*cos(t);
y=y_center + b*sin(t);
plot(x,y)

std_x = std(x_el);
std_y = std(y_el);
a = 2*std_x; % horizontal radius
b = 2*std_y; % vertical radius
x_center = x_shift; % x0,y0 ellipse centre coordinates
y_center = y_shift;
t = -pi:0.01:pi;
x = x_center + a*cos(t);
y = y_center + b*sin(t);
plot(x,y)
zlim([0 25])