%% House Keeping
clear all;
close all;
clc

%% Define constants

% system Component: Base
Kg = 48.4;  % Total Gear Ratio
Km = 0.0107; % Nm/amp - Motor Constant
Rm = 3.29; % ohms - Armature Resitance

% system Component: Rigid Arm
j_hub = .002; % Kg*m^2 - Base Inertia
j_load = .0015; %Kg*m^2 - Load Inertia of Bar
j = j_hub + j_load; %Kg*m^2 - Total Inertia

% system Component: Flexible Link
L = .45; % m - Link Length
Marm = .06; %Kg - Link Mass of ruler
j_arm = .004; %Kg*m^2 - Link Rigid Body Inertia
Mtip = .05; %Kg - Tip Mass
j_M = .01; % Kg*m^2 - Tip Mass Inertia
fc = 1.8; % Hz - Natural Frequency
j_L = j_arm + j_M; % Kg*m^2 - 
Karm = ((2*pi*fc)^2)*(j_L); % Flexible Link Stiffness

%% Enter Equation for rigid arm
% Assign coeffcient 
Kd = -.037;
Kp = 3.23;
thetaD = .2;

% Natural Frequency
Wn = sqrt(Kp*Kg*Km/(j*Rm));

% Zeta
z = (((Kg^2) * (Km^2)) + (Kd * Kg * Km))/(2*sqrt(Kp*Kg*Km*j*Rm));

% Set up and build transfer function
num = Wn^2;
den = [ 1 (2*z*Wn) Wn^2 ];

sysTF = tf(num,den);

[x t] = step(sysTF);

x = 2*thetaD*x;

dx = diff(x)/.0003;

Vin = Kp.*(2*thetaD - x(1:end-1)) + Kd.*(-dx);

maxV = max(Vin);

figure(1); clf;
hold on
y = ones(1,length(t))*(2*thetaD);
plot(t,y);
plot(t,x);
