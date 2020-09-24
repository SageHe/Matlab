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
Kp = 3.23;
Kd = -.037;
thetaD = .2;

% Natural Frequency
Wn = sqrt(Kp*Kg*Km/(j*Rm));

% Zeta
z = (((Kg^2) * (Km^2)) + (Kd * Kg * Km))/(2*sqrt(Kp*Kg*Km*j*Rm));

% Set up and build transfer function
num = Wn^2;
den = [ 1 (2*z*Wn) Wn^2 ];

sysTF = tf(num,den);

[theta t] = step(sysTF);

theta = 2*thetaD*theta;

dx = diff(theta)./diff(t);

Vin = Kp.*(2*thetaD - theta(1:end-1)) + Kd.*(-dx);

[ ~ , maxVi]= max(abs(Vin));
maxV = Vin(maxVi);

figure;
hold on
y = ones(1,length(t))*(2*thetaD);
plot(t,y);
plot(t,theta);

%% Enter equations for flexible arm
% K1 = 1.667;
% K2 = -.5556;
% K3 = 0.1667;
% K4 = 0.1667;

% Enter Equation for flexible arm
%Specify desired arm position in radians
thetaD_f = .2;
%Preallocate empty gains vector
gains = [];
%Create combinations of gains to test in theoretical model
for K1 = linspace(0,5,10) 
    for K2 = linspace(-5,0,15)
        for K3 = linspace(0,1.5,10)
            for K4 = linspace(0,1.5,10)    
p1 = - (Kg^2 * Km^2)/(j_hub*Rm);
p2 = (Kg^2 * Km^2 * L)/(j_hub*Rm);
q1 = Karm/(L*j_hub);
q2 = - (Karm*(j_hub + j_L))/(j_L*j_hub);
r1 = (Kg*Km)/(j_hub*Rm);
r2 = -(Kg*Km*L)/(j_hub*Rm);

lam3 = -p1 + K3*r1 + K4*r2;
lam2 = -q2 + K1*r1 + K2*r2 + K4*(p2*r1 - r2*p1);
lam1 = p1*q2 - q1*p2 + K3*(q1*r2 - r1*q2) + K2*(p2*r1 - r2*p1);
lam0 = K1*(q1*r2 - r1*q2);

% Position transfer function
num_f1 = [ K1*r1 0 K1*(q1*r2 - r1*q2)];
den_f1 = [ 1 lam3 lam2 lam1 lam0];
sysTF_f1 = tf(num_f1,den_f1);
[thetaF t_f1] = step(sysTF_f1);

thetaF = 2*thetaD*thetaF; %Should this be thetaD_f instead of thetaD

% figure
% plot(t_f1,thetaF)
% title('Arm')
% xlabel('Time(s)')
% ylabel('Arm Position')

% Displacement transfer function
num_f2 = [ K2*r2 K2*(p1*r2 - p2*r1) 0 ];
den_f2 = [ 1 lam3 lam2 lam1 lam0];
sysTF_f2 = tf(num_f2,den_f2);
[tipDis t_Dis] = step(sysTF_f2);

% figure
% plot(t_Dis,tipDis)
% title('Tip')
% xlabel('Time(s)')
% ylabel('Deflection (m)')
%Determine if current combination of gains satisfy specified criterion
if max(t_f1 > 1)
    A = find(t_f1 < 1);
    A = A(end);
    if ((max(thetaF) < (2*thetaD_f + .005)) && (thetaF(A) >= 2*.95*thetaD_f) && (max(abs(tipDis)) < .01))
        gains = [gains;[K1 K2 K3 K4]];
    end
end
            end
        end
    end
end