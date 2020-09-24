function [dydt, F] = rocketODE(t,y)
% Function is used as part of ode45 function to determine flight trajectory
% of bottle rocket. Considers multistage rocket initialized as water
% generated thrust and progressing to air generated thrust before
% terminating in the ballistic phase.
%
% Inputs: 
% t - time
% y - [velocity_x, velocity_y, distance, height, air volume, mass water, mass air]
%
% Outputs:
% dydt - [velocity_x, velocity_y, distance, height, air volume, mass water, mass air]
% F - thrust
%
% Author: 212
% Created: 12/2/17
% Modified: 12/2/17



% Initialize dydt output
dydt = zeros(7,1);
 
%% Constants 

% Fixed values
R = 287; % R [J/(kg*K)]
C_D = 0.5; % Drag coefficient
C_d = 0.8; % Discharge coefficient
rho_w = 1000; % Water density [kg/m^3]
gamma = 1.4; % Ratio of specific heats for air

% Environmental
rho_a = 0.961; % Ambient air density [kg/m^3]
P_amb = 83426.56; % Ambient pressure [Pa]
T_air_i = 300; % Initial air tempertature [K]
stand = 0.5; % Stand length [m]

% Bottle dimensions
v_bottle = 0.002; % Volume of bottle [m^3]
A_b = .25*pi*0.105^2; % Cross sectional area of bottle [m^2]
A_t = .25*pi*0.021^2; % Cross sectional area of throat [m^2]

% Launch variables
theta = 30; % Launch angle [degrees]
P_air_i = 428164.4; % Initial air pressure, gage + ambient [Pa]
v_water_i = 0.001; % Initial water volume [m^3]

% Calculable 
v_air_i = v_bottle-v_water_i; % Initial air volume [m^3]
P_end = P_air_i*(v_air_i/v_bottle)^gamma; % Pressure end of water
T_end = T_air_i*(v_air_i/v_bottle)^(gamma-1); % Temp end of water
M_air_i = (P_air_i/(R*T_air_i))*v_air_i; % Initial mass of air
M_air_2L = (P_amb*0.002)/(R*T_air_i); % Final mass of air
 
 
%% Variables
% Velocity X [m/s]
V_x = y(1);
 
% Velocity Z [m/s]
V_z = y(2);
 
% Heading
if y(3) < stand*cosd(theta) % On stand
    h_x = cosd(theta);
    h_z = sind(theta);
else                        % Off stand
    h_x = V_x/sqrt(V_x^2+V_z^2);
    h_z = V_z/sqrt(V_x^2+V_z^2);
end
 
% Volume of air
v_air = y(5);
 
% Water mass
M_w = y(6);
if M_w < 0
    M_w = 0;
end
 
% Air mass
M_air = y(7);
if M_air < 0
    M_air = 0;
end
 
% Total mass
M = 0.15 + M_w + M_air;
 
if M_w > 0 % Water phase
    P_air = ((v_air_i/v_air)^gamma)*P_air_i;
    V_e = sqrt((2*(P_air-P_amb))/rho_w);
    dydt(5) = C_d*A_t*V_e; % Change in air volume
    
elseif M_air > M_air_2L  % Air Phase
    P_end = P_air_i*(v_air_i/v_bottle)^gamma;
    T_end = T_air_i*(v_air_i/v_bottle)^(gamma-1);
    P = ((M_air/M_air_i)^gamma)*P_end;
    if P < P_amb
        P = P_amb;
    end
    rho = M_air/v_bottle;
    T = P/(rho*R);
    P_crit = P*(2/(gamma+1))^(gamma/(gamma-1));
    
    if P_crit > P_amb % Choked
        T_e = (2/(gamma+1))*T;
        P_e = P_crit;
        rho_e = P_e/(R*T_e);
        V_e = sqrt(gamma*R*T_e);
    else % Not chocked
        P_e = P_amb;
        Mach = sqrt((2/(gamma-1))*((P/P_amb)^((gamma-1)/gamma)-1));
        T_e = T/(1+((gamma+1)/2)*Mach^2);
        rho_e = P_amb/(R*T_e);
        V_e = Mach*sqrt(gamma*R*T_e);    
    end
    dydt(5) = 0; % Change in air volume
else % Ballistic Phase
    dydt(5) = 0; % Change in air volume
end
 
 
%% Acceleration Components
% Gravity
g_x = 0; % Force of gravity in x direction [m/s^2]
g_z = -9.81; % Force of gravity in z direction [m/s^2]
 
% Drag
D_x = h_x*(rho_a/2)*(V_x^2+V_z^2)*C_D*A_b;
D_z = h_z*(rho_a/2)*(V_x^2+V_z^2)*C_D*A_b;

% Thrust
if M_w > 0 % Water phase
    F = 2*C_d*(P_air-P_amb)*A_t;
    F_x = F * h_x;
    F_z = F * h_z;
elseif M_air > M_air_2L % air phase
    M_flow_air = C_d*rho_e*A_t*V_e;
    F = M_flow_air*V_e+(P_e-P_amb)*A_t;
    F_x = F * h_x;
    F_z = F * h_z;
else % Ballistic phase
    F_x = 0;
    F_z = 0;
end
 
%% Outputs
% Calculate next velocity_x (put acceleration_x in here)
dydt(1) = F_x/M - D_x/M + g_x;
 
% Calculate next velocity_z (put acceleration_z in here)
dydt(2) = F_z/M - D_z/M + g_z;
 
% Distance
dydt(3) = V_x;
 
% Height
dydt(4) = V_z;
 
 
% Mass flow
if M_w > 0 % Water phase
    dydt(6) = -C_d * A_t * sqrt(2*rho_w*(P_air-P_amb));
    dydt(7) = 0;
elseif M_air > M_air_2L % air phase
    dydt(6) = 0;
    dydt(7) = -M_flow_air;
else     % Ballistic phase
    dydt(6) = 0;
    dydt(7) = 0;
end
end
