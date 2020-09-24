close all
clear all
clc

t(1) = 0; % [s] initial time
rb(1) = 0; % [m] initial burn grain displacement
T_claimed = csvread('Thrust.csv');

%% INPUTS
cstar_eff = ; % [-], cstar efficiency
t_step = ; % [s] time step
P_atm = ; % [Pa] ambient pressure
a = ; % [-] burn rate coefficient
n = ; % [-] burn rate exponent
cstar = ; % [m/s] characteristic velocity
h_grain = ; % [in] motor grain height
r_grain = ; % [in] motor grain radius
r_throat = ; % [in] throat radius
r_exit = ; % [in] exit radius
Mass = ; % [kg]

%% CONVERSIONS
h_grain = h_grain*0.0254; % [m] motor grain height
r_grain = r_grain*0.0254; % [m] motor grain radius
r_throat = r_throat*0.0254; % [m] throat radius
r_exit = r_exit*0.0254; % [m] exit radius

%% QUANTITY CALCULATIONS
Vol = h_grain*r_grain^2*pi(); % [m^3]
rho_p = Mass/Vol; % [kg/m^3]
A_throat = pi()*(r_throat)^2; % [m^2]
A_exit = pi()*(r_exit)^2; % [m^2]
AR_sup = A_exit/A_throat; % supersonic area ratio

V_burn = 0; % [m^3]
j = 1;
while rb < r_grain && rb < h_grain % while there is unburned grain remaining
    [A_burn(j), V_burn(j+1)] = burn_geometry(r_grain,h_grain,rb); % [m] burn area, burn cavity volume
    Pc(j) = ((a * rho_p * A_burn(j) * cstar) / (A_throat)).^((1)/(1-n))/1e6; % [MPa] chamber pressure
    burn_rate(j) = a*(Pc(j)*10^6)^n; % [m/s] burn rate
    rb = rb + burn_rate(j) * t_step; % [m] updates burn displacement
    
    % delta_Vol = ; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup, delta_Vol);
    cstar = cstar*cstar_eff; % [m/s]
    if j == 1
        t(j) = t_step;
    else
        t(j) = t(j-1) + t_step;
    end
    j = j+1;
end

