%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Charles Puskar
% ASEN 3113 Section 012
% Design Lab
%
% Created: November 29, 2018
% Modified: November 29, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% GIVENS

r_earth = 6378; % km
r_orbit = 42164; % km
e = 0.0167;
alpha_ir = 0.85;
ep_ir = 0.85;
alpha_solar = 0.2;
T_op_min = 20 + 273.15; % K
T_op_max = 30 + 273.15; % K
T_surv = -40 + 273.15; % K
G_solar = 1361; % W/m^2
G_winter = 88; % W/m^s
G_summer = 63; % W/m^s
G_equi = (G_winter + G_summer)/2; % W/m^2
G_ecl = 11; % W/m^s
i_sols = 23.5; % degrees
i_equi = 0; % degrees
instr_power = 20; % W
boltz = 5.67e-8; % W/m^2*K^4
G_adjust_winter = (1-e)^(-2);
G_adjust_summer = (1+e)^(-2);

%% MINIMUM AREA

A_equi = instr_power/(ep_ir*boltz*(T_op_max^4) - alpha_solar*G_solar*cosd(i_equi) - alpha_ir*G_equi); % m^2
A_winter = instr_power/(ep_ir*boltz*(T_op_max^4) - alpha_solar*G_solar*G_adjust_winter*cosd(i_sols) - alpha_ir*G_winter); % m^2

A_s = max(A_equi,A_winter);

%% OPERATIONAL HEATER POWER
t = linspace(0,24,24*3600+1);

% winter solstice
orbit_adjust_winter = (sin(2.*pi.*t./24) >= 0).*sin(2.*pi.*t./24);
op_power_winter = ep_ir.*boltz.*A_s.*T_op_min.^4 - alpha_solar.*A_s.*G_solar.*G_adjust_winter.*cosd(i_sols).*orbit_adjust_winter - instr_power - alpha_ir.*A_s.*G_winter;
op_power_winter(op_power_winter<0) = 0;
T_op_winter = ((alpha_solar.*A_s.*G_solar.*G_adjust_winter.*cosd(i_sols).*orbit_adjust_winter + instr_power + alpha_ir.*A_s.*G_winter)/(ep_ir.*boltz.*A_s)).^(1/4);

% summer solstice
orbit_adjust_summer = (sin(2.*pi.*t./24)>= 0).*sin(2.*pi.*t./24);
op_power_summer = ep_ir.*boltz.*A_s.*T_op_min.^4 - alpha_solar.*A_s.*G_solar.*G_adjust_summer.*cosd(i_sols).*orbit_adjust_summer - instr_power - alpha_ir.*A_s.*G_summer;
op_power_summer(op_power_summer<0) = 0;
T_op_summer = ((alpha_solar.*A_s.*G_solar.*G_adjust_summer.*cosd(i_sols).*orbit_adjust_summer + instr_power + alpha_ir.*A_s.*G_summer)/(ep_ir.*boltz.*A_s)).^(1/4);

% equinox
eclipe_half = asin(r_earth/r_orbit);
eclipse_start = pi-eclipe_half;
eclipse_end = pi+eclipe_half;
eclipse_adjust = (((2.*pi.*t/24) <= eclipse_start)+((2.*pi.*t/24) >= eclipse_end));
orbit_adjust_equi = eclipse_adjust.*(sin(2.*pi.*t./24)>= 0).*sin(2.*pi.*t./24);
op_power_equi = ep_ir.*boltz.*A_s.*T_op_min.^4 - alpha_solar.*A_s.*G_solar.*cosd(i_equi).*orbit_adjust_equi - instr_power - alpha_ir.*A_s.*G_equi.*eclipse_adjust - alpha_ir.*A_s.*G_ecl.*(~eclipse_adjust);
op_power_equi(op_power_equi<0) = 0;
T_op_equi = ((alpha_solar.*A_s.*G_solar.*cosd(i_equi).*orbit_adjust_equi + instr_power + alpha_ir.*A_s.*G_equi.*eclipse_adjust + alpha_ir.*A_s.*G_ecl.*(~eclipse_adjust))/(ep_ir.*boltz.*A_s)).^(1/4);

% plots
op_power_winter_fig = figure;
hold on;
box on;
grid minor;
yyaxis left
plot(t,op_power_winter,'linewidth',1)
ylabel('Power [W]','interpreter','latex','fontsize',14)
yyaxis right
plot(t,T_op_winter-273.15,'linewidth',1)
plot(t,T_op_min.*ones(size(t))-273.15)
ylabel('Temperature [C$^\circ$]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Operational Heater Power (Winter Solstice)','interpreter','latex','fontsize',16)
legend({'Operational Power','Unheated Temperature','Minimum Operational Temperature'},'interpreter','latex','location','southeast');
saveas(op_power_winter_fig,'op_power_winter.jpg');

op_power_summer_fig = figure;
hold on;
box on;
grid minor;
yyaxis left
plot(t,op_power_summer,'linewidth',1)
ylabel('Power [W]','interpreter','latex','fontsize',14)
yyaxis right
plot(t,T_op_summer-273.15,'linewidth',1)
plot(t,T_op_min.*ones(size(t))-273.15)
ylabel('Temperature [C$^\circ$]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Operational Heater Power (Summer Solstice)','interpreter','latex','fontsize',16)
legend({'Operational Power','Unheated Temperature','Minimum Operational Temperature'},'interpreter','latex','location','southeast');
saveas(op_power_summer_fig,'op_power_summer.jpg');

op_power_equi_fig = figure;
hold on;
box on;
grid minor;
yyaxis left
plot(t,op_power_equi,'linewidth',1)
ylabel('Power [W]','interpreter','latex','fontsize',14)
yyaxis right
plot(t,T_op_equi-273.15,'linewidth',1)
plot(t,T_op_min.*ones(size(t))-273.15)
ylabel('Temperature [C$^\circ$]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Operational Heater Power (Equinox)','interpreter','latex','fontsize',16)
legend({'Operational Power','Unheated Temperature','Minimum Operational Temperature'},'interpreter','latex','location','southeast');
saveas(op_power_equi_fig,'op_power_equi.jpg');

%% SURVIVAL HEATER POWER

% winter solstice
surv_power_winter = ep_ir.*boltz.*A_s.*T_surv.^4 - alpha_solar.*A_s.*G_solar.*G_adjust_winter.*cosd(i_sols).*orbit_adjust_winter - alpha_ir.*A_s.*G_winter;
surv_power_winter(surv_power_winter<0) = 0;
T_surv_winter = ((alpha_solar.*A_s.*G_solar.*G_adjust_winter.*cosd(i_sols).*orbit_adjust_winter + alpha_ir.*A_s.*G_winter)/(ep_ir.*boltz.*A_s)).^(1/4);

% summer solstice
surv_power_summer = ep_ir.*boltz.*A_s.*T_surv.^4 - alpha_solar.*A_s.*G_solar.*G_adjust_summer.*cosd(i_sols).*orbit_adjust_summer - alpha_ir.*A_s.*G_summer;
surv_power_summer(surv_power_summer<0) = 0;
T_surv_summer = ((alpha_solar.*A_s.*G_solar.*G_adjust_summer.*cosd(i_sols).*orbit_adjust_summer + alpha_ir.*A_s.*G_summer)/(ep_ir.*boltz.*A_s)).^(1/4);

% equinox
surv_power_equi = ep_ir.*boltz.*A_s.*T_surv.^4 - alpha_solar.*A_s.*G_solar.*cosd(i_equi).*orbit_adjust_equi - alpha_ir.*A_s.*G_equi.*eclipse_adjust - alpha_ir.*A_s.*G_ecl.*(~eclipse_adjust);
surv_power_equi(surv_power_equi<0) = 0;
T_surv_equi = ((alpha_solar.*A_s.*G_solar.*cosd(i_equi).*orbit_adjust_equi + alpha_ir.*A_s.*G_equi.*eclipse_adjust + alpha_ir.*A_s.*G_ecl.*(~eclipse_adjust))/(ep_ir.*boltz.*A_s)).^(1/4);

% plots
surv_power_winter_fig = figure;
hold on;
box on;
grid minor;
yyaxis left
plot(t,surv_power_winter,'linewidth',1)
ylabel('Power [W]','interpreter','latex','fontsize',14)
yyaxis right
plot(t,T_surv_winter-273.15,'linewidth',1)
plot(t,T_surv.*ones(size(t))-273.15)
ylabel('Temperature [C$^\circ$]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Survival Heater Power (Winter Solstice)','interpreter','latex','fontsize',16)
legend({'Survival Power','Unheated Temperature','Minimum Survival Temperature'},'interpreter','latex','location','southeast');
saveas(surv_power_winter_fig,'surv_power_winter.jpg');

surv_power_summer_fig = figure;
hold on;
box on;
grid minor;
yyaxis left
plot(t,surv_power_summer,'linewidth',1)
ylabel('Power [W]','interpreter','latex','fontsize',14)
yyaxis right
plot(t,T_surv_summer-273.15,'linewidth',1)
plot(t,T_surv.*ones(size(t))-273.15)
ylabel('Temperature [C$^\circ$]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Survival Heater Power (Summer Solstice)','interpreter','latex','fontsize',16)
legend({'Survival Power','Unheated Temperature','Minimum Survival Temperature'},'interpreter','latex','location','southeast');
saveas(surv_power_summer_fig,'surv_power_summer.jpg');

surv_power_equi_fig = figure;
hold on;
box on;
grid minor;
yyaxis left
plot(t,surv_power_equi,'linewidth',1)
ylabel('Power [W]','interpreter','latex','fontsize',14)
yyaxis right
plot(t,T_surv_equi-273.15,'linewidth',1)
plot(t,T_surv.*ones(size(t))-273.15)
ylabel('Temperature [C$^\circ$]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Survival Heater Power (Equinox)','interpreter','latex','fontsize',16)
legend({'Survival Power','Unheated Temperature','Minimum Survival Temperature'},'interpreter','latex','location','southeast');
saveas(surv_power_equi_fig,'surv_power_equi.jpg');

%% SOLAR RADIATION FLUX

solar_winter = alpha_solar.*A_s.*G_solar.*G_adjust_winter.*cosd(i_sols).*orbit_adjust_winter;
solar_summer = alpha_solar.*A_s.*G_solar.*G_adjust_summer.*cosd(i_sols).*orbit_adjust_summer;
solar_equi = alpha_solar.*A_s.*G_solar.*cosd(i_equi).*orbit_adjust_equi;

solar_fig = figure;
hold on;
box on;
grid minor;
plot(t,solar_winter,'linewidth',1)
plot(t,solar_summer,'linewidth',1)
plot(t,solar_equi,'linewidth',1)
ylabel('Solar Radiation Flux [W]','interpreter','latex','fontsize',14)
xlabel('Time [hrs]','interpreter','latex','fontsize',14);
title('Solar Radiation Flux at Solstices and Equinox','interpreter','latex','fontsize',16)
legend({'Winter','Summer','Equinox'},'interpreter','latex','location','northeast');
saveas(surv_power_equi_fig,'solar_flux.jpg');
