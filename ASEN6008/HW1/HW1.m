clear all;close all;clc
%% Problem 1
%recalc Va,Vm, need full mars mu
A_tr = (1.49598023e8 + 2.27939186e8)/2;
Va = sqrt(((2*1.32712428e11)/(2.27939186e8)) - (1.32712428e11/A_tr));
Vp = sqrt(((2*1.32712428e11)/(1.49598023e8)) - (1.32712428e11/A_tr));
Vm = sqrt(1.32712428e11/2.27939186e8);

Vinf = Va - Vm;
ahyp = (-4.305e4)/Vinf^2;
Vp_hyp = sqrt(((2*4.305e4)/(3797.2)) - (4.305e4/ahyp));
V400mars = sqrt(4.305e4/3797.2);
dv2 = V400mars - Vp_hyp;

P = 2*pi*sqrt((A_tr^3)/(1.32712428e11));
TOF = .5*P;
% Calculate impulsive dv to raise apoapsis to 12000 km
a = 7191.938817629013;
e = 0.02454974900598137;
ra_i = a*(1 + e);
rp_i = a*(1 - e);
rp_f = rp_i;
ra_f = 12000;
ef = (ra_f - rp_f)/(ra_f + rp_f);
mu = 3.98600433e5;
v_i = sqrt(mu*(2/rp_i-1/a));
a_f = ra_f/(1+ef);
v_f = sqrt(mu*(2/rp_f-1/a_f));
dv = v_f - v_i;
%% Problem 2
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
AU = 1.49597870e8; %km
a_e = 1.49598023e8; %km

Rf_m = [-578441.618274359 2.27938449869731e8 0];
Vf_m = [-24.1281802482527 -0.0612303173808154 0];
tspan = [TOF 0];
X_mars = [Rf_m Vf_m];
[t,y] = ode45(@(t,X_mars) propmars(t,X_mars),tspan,X_mars,opts);
figure
plot3(y(:,1),y(:,2),y(:,3))
grid on
grid minor
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2.5e8 2.5e8])
ylim([-2.5e8 2.5e8])
title('Mars Orbit over TOF')

R0_sc = [0 -a_e 0];
V0_sc = [Vp 0 0];
R0_e = [-578411.002878924 -1.49596751684464e8 0]; %km
V0_e = [29.7830732658560 -0.115161262358529 0]; %km/s
R0_m = y(end,1:3);
V0_m = y(end,4:6);
clear t y
tspan = [0 TOF];
% X_sc = [R0_sc V0_sc];
% Xe = [R0_e V0_e];
% Xm = [R0_m V0_m];
X = [R0_sc V0_sc R0_e V0_e R0_m V0_m R0_sc V0_sc];
[t,y] = ode45(@(t,X) twobodyode(t,X),tspan,X,opts);
figure 
% hold on
% SC = plot3(y(:,1),y(:,2),y(:,3));
% Earth = plot3(y(:,7),y(:,8),y(:,9));
% Mars = plot3(y(:,13),y(:,14),y(:,15));

plot3(y(:,1),y(:,2),y(:,3),y(:,7),y(:,8),y(:,9),y(:,13),y(:,14),y(:,15))
grid on
grid minor
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
xlim([-2.5e8 2.5e8])
ylim([-2.5e8 2.5e8])
legend('SC Orbit','Earth Orbit','Mars Orbit')
title('Spacecraft Orbit with Two Body EOMS')

figure
plot3(y(:,19),y(:,20),y(:,21),y(:,7),y(:,8),y(:,9),y(:,13),y(:,14),y(:,15))
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
xlim([-2.5e8 2.5e8])
ylim([-2.5e8 2.5e8])
grid on
grid minor
legend('SC Perturbed Orbit','Earth Orbit','Mars Orbit')
title('SC Orbit with Two Body EOMS and Perturbations')

figure
subplot(2,1,1)
hold on
plot(t,(y(:,1) - y(:,19)))
plot(t,(y(:,2) - y(:,20)))
xlabel('Time (s)')
ylabel('Position Differences (km)')
grid on
grid minor
legend('X-Position Differences','Y-Position Differences')
subplot(2,1,2)
hold on
plot(t,(y(:,4) - y(:,22)))
plot(t,(y(:,5) - y(:,23)))
xlabel('Time (s)')
ylabel('Velocity Differences (km/s)')
grid on 
grid minor
legend('X-Velocity Differences','Y-Velocity Differences')
sgtitle('Idealized Hohmann Transfer VS Perturbed S/C Transfer')