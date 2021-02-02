clear all;close all;clc
%% Problem 2
Truth = load('HW1_p2a_traj.txt');
Truth(:,2:7) = Truth(:,2:7)/1000;

mu = 3.986004415e5;
T = 2*pi*sqrt((10000^3)/mu);
tspan = [0:60:86400];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[R,V] = calcposvel(10000,0.001,40,80,40,0);
J2 = 1082.63e-6;
Phi = eye(7);
Phi_flat = reshape(Phi,49,1);
Z = [R;V;J2;Phi_flat];
[t,y] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);

figure
plot3(y(:,1),y(:,2),y(:,3))
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Satellite Orbit Reference Trajectory')

figure
plot3(Truth(:,2),Truth(:,3),Truth(:,4))
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Satellite Orbit Truth Trajectory')
% Truth data goes to 86400 seconds, equivalent to index 27236 of ref. traj.
figure
subplot(2,1,1)
hold on
plot(t,(y(:,1) - Truth(:,2)))
plot(t,y(:,2) - Truth(:,3))
plot(t,y(:,3) - Truth(:,4))
grid on
grid minor
xlabel('Time (s)')
ylabel('r_{ref.} - r_{truth} [km]')
legend('r1','r2','r3')
subplot(2,1,2)
hold on
plot(t,(y(:,4) - Truth(:,5)))
plot(t,(y(:,5) - Truth(:,6)))
plot(t,(y(:,6) - Truth(:,7)))
grid on
grid minor
xlabel('Time(s)')
ylabel('v_{ref} - v_{truth} [km/2]')
legend('v1','v2','v3')
sgtitle('Truth VS Reference States Over Time')

clear t y

deltaX = [1 0 0 0 0.01 0 0]';

mu = 3.986004415e5;
T = 2*pi*sqrt((10000^3)/mu);
tspan = [0 15*T];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[R,V] = calcposvel(10000,0.001,40,80,40,0);
J2 = 1082.63e-6;
Phi = eye(7);
Phi_flat = reshape(Phi,49,1);
Z = [R;V;J2;Phi_flat];
[t,y] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);

R = R + [1;0;0];
V = V + [0;0.01;0];
Z = [R;V;J2;Phi_flat];
[tpert,ypert] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),t,Z,opts);
dx_t = zeros(7,numel(t));
for i = 1:numel(t)
    Phi = reshape(y(i,8:end),7,7);
    dx_t(:,i) = Phi*deltaX;
end
dx_t = dx_t';
figure
subplot(2,1,1)
hold on
plot(t,(y(:,1) - ypert(:,1)))
plot(t,(y(:,2) - ypert(:,2)))
plot(t,(y(:,3) - ypert(:,3)))
plot(t,dx_t(:,1))
plot(t,dx_t(:,2))
plot(t,dx_t(:,3))
xlabel('Time (s)')
ylabel('\delta r [km]')
grid on
grid minor
legend('\delta r1','\delta r2','\delta r3','\delta r1_{\Phi}','\delta r2_{\Phi}','\delta r3_{\Phi}')
subplot(2,1,2)
hold on
plot(t,(y(:,4) - ypert(:,4)))
plot(t,(y(:,5) - ypert(:,5)))
plot(t,(y(:,6) - ypert(:,6)))
plot(t,dx_t(:,4))
plot(t,dx_t(:,5))
plot(t,dx_t(:,6))
xlabel('Time (s)')
ylabel('\delta v [km/s]')
grid on
grid minor
legend('\delta v1','\delta v2','\delta v3','\delta v1_{\Phi}','\delta v2_{\Phi}','\delta v3_{\Phi}')
sgtitle('Deviation Vectors VS Time')

for i = 1:size(dx_t,1)
    dx_t_norm_pos(i) = norm(y(i,1:3) - ypert(i,1:3));
    dx_t_norm_vel(i) = norm(y(i,4:6) - ypert(i,4:6));
    dx_t_norm_pos_stm(i) = norm(dx_t(i,1:3));
    dx_t_norm_vel_stm(i) = norm(dx_t(i,4:6));
end
figure
subplot(2,1,1)
hold on
plot(t,dx_t_norm_pos)
plot(t,dx_t_norm_pos_stm)
grid on
grid minor
xlabel('Time (s)')
ylabel('||\delta r|| [km]')
legend('||r_{ref} - r_{pert}||','||\delta r_{\Phi}||')
subplot(2,1,2)
hold on
plot(t,dx_t_norm_vel)
plot(t,dx_t_norm_vel_stm)
grid on
grid minor
xlabel('Time (s)')
ylabel('||\delta v|| [km/s]')
legend('||v_{ref} - v_{pert}||','||\delta v_{\Phi}||')
sgtitle('Deviation Vectors VS Time')

pertdiff = ypert(:,1:6) - y(:,1:6);
dev_diff = pertdiff(:,1:6) - dx_t(:,1:6);
pos_normvec = vecnorm(dev_diff(:,1:3)');
vel_normvec = vecnorm(dev_diff(:,4:6)');

% figure
% subplot(2,1,1)
% hold on
% plot(t,dev_diff(:,1))
% plot(t,dev_diff(:,2))
% plot(t,dev_diff(:,3))
% grid on
% grid minor
% xlabel('Time (s)')
% ylabel('\delta r_{int} - \delta r_{\Phi} [km]')
% legend('\delta r1','\delta r2','\delta r3')
% subplot(2,1,2)
% hold on
% plot(t,dev_diff(:,4))
% plot(t,dev_diff(:,5))
% plot(t,dev_diff(:,6))
% grid on
% grid minor
% xlabel('Time (s)')
% ylabel('\delta v_{int} - \delta v_{\Phi} [km/s]')
% legend('\delta v1','\delta v2','\delta v3')
% sgtitle('Deviation Difference VS Time')

figure
subplot(2,1,1)
plot(t,pos_normvec)
xlabel('Time (s)')
ylabel('||\delta r_{int} - \delta r_{\Phi}|| [km]')
grid on
grid minor
subplot(2,1,2)
plot(t,vel_normvec)
xlabel('Time (s)')
ylabel('||\delta v_{int} - \delta v_{\Phi}|| [km/s]')
grid on
grid minor
sgtitle('Deviation Vector Difference VS Time')
%% Problem 3

%% Problem 4
clear all; 

mu = 3.986004415e5;
T = 2*pi*sqrt((10000^3)/mu);
tspan = [0:10:15*T];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[R,V] = calcposvel(10000,0.001,40,80,40,0);
J2 = 1082.63e-6;
Phi = eye(7);
Phi_flat = reshape(Phi,49,1);
Z = [R;V;J2;Phi_flat];
[t,y] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);

%Convert ECI satellite coords. to ECEF using basic DCM
We = 7.2921158553e-5; %Earth rotation rate, rads/s
phi0 = deg2rad(122);
for i = 1:numel(t)
    phi = phi0 + We*t(i);
    C = [cos(-phi) -sin(-phi) 0;sin(-phi) cos(-phi) 0;0 0 1];
    ECEFpos(:,i) = C*y(i,1:3)';
    ECEFvel(:,i) = C*y(i,4:6)';
end
ECEFpos = ECEFpos';
ECEFvel = ECEFvel';
figure
plot3(ECEFpos(:,1),ECEFpos(:,2),ECEFpos(:,3));
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('ECEF Trajectory Over 15 Orbits')

%Calculate station positions in ECEF 
Re = 6378; %earth radius in km
stat1pos_ecef = [Re*cosd(-35.398333)*cosd(148.981944);Re*cosd(-35.398333)*sind(148.981944);Re*sind(-35.398333)];
stat2pos_ecef = [Re*cosd(40.427222)*cosd(355.749444);Re*cosd(40.427222)*sind(355.749444);Re*sind(40.427222)];
stat3pos_ecef = [Re*cosd(35.247164)*cosd(243.205);Re*cosd(35.247164)*sind(243.205);Re*sind(35.247164)];

%calculate ECEF to ENU transformation matrices for three stations
stat_coords = [-35.398333 148.981944;40.427222 355.749444;35.247164 243.205];
C_ecef2enu_stat1 = [-sind(stat_coords(1,2)) cosd(stat_coords(1,2)) 0;...
                    -sind(stat_coords(1,1))*cosd(stat_coords(1,2)) -sind(stat_coords(1,1))*sind(stat_coords(1,2)) cosd(stat_coords(1,1));...
                    cosd(stat_coords(1,1))*cosd(stat_coords(1,2)) cosd(stat_coords(1,1))*sind(stat_coords(1,2)) sind(stat_coords(1,1))];
                
C_ecef2enu_stat2 = [-sind(stat_coords(2,2)) cosd(stat_coords(2,2)) 0;...
                    -sind(stat_coords(2,1))*cosd(stat_coords(2,2)) -sind(stat_coords(2,1))*sind(stat_coords(2,2)) cosd(stat_coords(2,1));...
                    cosd(stat_coords(2,1))*cosd(stat_coords(2,2)) cosd(stat_coords(2,1))*sind(stat_coords(2,2)) sind(stat_coords(2,1))];    
 
C_ecef2enu_stat3 = [-sind(stat_coords(3,2)) cosd(stat_coords(3,2)) 0;...
                    -sind(stat_coords(3,1))*cosd(stat_coords(3,2)) -sind(stat_coords(3,1))*sind(stat_coords(3,2)) cosd(stat_coords(3,1));...
                    cosd(stat_coords(3,1))*cosd(stat_coords(3,2)) cosd(stat_coords(3,1))*sind(stat_coords(3,2)) sind(stat_coords(3,1))]; 
rho = zeros(numel(t),3,3);    
El = zeros(numel(t),3);
f_t_ref = 8.44e9;
C = 2.99792458e5; %m/s
tmeas = [];
figure
subplot(3,1,1)
hold on
grid on
grid minor
for i = 1:numel(t)
    rho(i,:,1) = C_ecef2enu_stat1*(ECEFpos(i,:)' - stat1pos_ecef);
    rho(i,:,2) = C_ecef2enu_stat2*(ECEFpos(i,:)' - stat2pos_ecef);
    rho(i,:,3) = C_ecef2enu_stat3*(ECEFpos(i,:)' - stat3pos_ecef);
    ENUvel(i,:,1) = C_ecef2enu_stat1*ECEFvel(i,:)';
    ENUvel(i,:,2) = C_ecef2enu_stat2*ECEFvel(i,:)';
    ENUvel(i,:,3) = C_ecef2enu_stat3*ECEFvel(i,:)';
    El(i,1) = asind((rho(i,3,1))/norm(rho(i,:,1)));
    El(i,2) = asind((rho(i,3,2))/norm(rho(i,:,2)));
    El(i,3) = asind((rho(i,3,3))/norm(rho(i,:,3)));
    range1(i) = norm(rho(i,:,1));
    range2(i) = norm(rho(i,:,2));
    range3(i) = norm(rho(i,:,3));
    rr1(i) = (dot(rho(i,:,1),ENUvel(i,:,1))/norm(rho(i,:,1)));
    rr2(i) = (dot(rho(i,:,2),ENUvel(i,:,2))/norm(rho(i,:,2)));
    rr3(i) = (dot(rho(i,:,3),ENUvel(i,:,3))/norm(rho(i,:,3)));
    
    RU1(i) = (221/749)*(range1(i)/C)*f_t_ref;
    RU2(i) = (221/749)*(range2(i)/C)*f_t_ref;
    RU3(i) = (221/749)*(range3(i)/C)*f_t_ref;
    fshift1(i) = ((-2*rr1(i))/C)*f_t_ref;
    fshift2(i) = ((-2*rr2(i))/C)*f_t_ref;
    fshift3(i) = ((-2*rr3(i))/C)*f_t_ref;
    if El(i,1) <= 10
%         plot(t(i),El(i,1),'b*')
       El(i,1) = NaN; 
       range1(i) = NaN;
       rr1(i) = NaN;
       RU1(i) = NaN;
       fshift1(i) = NaN;
    end
    if El(i,2) <= 10
%         plot(t(i),El(i,2),'g*')
       El(i,2) = NaN; 
       range2(i) = NaN;
       rr2(i) = NaN;
       RU2(i) = NaN;
       fshift2(i) = NaN;
    end
    if El(i,3) <= 10
%         plot(t(i),El(i,3),'r*')
       El(i,3) = NaN; 
       range3(i) = NaN; 
       rr3(i) = NaN;
       RU3(i) = NaN;
       fshift3(i) = NaN;
    end
    if El(i,3) > 10
        tmeas = [tmeas t(i)];
    end
end
% vis = El > 10;
% El = El(vis);
plot(t,El(:,1),'b')
plot(t,El(:,2),'g')
plot(t,El(:,3),'r')
xlabel('Time (s)')
ylabel('Elevation Angle (degrees)')
legend('Station 1','Station 2','Station 3')
subplot(3,1,2)
hold on
plot(t,range1,'b')
plot(t,range2,'g')
plot(t,range3,'r')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range \rho, [km]')
subplot(3,1,3)
hold on
plot(t,rr1,'b')
plot(t,rr2,'g')
plot(t,rr3,'r')
xlabel('Time (s)')
ylabel('Range Rate $\dot{\rho}$, [km/s]','Interpreter','latex')
grid on
grid minor
sgtitle('Ideal Measurements for 15 Orbits')
% Investigate what data looks like in actual units that come from DSN
figure
subplot(2,1,1)
hold on
plot(t,RU1)
plot(t,RU2)
plot(t,RU3)
xlabel('Time (s)')
ylabel('Range Units, RU')
grid on 
grid minor
legend('Station 1','Station 2','Station 3')
subplot(2,1,2)
hold on
plot(t,fshift1)
plot(t,fshift2)
plot(t,fshift3)
xlabel('Time (s)')
ylabel('Doppler Shift f_{shift}, [Hz]')
grid on
grid minor
sgtitle('DSN Measurements for 15 Orbits')
% Add Gaussian noise with standard dev. of 0.5 mm/s to range-rate plot
noise = (0.5e-6*randn(1,numel(rr1)));
sn1 = rr1 + noise;
sn2 = rr2 + noise;
sn3 = rr3 + noise;
figure
subplot(2,1,1)
hold on
plot(t,sn1)
plot(t,sn2)
plot(t,sn3)
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Rate $\dot{\rho}$, [km/s]','Interpreter','latex')
title('Range Rate with Added Noise, \sigma = 0.5 mm/s')
legend('Station 1','Station 2','Station 3')
subplot(2,1,2)
hold on
plot(t,sn1 - rr1)
plot(t,sn2 - rr2)
plot(t,sn3 - rr3)
grid on
grid minor
xlabel('Time (s)')
ylabel('$\dot{\rho}_{noisy} - \dot{\rho}_{ideal}$ [km/s]','Interpreter','latex')
title('Range Rate Measurement Difference, \sigma = 0.5 mm/s')