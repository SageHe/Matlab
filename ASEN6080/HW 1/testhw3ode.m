%% Problem 4
clear all; 

mu = 3.986004415e5;
T = 2*pi*sqrt((10000^3)/mu);
tspan = [0:10:15*T];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[R,V] = calcposvel(10000,0.001,40,80,40,0);
J2 = 1082.63e-6;
Phi = eye(6);
Phi_flat = reshape(Phi,36,1);
Z = [R;V;Phi_flat];
[t,y] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

%Convert ECEF coords. to ECI using basic DCM
We = 7.2921158553e-5; %Earth rotation rate, rads/s
phi0 = deg2rad(122);
for i = 1:numel(t)
    phi = phi0 + We*t(i);
    C_ECI(:,:,i) = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1];
    ECEFpos(:,i) = C_ECI(:,:,i)'*y(i,1:3)';
    ECEFvel(:,i) = C_ECI(:,:,i)'*y(i,4:6)';
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
statspos_ecef = [stat1pos_ecef';stat2pos_ecef';stat3pos_ecef'];
%calculate ECEF to ENU transformation matrices for three stations
stat_coords = [-35.398333 148.981944;40.427222 355.749444;35.247164 243.205];
% C_ecef2enu_stat1 = [-sind(stat_coords(1,2)) cosd(stat_coords(1,2)) 0;...
%                     -sind(stat_coords(1,1))*cosd(stat_coords(1,2)) -sind(stat_coords(1,1))*sind(stat_coords(1,2)) cosd(stat_coords(1,1));...
%                     cosd(stat_coords(1,1))*cosd(stat_coords(1,2)) cosd(stat_coords(1,1))*sind(stat_coords(1,2)) sind(stat_coords(1,1))];
%                 
% C_ecef2enu_stat2 = [-sind(stat_coords(2,2)) cosd(stat_coords(2,2)) 0;...
%                     -sind(stat_coords(2,1))*cosd(stat_coords(2,2)) -sind(stat_coords(2,1))*sind(stat_coords(2,2)) cosd(stat_coords(2,1));...
%                     cosd(stat_coords(2,1))*cosd(stat_coords(2,2)) cosd(stat_coords(2,1))*sind(stat_coords(2,2)) sind(stat_coords(2,1))];    
%  
% C_ecef2enu_stat3 = [-sind(stat_coords(3,2)) cosd(stat_coords(3,2)) 0;...
%                     -sind(stat_coords(3,1))*cosd(stat_coords(3,2)) -sind(stat_coords(3,1))*sind(stat_coords(3,2)) cosd(stat_coords(3,1));...
%                     cosd(stat_coords(3,1))*cosd(stat_coords(3,2)) cosd(stat_coords(3,1))*sind(stat_coords(3,2)) sind(stat_coords(3,1))]; 
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
%     rho(i,:,1) = C_ecef2enu_stat1*(ECEFpos(i,:)' - stat1pos_ecef);
%     rho(i,:,2) = C_ecef2enu_stat2*(ECEFpos(i,:)' - stat2pos_ecef);
%     rho(i,:,3) = C_ecef2enu_stat3*(ECEFpos(i,:)' - stat3pos_ecef);
%     ENUvel(i,:,1) = C_ecef2enu_stat1*ECEFvel(i,:)';
%     ENUvel(i,:,2) = C_ecef2enu_stat2*ECEFvel(i,:)';
%     ENUvel(i,:,3) = C_ecef2enu_stat3*ECEFvel(i,:)';
%     El(i,1) = asind((rho(i,3,1))/norm(rho(i,:,1)));
%     El(i,2) = asind((rho(i,3,2))/norm(rho(i,:,2)));
%     El(i,3) = asind((rho(i,3,3))/norm(rho(i,:,3)));
%     range1(i) = norm(rho(i,:,1));
%     range2(i) = norm(rho(i,:,2));
%     range3(i) = norm(rho(i,:,3));
%     rr1(i) = (dot(rho(i,:,1),ENUvel(i,:,1))/norm(rho(i,:,1)));
%     rr2(i) = (dot(rho(i,:,2),ENUvel(i,:,2))/norm(rho(i,:,2)));
%     rr3(i) = (dot(rho(i,:,3),ENUvel(i,:,3))/norm(rho(i,:,3)));
%     rho(i,:,1) = y(i,1:3)' - C_ECI(:,:,i)*stat1pos_ecef;
%     rho(i,:,2) = y(i,1:3)' - C_ECI(:,:,i)*stat2pos_ecef;
%     rho(i,:,3) = y(i,1:3)' - C_ECI(:,:,i)*stat3pos_ecef;
%     ECIvel(i,:,1) = C_ECI(:,:,i)*(cross([0 0 We],[stat1pos_ecef]))';
%     ECIvel(i,:,2) = C_ECI(:,:,i)*(cross([0 0 We],[stat2pos_ecef]))';
%     ECIvel(i,:,3) = C_ECI(:,:,i)*(cross([0 0 We],[stat3pos_ecef]))';
%     El(i,1) = asind((rho(i,3,1))/norm(rho(i,:,1)));
%     El(i,2) = asind((rho(i,3,2))/norm(rho(i,:,2)));
%     El(i,3) = asind((rho(i,3,3))/norm(rho(i,:,3)));
%     range1(i) = norm(rho(i,:,1));
%     range2(i) = norm(rho(i,:,2));
%     range3(i) = norm(rho(i,:,3));
%     rr1(i) = (dot(rho(i,:,1),(y(i,4:6) - ECIvel(i,:,1)))/norm(rho(i,:,1)));
%     rr2(i) = (dot(rho(i,:,2),(y(i,4:6) - ECIvel(i,:,2)))/norm(rho(i,:,2)));
%     rr3(i) = (dot(rho(i,:,3),(y(i,4:6) - ECIvel(i,:,3)))/norm(rho(i,:,3)));
%     
%     RU1(i) = (221/749)*(range1(i)/C)*f_t_ref;
%     RU2(i) = (221/749)*(range2(i)/C)*f_t_ref;
%     RU3(i) = (221/749)*(range3(i)/C)*f_t_ref;
%     fshift1(i) = ((-2*rr1(i))/C)*f_t_ref;
%     fshift2(i) = ((-2*rr2(i))/C)*f_t_ref;
%     fshift3(i) = ((-2*rr3(i))/C)*f_t_ref;
    Xs_ECI = stats_ECI(statspos_ecef,t(i),phi0);
    
    [Y(i,:),El(i,:)] = predictmeas(y(i,1:6),Xs_ECI); %here y is the propagated sc state
%     El(i,1) = asind
% 
%     if El(i,1) <= 10
% %         plot(t(i),El(i,1),'b*')
%        El(i,1) = NaN; 
%        range1(i) = NaN;
%        rr1(i) = NaN;
%        RU1(i) = NaN;
%        fshift1(i) = NaN;
%     end
%     if El(i,2) <= 10
% %         plot(t(i),El(i,2),'g*')
%        El(i,2) = NaN; 
%        range2(i) = NaN;
%        rr2(i) = NaN;
%        RU2(i) = NaN;
%        fshift2(i) = NaN;
%     end
%     if El(i,3) <= 10
% %         plot(t(i),El(i,3),'r*')
%        El(i,3) = NaN; 
%        range3(i) = NaN; 
%        rr3(i) = NaN;
%        RU3(i) = NaN;
%        fshift3(i) = NaN;
%     end
%     if El(i,3) > 10
%         tmeas = [tmeas t(i)];
%     end
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
plot(t,Y(:,1),'b')
plot(t,Y(:,3),'g')
plot(t,Y(:,5),'r')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range \rho, [km]')
subplot(3,1,3)
hold on
plot(t,Y(:,2),'b')
plot(t,Y(:,4),'g')
plot(t,Y(:,6),'r')
xlabel('Time (s)')
ylabel('Range Rate $\dot{\rho}$, [km/s]','Interpreter','latex')
grid on
grid minor
sgtitle('Ideal Measurements for 15 Orbits')
% % Investigate what data looks like in actual units that come from DSN
% figure
% subplot(2,1,1)
% hold on
% plot(t,RU1)
% plot(t,RU2)
% plot(t,RU3)
% xlabel('Time (s)')
% ylabel('Range Units, RU')
% grid on 
% grid minor
% legend('Station 1','Station 2','Station 3')
% subplot(2,1,2)
% hold on
% plot(t,fshift1)
% plot(t,fshift2)
% plot(t,fshift3)
% xlabel('Time (s)')
% ylabel('Doppler Shift f_{shift}, [Hz]')
% grid on
% grid minor
% sgtitle('DSN Measurements for 15 Orbits')

% Add Gaussian noise with standard dev. of 0.5 mm/s to range-rate plot
noise = (0.5e-6*randn(1,size(Y,1)));
noise = noise';
sn1 = Y(:,2) + noise;
sn2 = Y(:,4) + noise;
sn3 = Y(:,6) + noise;
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
plot(t,sn1 - Y(:,2))
plot(t,sn2 - Y(:,4))
plot(t,sn3 - Y(:,6))
grid on
grid minor
xlabel('Time (s)')
ylabel('$\dot{\rho}_{noisy} - \dot{\rho}_{ideal}$ [km/s]','Interpreter','latex')
title('Range Rate Measurement Difference, \sigma = 0.5 mm/s')