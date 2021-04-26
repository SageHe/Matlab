clear all;close all;clc

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('Project2_Prob2_truth_traj_50days.mat')

tspan = Tt_50;
state = Xt_50(1,1:7);
n = 7;

Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);

Z = [state';Phi_flat];

[time,X] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);

figure
plot3(X(:,1),X(:,2),X(:,3))
grid on
grid minor
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km')

figure
subplot(3,1,1)
plot(time,(X(:,1) - Xt_50(:,1)))
xlabel('Time (s)')
ylabel('X Pos. error (km)')
grid on
grid minor
subplot(3,1,2)
plot(time,(X(:,2) - Xt_50(:,2)))
xlabel('Time (s)')
ylabel('Y Pos. error (km)')
grid on
grid minor
subplot(3,1,3)
plot(time,(X(:,3) - Xt_50(:,3)))
xlabel('Time (s)')
ylabel('Z Pos. error (km)')
grid on
grid minor

figure
subplot(3,1,1)
plot(time,(X(:,4) - Xt_50(:,4)))
xlabel('Time (s)')
ylabel('X Vel. error (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(time,(X(:,5) - Xt_50(:,5)))
xlabel('Time (s)')
ylabel('Y Vel. error (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(time,(X(:,6) - Xt_50(:,6)))
xlabel('Time (s)')
ylabel('Z Vel. error (km/s)')
grid on
grid minor
%% 
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

load('2a_meas')
load('Project2_Prob2_truth_traj_50days.mat')

y = refmat_data(:,2:7);

Re = 6378.1363; %earth radius in km
theta0 = deg2rad(0);
h1 = 0.691750; h2 = 0.834539; h3 = 1.07114904;

stat1pos_ecef = [(Re+h1)*cosd(-35.398333)*cosd(148.981944);(Re+h1)*cosd(-35.398333)*sind(148.981944);(Re+h1)*sind(-35.398333)]; 
stat2pos_ecef = [(Re+h2)*cosd(40.427222)*cosd(-355.749444);(Re+h2)*cosd(40.427222)*sind(-355.749444);(Re+h2)*sind(40.427222)];
stat3pos_ecef = [(Re+h3)*cosd(35.247164)*cosd(243.205);(Re+h3)*cosd(35.247164)*sind(243.205);(Re+h3)*sind(35.247164)];
statspos_ecef = [stat1pos_ecef';stat2pos_ecef';stat3pos_ecef'];

tspan = refmat_data(:,1);  %playing with Cd value showed Cd = 1.000045 resulted in Gaussian residuals over entire 200 day timespan
state = Xt_50(1,1:7);
% state(7) = 1.000045;
n = 7;

Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);

Z = [state';Phi_flat];

[time,X] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);

for i = 1:numel(time)
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,time(i),theta0);
    
    statnum = find(~isnan(yi));
    statnum = statnum(2)/2;
    
    if statnum == 1
        [hi(i,:)] = predictmeas(X(i,1:6),Xs_ECI(1,:));
    elseif statnum == 2
        [hi(i,:)] = predictmeas(X(i,1:6),Xs_ECI(2,:));
    elseif statnum == 3
        [hi(i,:)] = predictmeas(X(i,1:6),Xs_ECI(3,:));
    end
    inds = find(~isnan(yi));
    ri(i,:) = hi(i,:) - yi(inds);
end

figure
subplot(2,1,1)
plot((time/3600/24),ri(:,1),'*')
grid on
grid minor
xlabel('Time (days)')
ylabel('Range Res. (km)')
subplot(2,1,2)
plot((time/3600/24),ri(:,2),'*')
grid on
grid minor
xlabel('Time (days)')
ylabel('Range Rate Res. (km/s)')