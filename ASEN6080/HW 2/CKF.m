%% Homework 2 - Orbit estimation using CKF, EKF, and Batch filtering -- CKF 
clear all;close all;clc
%load in measurements from HW1
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('HW1_Nom_Measurements')
Truth = load('HW1_p2a_traj.txt');
Truth(:,2:7) = Truth(:,2:7)/1000;
rangemeas(14319:end,:) = [];
rrmeas(14319:end,:) = [];
t(14319:end) = [];
% load('HW1_Noisy_Measurements')
J2 = 1082.63e-6;
%Calculate station positions in ECEF 
Re = 6378; %earth radius in km
theta0 = deg2rad(122);
stat1pos_ecef = [Re*cosd(-35.398333)*cosd(148.981944);Re*cosd(-35.398333)*sind(148.981944);Re*sind(-35.398333)];
stat2pos_ecef = [Re*cosd(40.427222)*cosd(355.749444);Re*cosd(40.427222)*sind(355.749444);Re*sind(40.427222)];
stat3pos_ecef = [Re*cosd(35.247164)*cosd(243.205);Re*cosd(35.247164)*sind(243.205);Re*sind(35.247164)];
statspos_ecef = [stat1pos_ecef';stat2pos_ecef';stat3pos_ecef'];
%Add Gaussian noise with specified standard deviations
sigma_rho = 1e-3;
sigma_rhodot = 1e-6;
rho_noise = sigma_rho*randn(size(rangemeas,1),1);
rhodot_noise = sigma_rhodot*randn(size(rangemeas,1),1);
sigma_r = 1;
sigma_v = 1e-3;

rangemeasnoisy = rangemeas + rho_noise;
rrmeasnoisy = rrmeas + rhodot_noise;

%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2 0]);
R = diag([sigma_rho sigma_rhodot sigma_rho sigma_rhodot sigma_rho sigma_rhodot]);
%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
x0 = [r;v;J2];
x_ti(1,:) = x0;
% dx = zeros(7,1);
% dx = [.5 .5 .5 1e-7 1e-7 1e-7 0]';
% dx = [-.8949639 -0.266738025 -0.667648229 0.0012933487 -0.00051752255 0.00015225489 0]';
dx = [0.5 -0.5 0.5 0.5e-3 -0.5e-3 0.5e-3 0]';
dt = t(2) - t(1);


% x0 = x0 + dx;
Phi = eye(size(x0,1));
Phi_flat = reshape(Phi,size(x0,1)^2,1);
Z = [x0;Phi_flat];
tspan = [0:dt:t(end)];
[~,X_true] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);

x_hat(1,:) = x0 + dx;

%Initialize Filter
ti_minus = t(1,1);
x_t_im1 = x0 + dx; %+ [37.846044185977355 29.142452212735556 -27.668907935023171 -0.012751280590239 0.029774403648732 0.014338968404368 0]';
Pi_m1 = P0;
Deltax_i_m1_plus = zeros(7,1);
numstats = size(rangemeas,2);
% y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
[t_measstop,t_measstart] = parsemeas(t,rangemeas);
gapcount = 1;
% Read next observation
i = 2;
while i <= numel(t)
    ti = t(i);
%     statnum = find(~isnan(parsedrange(i,2:4)));
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    n = size(x0,1);
    %Propagate from t_i-1 to t_i
%     dt = ti - ti_minus;
if isnan(yi)
    Phi = eye(n);
    Phi_flat = reshape(Phi,size(Phi,1)^2,1);
    Z = [x_t_im1;Phi_flat];
    tspan = [t(t_measstop(gapcount)):dt:t(t_measstart(gapcount))];
    [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    x_ti(t_measstop(gapcount):t_measstart(gapcount),:) = x(:,1:7);
    x_hat(t_measstop(gapcount):t_measstart(gapcount),:) = x(:,1:7);
    Phi = reshape(x(end,n+1:end),n,n);
    Deltaxi_minus = Phi*Deltax_i_m1_plus;
%     Deltaxi_plus = Deltaxi_minus;
    Pi_minus = Phi*Pi_m1*Phi';
%     Pi_plus = Pi_minus;
    
    i = t_measstart(gapcount);
    gapcount = gapcount + 1;
else
    if isnan(y(i-1,:))
        ti_minus = t(i-1);
    end
    Phi = eye(n);
    Phi_flat = reshape(Phi,size(Phi,1)^2,1);
    Z = [x_t_im1;Phi_flat];
    tspan = [ti_minus:dt/2:ti];
    [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    x_ti(i,:) = x(end,1:7);
    %Time Update
    Phi = reshape(x(end,n+1:end),n,n);
    Deltaxi_minus = Phi*Deltax_i_m1_plus;
    Pi_minus = Phi*Pi_m1*Phi';
    %Process Observations
    [hi(i,:),~] = predictmeas(x(end,1:6),Xs_ECI);
    ri(i,:) = yi - hi(i,:);
    
    Hi = zeros(2*size(statspos_ecef,1),n);
    for k = 1:size(statspos_ecef,1)
        if ~isnan(ri(i,2*k))
            Hi(2*k-1:2*k,1:6) = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(k,:));
        elseif isnan(ri(i,2*k))
            ri(i,2*k-1:2*k) = [0 0];
        end
%         Hi = Htilde_sc_rho_rhod(x(end,1:6),[Rs_ECI Vs_ECI]);
    end
    Ki = Pi_minus*Hi'*inv(Hi*Pi_minus*Hi' + R);
    %Measurement update
    Deltaxi_plus = Deltaxi_minus + Ki*(ri(i,:)' - Hi*Deltaxi_minus);
    Pi_plus = (eye(7) - Ki*Hi)*Pi_minus*(eye(7) - Ki*Hi)' + Ki*R*Ki';
    
    i = i+1;
    %set up for next minor iteration
end
    ti_minus = ti;
    x_t_im1 = x_ti(i-1,:)';
    Pi_mi = Pi_plus;
    Deltax_i_m1_plus = Deltaxi_plus;
    x_hat(i-1,:) = x_ti(i-1,:) + Deltaxi_plus';
    
    
end

Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [(x0 + dx);Phi_flat];
tspan = [0 t(end)];
[~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);

Phi = reshape(x(end,n+1:end),n,n);
Deltax0_plus = inv(Phi)*Deltaxi_plus;
x0 = x0 + Deltax0_plus;
state_error = x_hat - X_true(:,1:7);
figure
subplot(6,1,1)
plot(t,state_error(:,1))
% xlabel('Time (s)')
% ylabel('state error')
subplot(6,1,2)
plot(t,state_error(:,2))
% xlabel('Time (s)')
% ylabel('y state error')
subplot(6,1,3)
plot(t,state_error(:,3))
subplot(6,1,4)
plot(t,state_error(:,4))
subplot(6,1,5)
plot(t,state_error(:,5))
subplot(6,1,6)
plot(t,state_error(:,6))
sgtitle('Reference trajectory errors')

figure
subplot(6,1,1)
plot(t,x_hat(:,1))
% xlabel('Time (s)')
% ylabel('state error')
subplot(6,1,2)
plot(t,x_hat(:,2))
% xlabel('Time (s)')
% ylabel('y state error')
subplot(6,1,3)
plot(t,x_hat(:,3))
subplot(6,1,4)
plot(t,x_hat(:,4))
subplot(6,1,5)
plot(t,x_hat(:,5))
subplot(6,1,6)
plot(t,x_hat(:,6))