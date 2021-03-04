%% Homework 2 - Orbit estimation using CKF, EKF, and Batch filtering -- EKF 
clear all;close all;clc
tic
%load in measurements from HW1
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
% load('HW1_Nom_Measurements')
% load('HW1_Nom_Measurements_V2')
% load('HW3_test_measurements')
load('J3_Measurements')
rangemeas(14319:end,:) = []; %Clean up end of measurement file that has no measurements
rrmeas(14319:end,:) = [];
t(14319:end) = [];
truth_state(14319:end,:) = [];
% load('HW1_Noisy_Measurements')
J2 = 1082.63e-6;
mu = 3.986004415e5;
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
sigma_a = 1e-12;

rangemeasnoisy = rangemeas + rho_noise;
rrmeasnoisy = rrmeas + rhodot_noise;
y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
% y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)];

%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2 sigma_a^2 sigma_a^2 sigma_a^2]);
R = diag([1e-6 1e-12]);
%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
P = 2*pi*sqrt((10000^3)/mu);
tau = P/30;
x0 = [r;v;0;0;0];
n = size(x0,1);
dt = t(2) - t(1);

% Phi = eye(size(x0,1));
% Phi_flat = reshape(Phi,size(x0,1)^2,1);
% Z = [x0;Phi_flat];
% tspan = [0:dt:t(end)];
% [~,X_true] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

% dx = zeros(7,1);
% dx = [.5 .5 .5 1e-7 1e-7 1e-7 0]';
% dx = [-.8949639 -0.266738025 -0.667648229 0.0012933487 -0.00051752255 0.00015225489 0]';
dx = [0.5 -0.5 0.5 0.5e-7 -0.5e-7 0.5e-7 0 0 0]';
% j = 1;
% initialize filter
x_hat0(1,:) = x0 + dx;
x_hat0 = x_hat0';
x_hat(1,:) = x_hat0;
ti_m1 = t(1);
P_plus(:,:,1) = P0;
ri = NaN(numel(t),2);

mu_a = zeros(1,3);
sigma_a = 1e-7*ones(1,3);
v = mvnrnd(mu_a,sigma_a,numel(t));
% read next observation
for i = 2:numel(t)
    ti = t(i);
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    % propagate from t_i-1 to t_i
    Phi = eye(n);
    Phi_flat = reshape(Phi,n^2,1);
    Z = [x_hat(i-1,:)';Phi_flat];
    tspan = [0 10];
    [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_DMC_ODE(t,Z,tau,v(i,:)),tspan,Z,opts);
    x_hat(i,:) = x(end,1:n); % technically x_hat_minus I think
    Phi = reshape(x(end,n+1:end),n,n);
    % time update
    %calculate discrete-time P.N
    wvec = x_hat(7:9);
    Q_DT = calc_DT_PN_DMC(dt,wvec);

    P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi';% + Q_DT;
    % compute obs sensitiviteis and Kalman gain
    statnum = find(~isnan(yi));
    if ~isempty(statnum)
        statnum = statnum(2)/2;
    end
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if statnum == 1%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:6) = Htilde_sc_rho_rhod(x_hat(i,1:6),Xs_ECI(1,:));
           [hi(i,:)] = predictmeas(x_hat(i,1:6),Xs_ECI(1,:));
        elseif statnum == 2%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:6) = Htilde_sc_rho_rhod(x_hat(i,1:6),Xs_ECI(2,:));
           [hi(i,:)] = predictmeas(x_hat(i,1:6),Xs_ECI(2,:));
        elseif statnum == 3 
            Hi(:,1:6) = Htilde_sc_rho_rhod(x_hat(i,1:6),Xs_ECI(3,:));
           [hi(i,:)] = predictmeas(x_hat(i,1:6),Xs_ECI(3,:));
        end
        inds = find(~isnan(y(i,:)));
        yi = yi(1,inds);
        
        ri(i,:) = yi - hi(i,:);
        Ki = P_minus(:,:,i)*Hi'*pinv(Hi*P_minus(:,:,i)*Hi' + R);
        % measurement update
        x_hat(i,:) = x_hat(i,:) + (Ki*ri(i,:)')';
        P_plus(:,:,i) = (eye(n) - Ki*Hi)*P_minus(:,:,i)*(eye(n) - Ki*Hi)' + Ki*R*Ki';
    else
        P_plus(:,:,i) = P_minus(:,:,i);
    end
end

figure
subplot(3,1,1)
plot(t,x_hat(:,1))
subplot(3,1,2)
plot(t,x_hat(:,2))
subplot(3,1,3)
plot(t,x_hat(:,3))

figure
subplot(3,1,1)
plot(t,x_hat(:,1) - truth_state(:,1))
subplot(3,1,2)
plot(t,x_hat(:,2) - truth_state(:,2))
subplot(3,1,3)
plot(t,x_hat(:,3) - truth_state(:,3))

figure
subplot(2,1,1)
plot(t,ri(:,1),'*')
ylim([-.002 .002])
subplot(2,1,2)
plot(t,ri(:,2),'*')
ylim([-.5e-5 .5e-5])
toc