%% Homework 2 - Orbit estimation using CKF, EKF, and Batch filtering -- CKF 
clear all;close all;clc
%load in measurements from HW1
load('HW1_Nom_Measurements')
load('HW1_Noisy_Measurements')
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

% rangemeasnoisy = rangemeas + rho_noise;
% rrmeasnoisy = rrmeas + rhodot_noise;
%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2 0]);
R = diag([sigma_rho sigma_rhodot sigma_rho sigma_rhodot sigma_rho sigma_rhodot]);
%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
x0 = [r;v;J2];
dx = zeros(6,1);
% parsedrange = parsemeas(rangemeas,t);
% parsedrr = parsemeas(rrmeas,t);
% tmeas = parsedrange(:,1);

%Initialize Filter
ti_minus = t(1,1);
dt = t(2) - t(1);
x_t_im1 = x0;
Pi_m1 = P0;
Deltax_i_m1_plus = zeros(7,1);
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
numstats = size(rangemeas,2);
y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
% y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
% Read next observation
for i = 2:numel(t)
    clear hi ri
    ti = t(i);
%     statnum = find(~isnan(parsedrange(i,2:4)));
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    n = 7;
    %Propagate from t_i-1 to t_i
%     dt = ti - ti_minus;
    Phi = eye(n);
    Phi_flat = reshape(Phi,size(Phi,1)^2,1);
    Z = [x_t_im1;Phi_flat];
    tspan = [ti_minus:dt/2:ti];
    [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    x_ti = x(end,1:7);
    %Time Update
    Phi = reshape(x(end,n+1:end),n,n);
    Deltaxi_minus = Phi*Deltax_i_m1_plus;
    Pi_minus = Phi*Pi_m1*Phi';
    %Process Observations
    hi = predictmeas(x(end,1:6),Xs_ECI);
    ri = yi - hi;
    
    Hi = zeros(2*size(statspos_ecef,1),n);
    for k = 1:size(statspos_ecef,1)
        if ~isnan(ri(2*k))
            Hi(2*k-1:2*k,1:6) = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(k,:));
        elseif isnan(ri(2*k))
            ri(2*k-1:2*k) = [0 0];
        end
%         Hi = Htilde_sc_rho_rhod(x(end,1:6),[Rs_ECI Vs_ECI]);
    end
    Ki = Pi_minus*Hi'*inv(Hi*Pi_minus*Hi' + R);
    %Measurement update
%     Deltaxi_plus = Deltaxi_minus + Ki*(ri - Hi*Deltaxi_minus);
%     Pi_plus = (eye(6) - Ki*Hi)*Pi_minus*(eye(6) - Ki*Hi)' + Ki*R*Ki';
    %set up for next minor iteration
    ti_minus = ti;
    x_t_im1 = x_ti';
    
    
end