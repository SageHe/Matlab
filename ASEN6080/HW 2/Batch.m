clear all;close all;clc
%load in measurements from HW1
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('HW1_Nom_Measurements')
rangemeas(14319:end,:) = []; %Clean up end of measurement file that has no measurements
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
%Add noise to measuremts
rangemeasnoisy = rangemeas + rho_noise;
rrmeasnoisy = rrmeas + rhodot_noise;
%Arrange measurements in vector format for processing
% y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];


%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2 0]);
R = diag([sigma_rho sigma_rhodot sigma_rho sigma_rhodot sigma_rho sigma_rhodot]);
%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
x0 = [r;v;J2];
n = size(x0,1);
dt = t(2) - t(1);

dx = [0.8 0.6 -0.2 1e-6 1e-7 -6e-7 0]';
%initialize filter
j = 1;
xhat0(:,1) =  x0 + dx;
Deltax_minus(:,1) = zeros(7,1);
while j < 10
xhat(:,1) = xhat0(:,j);
Lambda(:,:,j) = inv(P0);
N = P0\Deltax_minus; %better alternative to inv(P0)*Deltax_minus


for i = 2:numel(t)
    %read next observation
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %integrate from t0 to ti
    Phi = eye(n);
    Phi_flat = reshape(Phi,n^2,1);
    Z = [xhat0;Phi_flat];
    tspan = [0:dt:t(i)];
    [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);

    
    
    
    
    
    
    