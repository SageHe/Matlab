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
R = diag([1e-3 1e-6]);%diag([sigma_rho sigma_rhodot sigma_rho sigma_rhodot sigma_rho sigma_rhodot]);
%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
x0 = [r;v;J2];
n = size(x0,1);
dt = t(2) - t(1);

dx = [0.8 0.6 -0.2 1e-6 1e-7 -1e-7 0]';
% dx = [0.5 -0.5 0.5 0.5e-3 -0.5e-3 0.5e-3 0]';
%initialize filter
j = 1;
xhat0(:,j) =  x0 + dx;
Deltax_minus(:,1) = zeros(7,1);
Deltax_plus(:,1) = ones(7,1);
while norm(Deltax_plus(:,end)) > 1e-7
% xhat(:,1) = xhat0(:,j);
Lambda(:,:,j) = pinv(P0);
N = pinv(P0)*Deltax_minus;%P0\Deltax_minus; %better alternative to inv(P0)*Deltax_minus
%propagate for entire measurement duration for each major iteration
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:7);
for i = 2:numel(t)
    %read next observation
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %integrate from t0 to ti
%     Phi = eye(n);
%     Phi_flat = reshape(Phi,n^2,1);
%     Z = [xhat0;Phi_flat];
%     tspan = [0:dt:t(i)];
%     [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    %determine station that produced meas
    statnum = find(~isnan(yi));
    if ~isempty(statnum)
        statnum = statnum(2)/2;
    end
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if statnum == 1%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:6) = Htilde_sc_rho_rhod(xhat(i,1:6),Xs_ECI(1,:));
           [hi(i,:),~] = predictmeas(xhat(i,1:6),Xs_ECI(1,:));
        elseif statnum == 2%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:6) = Htilde_sc_rho_rhod(xhat(i,1:6),Xs_ECI(2,:));
           [hi(i,:),~] = predictmeas(xhat(i,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:6) = Htilde_sc_rho_rhod(xhat(i,1:6),Xs_ECI(3,:));
           [hi(i,:),~] = predictmeas(xhat(i,1:6),Xs_ECI(3,:));
        end
        %accumulate observation
        ri(i,:) = yi(statnum*2-1:statnum*2) - hi(i,:);
        if isnan(hi(i,:))
            ri(i,:) = zeros(1,2);
        end
        Phi_flat_ti_t0 = x(i,n+1:end);
        Phi_ti_t0 = reshape(Phi_flat_ti_t0,n,n);
        Lambda(:,:,i) = Lambda(:,:,i-1) + (Hi*Phi_ti_t0)'*inv(R)*Hi*Phi_ti_t0;
        N(:,i) = N(:,i-1) + (Hi*Phi_ti_t0)'*inv(R)*ri(i,:)';
    else    %no measurements occured during the timestep
        Lambda(:,:,i) = Lambda(:,:,i-1);
        N(:,i) = N(:,i-1);
    end
end
%solve normal equations
Deltax_plus(:,j) = Lambda(:,:,i)\N(:,i);
P0_plus(:,:,j) = inv(Lambda(:,:,i));

xhat0(:,j+1) = xhat0(:,j) + Deltax_plus(:,j);
Deltax_minus(:,j+1) = Deltax_minus(:,j) - Deltax_plus(:,j);
j = j + 1;
end
    
for i = 1:numel(t)
    P_plus(:,:,i) = inv(Lambda(:,:,i)); %might need to implement square root method here due to singular matrix
end
    
    