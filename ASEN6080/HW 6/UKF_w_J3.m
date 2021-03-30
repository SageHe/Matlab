clear all;clc
tic
%load in measurements from HW1
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
% load('HW2_Nom_Measurements')
% load('y_noisy');
% load('J3_Measurements')
load('J3_Meas_Noised')
% rangemeas(14319:end,:) = []; %Clean up end of measurement file that has no measurements
% rrmeas(14319:end,:) = [];
rangemeasnoisy(14319:end,:) = []; %Clean up end of measurement file that has no measurements
rrmeasnoisy(14319:end,:) = [];
t(14319:end) = [];
truth_state(14319:end,:) = [];
% load('HW1_Noisy_Measurements')
mu = 3.986004415e5;
J2 = 1082.63e-6;
J3 = -2.5323e-6;
%process noise flag
PN = 1;
sigma = 5e-10;

muscale = 1e-4;
J2scale = 1e3;
J3scale = 1e6;
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
% rho_noise = sigma_rho*randn(size(rangemeas,1),1);
% rhodot_noise = sigma_rhodot*randn(size(rangemeas,1),1);
sigma_r = 1;
sigma_v = 1e-3;

% rangemeasnoisy = rangemeas + rho_noise;
% rrmeasnoisy = rrmeas + rhodot_noise;

%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2 1e-12 1e-12 1e-12]);
R_meas = diag([sigma_rho^2 sigma_rhodot^2]);

%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
x0 = [r;v;mu*muscale;J2*J2scale;J3*J3scale];
n = size(x0,1);
dt = t(2) - t(1);

% Phi = eye(size(x0,1));
% Phi_flat = reshape(Phi,size(x0,1)^2,1);
% Z = [x0;Phi_flat];
% tspan = [0:dt:t(end)];
% [~,X_true] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

% dx = zeros(7,1);
dx = [0.5 -0.5 0.5 0.5e-7 -0.5e-7 0.5e-7 0 0 0]';
% dx = [0.8 0.6 -0.2 1e-6 1e-7 -6e-7]';
%Initialize filter
x_hat0(1,:) = x0 + dx;
x_hat0 = x_hat0';
x_t_im1 = x_hat0;
x_hat_plus(:,1) = x_t_im1;
P_plus(:,:,1) = P0;

alpha = 1e-4;
beta = 2;
kappa = 0;
lambda = alpha^2*(kappa + n) - n;
gamma = sqrt(n + lambda);

Q_DT = zeros(n,n);
Qsigma = sigma*ones(1,3);
Q_DT(1:6,1:6) = calc_DT_PN(dt,Qsigma);

%compute weights
W_k_m(1) = lambda/(lambda + n);
W_k_m(2:2*n+1) = 1/(2*(lambda + n));

W_k_c(1) = (lambda/(lambda + n)) + (1 - alpha^2 + beta);
W_k_c(2:2*n+1) = 1/(2*(lambda + n));

ri = NaN(numel(t),6);

y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
for i = 2:numel(t)
    %read next observation
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %propagate sigma points
    X_k_plus(:,:,i-1) = zeros(2*n+1,n);
    S = sqrtm(P_plus(:,:,i-1));
    S = [S S];
    X_k_plus1 = x_hat_plus(:,i-1)';
    for j = 1:2*n
        if j <= n
            X_k_plus_loop(j,:) = x_hat_plus(:,i-1)' + gamma*S(:,j)';
        else 
            X_k_plus_loop(j,:) = x_hat_plus(:,i-1)' - gamma*S(:,j)';
        end
    end
    X_k_plus(:,:,i-1) = [X_k_plus1;X_k_plus_loop];
    Z = reshape(X_k_plus(:,:,i-1)',(n*(2*n+1)),1);
    tspan = [0:dt/2:10];
    [~,X_k_minus_flat] = ode45(@(t,Z) keplerJ2J3_ODE(t,Z,n),tspan,Z,opts);
    X_k_minus(:,:,i) = reshape(X_k_minus_flat(end,:),n,2*n+1)';
%     X_k_minus(:,:,i) = X_k_minus(:,:,i)';
    for j = 1:2*n+1
        X_k_minus_loop(j,:) = X_k_minus(j,:,i)*W_k_m(j);
    end
    x_hat_minus(:,i) = sum(X_k_minus_loop,1);
    P_minus(:,:,i) = zeros(n,n);
    for j = 1:2*n+1
%         P_minus_loop(:,:,j) = W_k_c(j)*((X_k_minus(j,:,i) - x_hat_minus(:,i))*(X_k_minus(j,:,i) - x_hat_minus(:,i))');
        P_minus(:,:,i) = P_minus(:,:,i) + W_k_c(j)*((X_k_minus(j,:,i)' - x_hat_minus(:,i))*(X_k_minus(j,:,i) - x_hat_minus(:,i)'));
    end
    
    if PN
        P_minus(:,:,i) = P_minus(:,:,i) + Q_DT;
    end

    %resample sigma points from calculated x_hat_minus and P_minus
    X_k_minus(:,:,i) = zeros(2*n+1,n);
    S = sqrtm(P_minus(:,:,i));
    S = [S S];
    X_k_minus1 = x_hat_minus(:,i)'; %CHECK IF X_HAT_MINUS IS NX1 OR 1XN
    for j = 1:2*n
        if j <=n
            X_k_minus_temp(j,:) = x_hat_minus(:,i)' + gamma*S(:,j)';
        else
            X_k_minus_temp(j,:) = x_hat_minus(:,i)' - gamma*S(:,j)';
        end
    end
    X_k_minus(:,:,i) = [X_k_minus1;X_k_minus_temp];
    
    %Measurement update step
    if ~all(isnan(yi))
        inds = find(~isnan(yi));
        for j = 1:2*n+1
            Y_k_all(j,:,i) = predictmeas(X_k_minus(j,:,i),Xs_ECI);
        end
        Y_k(:,:,i) = Y_k_all(:,inds,i);
        y_hat(i,:) = zeros(1,2);
        for j = 1:2*n+1
            y_hat(i,:) = y_hat(i,:) + W_k_m(j)*Y_k(j,:,i);
        end
        Pyy_minus(:,:,i) = zeros(2,2);
        Pxy_minus(:,:,i) = zeros(n,2);
        for j = 1:2*n+1
            Pyy_minus(:,:,i) = Pyy_minus(:,:,i) + W_k_c(j)*((Y_k(j,:,i) - y_hat(i,:))'*(Y_k(j,:,i) - y_hat(i,:)));
            Pxy_minus(:,:,i) = Pxy_minus(:,:,i) + W_k_c(j)*((X_k_minus(j,:,i)' - x_hat_minus(:,i))*(Y_k(j,:,i) - y_hat(i,:)));
        end
        Pyy_minus(:,:,i) = Pyy_minus(:,:,i) + R_meas;
        K(:,:,i) = Pxy_minus(:,:,i)*pinv(Pyy_minus(:,:,i));

        x_hat_plus(:,i) = x_hat_minus(:,i) + K(:,:,i)*(yi(inds) - y_hat(i,:))';
        P_plus(:,:,i) = P_minus(:,:,i) - K(:,:,i)*Pyy_minus(:,:,i)*K(:,:,i)';
        
        ri(i,inds) = yi(inds) - y_hat(i,:);
    else
        P_plus(:,:,i) = P_minus(:,:,i);
        x_hat_plus(:,i) = x_hat_minus(:,i);
    end
    i
    covbounds(i,:) = [3*sqrt(P_plus(1,1,i)) 3*sqrt(P_plus(2,2,i)) 3*sqrt(P_plus(3,3,i)) 3*sqrt(P_plus(4,4,i)) 3*sqrt(P_plus(5,5,i)) 3*sqrt(P_plus(6,6,i))]; 
    posnorm(i) = norm(x_hat_plus(1:3,i) - truth_state(i,1:3)');
    velnorm(i) = norm(x_hat_plus(4:6,i) - truth_state(i,4:6)');
    
end
toc
state_error = x_hat_plus(1:6,:) - truth_state(:,1:6)';
posRMS = sqrt((1/numel(t))*sum(posnorm.^2));
velRMS = sqrt((1/numel(t))*sum(velnorm.^2));

figure
subplot(3,1,1)
hold on
plot(t,state_error(1,:))
plot(t,covbounds(:,1),'--r')
plot(t,-covbounds(:,1),'--r')
xlabel('Time (s)')
ylabel('\delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
hold on
plot(t,state_error(2,:))
plot(t,covbounds(:,2),'--r')
plot(t,-covbounds(:,2),'--r')
xlabel('Time (s)')
ylabel('delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
hold on
plot(t,state_error(3,:))
plot(t,covbounds(:,3),'--r')
plot(t,-covbounds(:,3),'--r')
xlabel('Time (s)')
ylabel('\delta r_z (km)')
grid on
grid minor
sgtitle('UKF Position State Errors w/PN and J3')

figure
subplot(3,1,1)
hold on
plot(t,state_error(4,:))
plot(t,covbounds(:,4),'--r')
plot(t,-covbounds(:,4),'--r')
xlabel('Time (s)')
ylabel('\delta v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
hold on
plot(t,state_error(5,:))
plot(t,covbounds(:,5),'--r')
plot(t,-covbounds(:,5),'--r')
xlabel('Time (s)')
ylabel('delta v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
hold on
plot(t,state_error(6,:))
plot(t,covbounds(:,6),'--r')
plot(t,-covbounds(:,6),'--r')
xlabel('Time (s)')
ylabel('\delta v_z (km/s)')
grid on 
grid minor
sgtitle('UKF Velocity State Errors w/PN and J3')

figure
subplot(2,1,1)
hold on
plot(t,ri(:,1),'*')
plot(t,ri(:,3),'*')
plot(t,ri(:,5),'*')
xlabel('Time (s)')
ylabel('Range Residual (km)')
ylim([-0.01 0.01])
grid on
grid minor
subplot(2,1,2)
hold on
plot(t,ri(:,2),'*')
plot(t,ri(:,4),'*')
plot(t,ri(:,6),'*')
xlabel('Time (s)')
ylabel('Range Rate Residual (km/s)')
ylim([-1e-5 1e-5])
grid on
grid minor
sgtitle('UKF Residuals, w/ PN and J3')