% ASEN 6080
% Homework 2
% NL Batch-Lite Estimation

tic

clear

% load('y_nonoise.mat')
load('y_noisy.mat')
load('StationStates.mat')
1
load('Xtrue.mat')

X1eci = StationStates(:,1:6);
X2eci = StationStates(:,7:12);
X3eci = StationStates(:,13:18);
% y_data = y_nonoise(1:6,:); % select data set
% IDs = y_nonoise(7:9,:);    % station IDs
y_data = y_noisy(1:6,:); % select data set
IDs = y_noisy(7:9,:);    % station IDs
%===============================================
% load('y_nonoisej3.mat')
% load('y_noisyj3.mat')
% load('StationStatesj3.mat')
% load('Xtrue_j3.mat')
% 
% X1eci = StationStatesj3(:,1:6);
% X2eci = StationStatesj3(:,7:12);
% X3eci = StationStatesj3(:,13:18);
% % y_data = y_nonoisej3(1:6,:); % select data set
% % IDs = y_nonoisej3(7:9,:);    % station IDs
% y_data = y_noisyj3(1:6,:); % select data set
% IDs = y_noisyj3(7:9,:);    % station IDs
% xtrue = Xtrue_j3(:,1:6)';
%================================================

% Orbital elements/ICs
SMA = 10000;                    % km
ecc = 0.001;                    % dimensionless
inc = 40;                       % degrees
RAAN = 80;                      % degrees
arg = 40;                       % degrees
TA = 0;                         % degrees
mu = 3.986004415E5;             % km^3/s^2
period = 2*pi*sqrt(SMA^3/mu);   % period in seconds
j2 = 0.0010826269;              % dimensionless
dt = 10;                        % seconds
time = 0:dt:15*period;

% convert to Cartesian inertial coordinates
OEs = [SMA ecc inc RAAN arg TA];
angle_type = 'degrees';
x0 = kep2car(OEs, mu, angle_type);

n = 7;                          % number of states
Phi0 = eye(n);                  % initial STM
Phi0_flat = reshape(Phi0,n^2,1);% prepare STM to pass into ode45
X0 = [x0'; j2; Phi0_flat];      % prepare (n^2+n)x1 vector

P0 = diag([1 1 1 1e-6 1e-6 1e-6]); % a priori covariance
Ri = diag([1e-3 1e-12]);        % measurement noise
                                % sigma rho = 1 km
                                % sigma dhro = 1 mm/s
                            
dX0 = [0.8 0.6 -0.2 1e-6 1e-7 -6e-7]'; % a priori state deviation
%dX0 = zeros(6,1);
x0p = x0 + dX0';    % 1x6, initial "guess" of state

options = odeset('RelTol',1e-12,'AbsTol',1e-12); % set integration tolerances
Phi_ti_t0 = zeros(n-1,n-1,length(time));         % pre-allocate

% propagate over entire simulation
[t,X_wPhi] = ode45(@keplerJ2_wPhi_ODE,time, [x0p'; j2; Phi0_flat], options);
for kk = 1:length(X_wPhi)
    Phi_i_im1_7x7 = reshape(X_wPhi(kk,n+1:end),n,n);
    Phi_ti_t0(:,:,kk) = Phi_i_im1_7x7(1:6,1:6);         % store STMs for every time step
end

tol = 1e-7;
maxMajor = 20;
BatchSize = 14400;

% pre-allocate
dX_plus = zeros(n-1,maxMajor);
dX_minus = zeros(n-1,maxMajor);
normdX_plus = zeros(1,maxMajor);
Xhat0 = zeros(n-1,maxMajor);
P0_plus = zeros(n-1,n-1,maxMajor);

% initialize
jj = 1;
dX_minus(:,jj) = zeros(6,1);% redundant to initialize w/ zeros, but keeping consistent to block diagram
Xhat0(:,1) = x0p;

% begin major interations
while jj < maxMajor

    % pre-allocate
    Xhat = zeros(length(time),n-1);
    Xhat_wPhi = zeros(length(time),n^2+n);
    ri = zeros(2,length(time));
    Lambda = zeros(n-1,n-1,length(time));
    N = zeros(n-1,length(time));

    % initialize
    Xhat(1,:) = Xhat0(:,jj)';
    Xhat_wPhi(1,:) = X0';
    Lambda(:,:,jj) = inv(P0);
    N(:,jj) = P0\dX_minus(:,jj);

    % begin minor iterations loop
    for ii = 1:BatchSize

        % read next observation
        yi = y_data(:,ii+1);
        yi = yi(~isnan(yi));
        
        % propagate from ti-1 to ti
        [t,X_wPhi] = ode45(@keplerJ2_wPhi_ODE,[0 dt], [Xhat(ii,:)'; j2; Phi0_flat], options);
        Xhat_wPhi(ii+1,:) = X_wPhi(end,:);
        Xhat(ii+1,:) = Xhat_wPhi(ii+1,1:6);     

        % determine which station(s) is recording measurements, if any
        IDi = IDs(:,ii+1);      % extract station IDs that are recording
        IDi = IDi(~isnan(IDi)); % remove NaN
        if size(IDi,1) == 1     % if 1 measurement
            if IDi == 1         % if station 1 recording
                yhati = obs_prediction(Xhat(ii+1,:),X1eci(ii+1,:));  % meas pred
                Hi = Htilde_sc_rho_rhod(Xhat(ii+1,:),X1eci(ii+1,:)); % meas sensitivity
            elseif IDi == 2     % if station 2 recording
                yhati = obs_prediction(Xhat(ii+1,:),X2eci(ii+1,:));  % meas pred
                Hi = Htilde_sc_rho_rhod(Xhat(ii+1,:),X2eci(ii+1,:)); % meas sensitivity
            else                % if station 3 recording
                yhati = obs_prediction(Xhat(ii+1,:),X3eci(ii+1,:));  % meas pred
                Hi = Htilde_sc_rho_rhod(Xhat(ii+1,:),X3eci(ii+1,:)); % meas sensitivity
            end
            
            % accumulate observation
            ri(:,ii+1) = yi - yhati; % measurement residual
            Lambda(:,:,ii+1) = Lambda(:,:,ii) + (Hi*Phi_ti_t0(:,:,ii+1))'*((Ri)\(Hi))*Phi_ti_t0(:,:,ii+1);
            N(:,ii+1) = N(:,ii) + (Hi*Phi_ti_t0(:,:,ii+1))'*((Ri)\(ri(:,ii+1)));

        elseif size(IDi,2) == 2     % if 2 measurements
            for jj = 1:2
                if IDi(jj) == 1     % if station 1 recording
                    yhati = obs_prediction(Xhat(ii+1,:),X1eci(ii+1,:)); % meas pred
                    Hi = Htilde_sc_rho_rhod(Xhat(ii+1,:),X1eci(ii+1,:));% meas sensitivity
                elseif IDi(jj) == 2 % if station 2 recording
                    yhati = obs_prediction(Xhat(ii+1,:),X2eci(ii+1,:)); % meas pred
                    Hi = Htilde_sc_rho_rhod(Xhat(ii+1,:),X2eci(ii+1,:));% meas sensitivity
                else                % if station 3 recording
                    yhati = obs_prediction(Xhat(ii+1,:),X3eci(ii+1,:)); % meas pred
                    Hi = Htilde_sc_rho_rhod(Xhat(ii+1,:),X3eci(ii+1,:));% meas sensitivity
                end
                
                % accumulate observation
                ri(:,ii+1) = yi - yhati; % measurement residual
                Lambda(:,:,ii+1) = Lambda(:,:,ii) + (Hi*Phi_ti_t0(:,:,ii+1))'*((Ri)\(Hi))*Phi_ti_t0(:,:,ii+1);
                N(:,ii+1) = N(:,ii) + (Hi*Phi_ti_t0(:,:,ii+1))'*((Ri)\(ri(:,ii+1)));
                
            end

        else    % if no measurements
            Lambda(:,:,ii+1) = Lambda(:,:,ii);
            N(:,ii+1) = N(:,ii);
        end
    end
    
    % solve "normal equations"
    dX_plus(:,jj) = Lambda(:,:,ii+1)\N(:,ii+1);
    normdX_plus(jj) = norm(dX_plus(:,jj));
    P0_plus(:,:,jj) = inv(Lambda(:,:,ii+1));    % keep track of covariances of each initial guess
    
    % check to see if tolerance met; how large is the adjustment to the initial state?
    if normdX_plus(jj) < tol
        break
    end
    
    % if no convergence
    Xhat0(:,jj+1) = Xhat0(:,jj) + dX_plus(:,jj);        % adjust initial guess
    dX_minus(:,jj+1) = dX_minus(:,jj) - dX_plus(:,jj);  % keep track of adjustments
    jj = jj + 1;
end
toc
%%
[~,X_hat_wPhi] = ode45(@keplerJ2_wPhi_ODE,time,[Xhat0(:,2); j2; Phi0_flat], options);
X_hat = X_hat_wPhi(:,1:6)';

P_plus = zeros(n-1,n-1,length(time));
P_plus(:,:,1) = P0;
for kk = 1:length(time)
    P_plus(:,:,kk) = inv(Lambda(:,:,kk));
end

errors = X_hat-xtrue;

RMS_x = sqrt((sum(errors(1,:).^2,'all'))/length(errors))*1e3;   % m
RMS_y = sqrt((sum(errors(2,:).^2,'all'))/length(errors))*1e3;   % m
RMS_z = sqrt((sum(errors(3,:).^2,'all'))/length(errors))*1e3;   % m
RMS_dx = sqrt((sum(errors(4,:).^2,'all'))/length(errors))*1e6;  % mm/s
RMS_dy = sqrt((sum(errors(5,:).^2,'all'))/length(errors))*1e6;  % mm/s
RMS_dz = sqrt((sum(errors(6,:).^2,'all'))/length(errors))*1e6;  % mm/s

RMS_r = sqrt((sum(errors(1:3,:).^2,'all'))/(length(errors)*3))*1e3;
RMS_v = sqrt((sum(errors(4:6,:).^2,'all'))/(length(errors)*3))*1e6;

RMS = [RMS_x; RMS_y; RMS_z; RMS_dx; RMS_dy; RMS_dz; RMS_r; RMS_v]

%%
hours = time/3600;

figure(4)
subplot(2,1,1)
plot(hours,ri(1,:)*1000)
grid on
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
ylabel('Rho Residual (m)')
subplot(2,1,2)
plot(hours,ri(2,:)*10^6)
grid on
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
ylabel('dRho Residual (mm/s)')
sgtitle('Post-fit Residuals, Batch, Ignoring J3')

LWsig = 1.5;

pos = 1e2;
vel = 1e2;

figure(5)
ax1 = subplot(6,1,1);
plot(hours,3*sqrt(squeeze(P_plus(1,1,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(1,1,:)))*1e3,'r:','LineWidth',LWsig)
plot(hours,errors(1,:)*1e3,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltax (m)')
ylim([-pos pos])

ax2 = subplot(6,1,2);
plot(hours,3*sqrt(squeeze(P_plus(2,2,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(2,2,:)))*1e3,'r:','LineWidth',LWsig)
plot(hours,errors(2,:)*1e3,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltay (m)')
ylim([-pos pos])

ax3 = subplot(6,1,3);
plot(hours,3*sqrt(squeeze(P_plus(3,3,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(3,3,:)))*1e3,'r:','LineWidth',LWsig)
plot(hours,errors(3,:)*1e3,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltaz (m)')
ylim([-pos pos])

ax4 = subplot(6,1,4);
plot(hours,3*sqrt(squeeze(P_plus(4,4,:)))*1e6,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(4,4,:)))*1e6,'r:','LineWidth',LWsig)
plot(hours,errors(4,:)*1e6,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltaxdot (mm/s)')
ylim([-vel vel])

ax5 = subplot(6,1,5);
plot(hours,3*sqrt(squeeze(P_plus(5,5,:)))*1e6,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(5,5,:)))*1e6,'r:','LineWidth',LWsig)
plot(hours,errors(5,:)*1e6,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltady (mm/s)')
ylim([-vel vel])

ax6 = subplot(6,1,6);
plot(hours,3*sqrt(squeeze(P_plus(6,6,:)))*1e6,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(6,6,:)))*1e6,'r:','LineWidth',LWsig)
plot(hours,errors(6,:)*1e6,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltadz (mm/s)')
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')
ylim([-vel vel])
sgtitle('State Errors and 3\sigma Covariance Bounds, Batch, Ignoring J3')