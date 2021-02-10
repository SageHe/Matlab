% ASEN 6080
% Homework 2
% Classic Kalman Filter (CKF)

tic

clear

% load('y_nonoise.mat')
load('y_noisy.mat')
load('StationStates.mat')
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
x0 = kep2car(OEs, mu, angle_type);  % 1x6


n = 7;                          % number of states
Phi0 = eye(n);                  % initial STM
Phi0_flat = reshape(Phi0,n^2,1);% prepare STM to pass into ode45
X0 = [x0'; j2; Phi0_flat];      % prepare (n^2+n)x1 vector

%dx0 = [0.8 0.6 -0.2 1e-6 1e-7 -6e-7]'; % a priori state deviation
dx0 = [.5 -.5 .5 .5e-3 -.5e-3 .5e-3]';
%dx0 = zeros(6,1);
x0p = x0 + dx0'; % 1x6
X0p = [x0p'; j2; Phi0_flat];

P0 = diag([1 1 1 1e-6 1e-6 1e-6]);      % a priori covariance
Ri = diag([1e-3 1e-12]);    % measurement noise
                            % sigma rho = 1 km
                            % sigma dhro = 1 mm/s
                            
% set tolerances
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,Xp_wPhi] = ode45(@keplerJ2_wPhi_ODE,[0 time(end)], X0p, options);

maxMajor = 2;
tol = 1e-7;
NoMeas = length(time)-1;                            
jj = 1;

% pre-allocate
xNL0 = zeros(n-1,maxMajor);
Dx0_plus = zeros(n-1,maxMajor);

% initialize
xNL0(:,jj) = x0p';

% begin major iterations
while jj < maxMajor                            

    % pre-allocate
    Hi = zeros(2,6,length(time));
    Ki = zeros(6,2,length(time));
    xNL = zeros(length(time),n-1);
    xNL_wPhi = zeros(length(time),n^2+n);
    DX_plus = zeros(n-1,length(time));
    DX_minus = zeros(n-1,length(time));
    P_plus = zeros(n-1,n-1,length(time));
    P_minus = zeros(n-1,n-1,length(time));
    ri = zeros(2,length(time));

    % initialize
    DX_plus(1:6,1) = Dx0_plus(:,jj);    % redundant but keep consistent with block diagram
    P_plus(:,:,1) = P0;
    xNL(1,:) = xNL0(:,jj)';
    xNL_wPhi(1,:) = X0p';

    % begin minor iterations
    for ii = 1:NoMeas

        % read next observation
        yi = y_data(:,ii+1);
        yi = yi(~isnan(yi));

        % propagate from i-1 to i
        [~,X_wPhi] = ode45(@keplerJ2_wPhi_ODE,[0 dt], [xNL(ii,:)'; j2; Phi0_flat], options);
        xNL_wPhi(ii+1,:) = X_wPhi(end,:);
        xNL(ii+1,:) = xNL_wPhi(ii+1,1:6);
        Phi_im1_i = reshape(xNL_wPhi(ii+1,n+1:end),n,n);
        Phi_im1_i = Phi_im1_i(1:6,1:6);

        % time update
        DX_minus(:,ii+1) = Phi_im1_i*DX_plus(:,ii);
        P_minus(:,:,ii+1) = Phi_im1_i*P_plus(:,:,ii)*Phi_im1_i'; % plus P.N.

        % process observations
        % determine which station(s) is recording measurements, if any
        IDi = IDs(:,ii+1);      % extract station IDs that are recording
        IDi = IDi(~isnan(IDi)); % remove NaN
        if size(IDi,1) == 1     % if 1 measurement
            if IDi == 1         % if station 1 recording
                yhati = obs_prediction(xNL(ii+1,:),X1eci(ii+1,:));  % meas pred
                Hi(:,:,ii+1) = Htilde_sc_rho_rhod(xNL(ii+1,:),X1eci(ii+1,:)); % meas sensitivity
            elseif IDi == 2     % if station 2 recording
                yhati = obs_prediction(xNL(ii+1,:),X2eci(ii+1,:));  % meas pred
                Hi(:,:,ii+1) = Htilde_sc_rho_rhod(xNL(ii+1,:),X2eci(ii+1,:)); % meas sensitivity
            else                % if station 3 recording
                yhati = obs_prediction(xNL(ii+1,:),X3eci(ii+1,:));  % meas pred
                Hi(:,:,ii+1) = Htilde_sc_rho_rhod(xNL(ii+1,:),X3eci(ii+1,:)); % meas sensitivity
            end

            ri(:,ii+1) = yi - yhati; % measurement residual
            Ki(:,:,ii+1) = P_minus(:,:,ii+1)*Hi(:,:,ii+1)'/(Hi(:,:,ii+1)*P_minus(:,:,ii+1)*Hi(:,:,ii+1)' + Ri);   % Kalman gain

            % measurement updates
            DX_plus(:,ii+1) = DX_minus(:,ii+1) + Ki(:,:,ii+1)*(ri(:,ii+1) - Hi(:,:,ii+1)*DX_minus(:,ii+1));
            P_plus(:,:,ii+1) = (eye(6) - Ki(:,:,ii+1)*Hi(:,:,ii+1))*P_minus(:,:,ii+1)*(eye(6) - Ki(:,:,ii+1)*Hi(:,:,ii+1))' + Ki(:,:,ii+1)*Ri*Ki(:,:,ii+1)';

        elseif size(IDi,2) == 2     % if 2 measurements
            for kk = 1:2
                if IDi(kk) == 1     % if station 1 recording
                    yhati = obs_prediction(xNL(ii+1,:),X1eci(ii+1,:)); % meas pred
                    Hi(:,:,ii+1) = Htilde_sc_rho_rhod(xNL(ii+1,:),X1eci(ii+1,:));% meas sensitivity
                elseif IDi(kk) == 2 % if station 2 recording
                    yhati = obs_prediction(xNL(ii+1,:),X2eci(ii+1,:)); % meas pred
                    Hi(:,:,ii+1) = Htilde_sc_rho_rhod(xNL(ii+1,:),X2eci(ii+1,:));% meas sensitivity
                else                % if station 3 recording
                    yhati = obs_prediction(xNL(ii+1,:),X3eci(ii+1,:)); % meas pred
                    Hi(:,:,ii+1) = Htilde_sc_rho_rhod(xNL(ii+1,:),X3eci(ii+1,:));% meas sensitivity
                end

                ri(:,ii+1) = yi - yhati; % measurement residual
                Ki(:,:,ii+1) = P_minus(:,:,ii+1)*Hi(:,:,ii+1)'/(Hi(:,:,ii+1)*P_minus(:,:,ii+1)*Hi(:,:,ii+1)' + Ri);   % Kalman gain

                % measurement updates
                DX_plus(:,ii+1) = DX_minus(:,ii+1) + Ki(:,:,ii+1)*(ri(:,ii+1) - Hi(:,:,ii+1)*DX_minus(:,ii+1));
                P_plus(:,:,ii+1) = (eye(6) - Ki(:,:,ii+1)*Hi(:,:,ii+1))*P_minus(:,:,ii+1)*(eye(6) - Ki(:,:,ii+1)*Hi(:,:,ii+1))' + Ki(:,:,ii+1)*Ri*Ki(:,:,ii+1)';
            end

        else    % if no measurements
            DX_plus(:,ii+1) = DX_minus(:,ii+1);
            P_plus(:,:,ii+1) = P_minus(:,:,ii+1);
        end
    end

    % check for convergence
    if norm(DX_plus(:,ii+1)) < tol
        break
    end

    % if no convergence
    % propagate from 0 to i to get STM from 0 to i
    [~,X] = ode45(@keplerJ2_wPhi_ODE,[0 time(end)], [xNL0(:,jj); j2; Phi0_flat], options);
    xNL_wPhi_i0 = X(end,:);
    Phi_i_0 = reshape(xNL_wPhi_i0(1,n+1:end),n,n);
    Phi_i_0 = Phi_i_0(1:6,1:6);
    
    % invert STM to propagate from i to 0 to get new guess for Deltax_plus
    Dx0_plus(:,jj+1) = Phi_i_0\DX_plus(:,ii+1);
    xNL0(:,jj+1) = xNL0(:,jj) + Dx0_plus(:,jj+1);   % new initial NL state guess
    
    jj = jj + 1;
    
end

dx_hat = xNL0(:,jj) - xNL0(:,1);
disp('Estimated deviation from true initial state: ')
disp(dx_hat)
toc
%%

[~,X_hat_wPhi] = ode45(@keplerJ2_wPhi_ODE,time,[xNL0(:,jj); j2; Phi0_flat], options);
X_hat = X_hat_wPhi(:,1:6)';

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
pos = 1e2;
vel = 1e2;
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
sgtitle('Post-Fit Meaurement Residuals, CKF, Ignoring J3')

LWsig = 1.5;

figure(6)
ax1 = subplot(6,1,1);
plot(hours,3*sqrt(squeeze(P_plus(1,1,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(1,1,:)))*1e3,'r:','LineWidth',LWsig)
plot(hours,(X_hat(1,:)-xtrue(1,:))*1000,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltax (m)')
ylim([-pos,pos])

ax2 = subplot(6,1,2);
plot(hours,3*sqrt(squeeze(P_plus(2,2,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(2,2,:)))*1e3,'r:','LineWidth',LWsig)
plot(hours,(X_hat(2,:)-xtrue(2,:))*1000,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltay (m)')
ylim([-pos,pos])

ax3 = subplot(6,1,3);
plot(hours,3*sqrt(squeeze(P_plus(3,3,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(3,3,:)))*1e3,'r:','LineWidth',LWsig)
plot(hours,(X_hat(3,:)-xtrue(3,:))*1000,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltaz (m)')
ylim([-pos,pos])

ax4 = subplot(6,1,4);
plot(hours,3*sqrt(squeeze(P_plus(4,4,:)))*1e6,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(4,4,:)))*1e6,'r:','LineWidth',LWsig)
plot(hours,(X_hat(4,:)-xtrue(4,:))*1e6,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltaxdot (mm/s)')
ylim([-vel,vel])

ax5 = subplot(6,1,5);
plot(hours,3*sqrt(squeeze(P_plus(5,5,:)))*1e6,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(5,5,:)))*1e6,'r:','LineWidth',LWsig)
plot(hours,(X_hat(5,:)-xtrue(5,:))*1e6,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltaydot (mm/s)')
ylim([-vel,vel])

ax6 = subplot(6,1,6);
plot(hours,3*sqrt(squeeze(P_plus(6,6,:)))*1e6,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(6,6,:)))*1e6,'r:','LineWidth',LWsig)
plot(hours,(X_hat(6,:)-xtrue(6,:))*1e6,'b','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('\deltazdot (mm/s)')
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')
ylim([-vel,vel])
sgtitle('State Errors and 3\sigma Covariance Bounds, CKF, Ignoring J3')
%xlim([0 1.2])