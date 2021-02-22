% ASEN 6080
% Project 1
% NL Batch Estimation

tic

clear

load('observations.mat')

% pick to run in METERS or KILOMETERS
units_prompt = 'Enter 1 for meters, 0 for kilometers: ';
METER = input(units_prompt);

% pick to run station positions in ECEF or ECI
ECEF_prompt = 'Enter 1 for stations in ECEF, 0 for ECI: ';
ECEF = input(ECEF_prompt);

% pick which station to fix, or none
stationfix_prompt = 'Enter station number (101, 337, or 394) to be fixed, or 0 for none: ';
fixed_station = input(stationfix_prompt);

% choose one or both measurements
meas_prompt = 'Enter 2 to process both measurements, 1 for rho only, or 0 for rhodot only: ';
meas_mode = input(meas_prompt);

% other batch settings
tol = 1e-7;
maxMajor = 20;

% rearrange data into array with NaNs for timesteps w/o observations
dt = observations(2,1) - observations(1,1);
time = (0:dt:observations(end,1))';
obs_wNaN = NaN(length(time),4);
nn = 1;
for mm = 1:length(time)
    if time(mm) == observations(nn,1)
        obs_wNaN(mm,:) = observations(nn,:);
        nn = nn + 1;
    else
        obs_wNaN(mm,1) = time(mm);
    end
end

if METER == 1
% unpack data
IDs = obs_wNaN(:,2)';   % station IDs
ranges = obs_wNaN(:,3)';% range data, meters
RRs = obs_wNaN(:,4)';   % range rate data, m/s
y_data = [ranges; RRs];

% initial conditions, other given params
dt = time(2) - time(1);     % timestep, seconds
mu = 3.986004415E14;        % m^3/s^2
j2 = 1.082626925638815e-3;  % dimensionless
Re = 6378136.3;             % radius of Earth, meters
omegaE = 7.2921158553e-5;   % angular rate b/w ECI and ECEF, rad/s

H = 88667;                  % parameter for calculating rho, meters
rho0 = 3.614e-13;           % initial atmospheric density, kg/m^3
r0 = 700000.0 + Re;         % initial satellite radius, meters
A = 3.0;                    % cross-sectional area of s/c, m^2
m = 970;                    % mass of s/c, kg
CD = 2.0;                   % coefficient of drag, dimensionless

R1ecef = [-5127510.0 -3794160.0 0.0]';       % station 101 position, ECEF, meters
R2ecef0 = [3860910.0 3238490.0 3898094.0]';  % initial guess, station 337
R3ecef0 = [549505.0 -1380872.0 6182197.0]';  % initial guess, station 394

x0sc = [757700 5222607 4851500 2213.21...
        4678.34 -5371.3]';              % initial s/c state
dx0 = [100 -100 100 10 -10 10]';        % devation from given initial s/c state
%x0sc = x0sc + dx0;
    
n = 18;                                 % number of states
Phi0 = eye(n);                          % initial STM
Phi0_flat = reshape(Phi0,n^2,1);        % prepare STM to pass into ode45

X0 = [x0sc; mu; j2; CD; R1ecef;...      % 18x1 initial state vector
      R2ecef0; R3ecef0];                % station initial positions are equal in ECEF and ECI
  
X0_wPhi = [X0; Phi0_flat];              % attach Phi to state vector
  
P0 = diag([1e6 1e6 1e6 1e6 1e6 1e6 1e20...
           1e6 1e6 1e-10 1e-10 1e-10 1e6...
           1e6 1e6 1e6 1e6 1e6]);       % a priori covariance
       
Ri = diag([1e-4 1e-6]);                 % measurement noise
                                        % sigma rho = 1 cm
                                        % sigma dhro = 1 mm/s
else
% run in KILOMETERS
% unpack data
IDs = obs_wNaN(:,2)';        % station IDs
ranges = obs_wNaN(:,3)'/1000;% range data, km
RRs = obs_wNaN(:,4)'/1000;   % range rate data, km/s
y_data = [ranges; RRs];

% initial conditions, other given params
dt = time(2) - time(1);     % timestep, seconds
mu = 3.986004415E5;         % km^3/s^2
j2 = 1.082626925638815e-3;  % dimensionless
Re = 6378136.3/1000;        % radius of Earth, km
omegaE = 7.2921158553e-5;   % angular rate b/w ECI and ECEF, rad/s

H = 88667/1000;             % parameter for calculating rho, km
rho0 = 3.614e-13*1e9;       % initial atmospheric density, kg/km^3
r0 = 700000.0/1000 + Re;    % initial satellite radius, km
A = 3.0/1e6;                % cross-sectional area of s/c, km^2
m = 970;                    % mass of s/c, kg
CD = 2.0;                   % coefficient of drag, dimensionless

R1ecef = [-5127510.0 -3794160.0 0.0]'/1000;       % station 101 position, ECEF, km
R2ecef0 = [3860910.0 3238490.0 3898094.0]'/1000;  % initial guess, station 337
R3ecef0 = [549505.0 -1380872.0 6182197.0]'/1000;  % initial guess, station 394

x0sc = [757700 5222607 4851500 2213.21...
        4678.34 -5371.3]'/1000;         % initial s/c state
dx0 = [100 -100 100 10 -10 10]'/1000;   % devation from given initial s/c state
x0sc = x0sc + dx0;
    
n = 18;                                 % number of states
Phi0 = eye(n);                          % initial STM
Phi0_flat = reshape(Phi0,n^2,1);        % prepare STM to pass into ode45

X0 = [x0sc; mu; j2; CD; R1ecef;...      % 18x1 initial state vector
      R2ecef0; R3ecef0];                % station initial positions are equal in ECEF and ECI

X0_wPhi = [X0; Phi0_flat];              % attach Phi to state vector
  
P0 = diag([1 1 1 1 1 1 1e2...
           1e6 1e6 1e-16 1e-16 1e-16...
           1 1 1 1 1 1]);               % a priori covariance
       
       
% P0 = diag([1e-6 1e-6 1e-6 1e-10 1e-10 1e-10...
%            1e-6 1e-6 1e-5 5e-6 5e-6 5e-6...
%            5e-6 5e-6 5e-6 5e-6 5e-6 5e-6]);   % matching covariance

% P0 = diag([1e-6 1e-6 1e-6 1e-10 1e-10 1e-10...
%            1e-6 1e-6 1e-5 5e-6 5e-6 5e-6...
%            5e-6 5e-6 5e-6 5e-6 5e-6 5e-6])/1e20;   % too-small covariance
       
Ri = diag([1e-10 1e-12]);               % measurement noise
                                        % sigma rho = 1 cm
                                        % sigma dhro = 1 mm/s
end                                       
%=============================================================================================
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % set integration tolerances
Phi_ti_t0 = zeros(n,n,length(time));             % pre-allocate

% propagate over entire simulation
[~,X_wPhi_0] = ode45(@(t,X0_wPhi) proj1_ODE(t,X0_wPhi,H,A,m,rho0,r0,Re,ECEF),time, X0_wPhi, options);                     

for kk = 1:length(X_wPhi_0)
    Phi_i_im1_nxn = reshape(X_wPhi_0(kk,n+1:end),n,n);
    Phi_ti_t0(:,:,kk) = Phi_i_im1_nxn(:,:);         % store STMs for every time step
end

% pre-allocate
dX_plus = zeros(n,maxMajor);
dX_minus = zeros(n,maxMajor);
normdX_plus = zeros(1,maxMajor);
Xhat0 = zeros(n,maxMajor);
P0_plus = zeros(n,n,maxMajor);

% initialize
BatchSize = length(time)-1;
jj = 1;
dX_minus(:,jj) = zeros(n,1);% redundant to initialize w/ zeros, but keeping consistent to block diagram
Xhat0(:,1) = X0;

% begin major interations
while jj <= 4

    % pre-allocate
    ri = NaN(2,length(time));
    Hi = NaN(2,n,length(time));
    Lambda = zeros(n,n,length(time));
    N = zeros(n,length(time));

    % initialize
    Xhat(:,1) = [Xhat0(:,jj); Phi0_flat];
    Lambda(:,:,jj) = pinv(P0);
    N(:,jj) = pinv(P0)*dX_minus(:,jj);
    
    % propagate entire arc
    [~,X_wPhi] = ode45(@(t,Xhat) proj1_ODE(t,Xhat,H,A,m,rho0,r0,Re,ECEF),time,Xhat,options);
    Xhat = X_wPhi(:,1:n);

    % begin minor iterations
    for ii = 1:BatchSize
        if ECEF == 1
            % prepare ECI to ECEF rotation matrix
            tha = omegaE*360/2/pi*time(ii+1); % theta_GMST

            % transformation from ECI to ECEF
            C_eci2ecef = [cosd(tha) sind(tha) 0;
                         -sind(tha) cosd(tha) 0;
                             0          0     1];

            % transform station positions to ECI
            R1eci = C_eci2ecef'*R1ecef;
            R2eci = C_eci2ecef'*Xhat(ii+1,13:15)';
            R3eci = C_eci2ecef'*Xhat(ii+1,16:18)';

        else % if ECI
            R1eci = Xhat(ii+1,10:12)';
            R2eci = Xhat(ii+1,13:15)';
            R3eci = Xhat(ii+1,16:18)';
        end
        
        % read next observation
        yi = y_data(:,ii+1);
        yi = yi(~isnan(yi));

        % determine which station(s) is recording measurements, if any
        IDi = IDs(:,ii+1);      % extract station IDs that are recording
        IDi = IDi(~isnan(IDi+1)); % remove NaN
        if size(IDi,1) == 1     % if 1 measurement
            if IDi == 101       % if station 1 recording
                yhati = obs_prediction(Xhat(ii+1,1:6),R1eci);      % meas pred
                Hi_sc = Htilde_sc_rho_rhod(Xhat(ii+1,1:6),R1eci);  % meas sensitivity w.r.t 1st 9 elements of state
                Hi_obs = Htilde_obs_rho_rhod(Xhat(ii+1,1:6),R1eci);% measurement sensitivity w.r.t. recording station
                if fixed_station == 101
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,9)];
                else
                    Hi(:,:,ii+1) = [Hi_sc Hi_obs zeros(2,6)];
                end
            elseif IDi == 337   % if station 2 recording
                yhati = obs_prediction(Xhat(ii+1,1:6),R2eci);      % meas pred
                Hi_sc = Htilde_sc_rho_rhod(Xhat(ii+1,1:6),R2eci);  % meas sensitivity
                Hi_obs = Htilde_obs_rho_rhod(Xhat(ii+1,1:6),R2eci);% measurement sensitivity w.r.t. recording station
                if fixed_station == 337
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,9)];
                else
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,3) Hi_obs zeros(2,3)];
                end 
            else                % if station 3 recording
                yhati = obs_prediction(Xhat(ii+1,1:6),R3eci);      % meas pred
                Hi_sc = Htilde_sc_rho_rhod(Xhat(ii+1,1:6),R3eci);  % meas sensitivity
                Hi_obs = Htilde_obs_rho_rhod(Xhat(ii+1,1:6),R3eci);% measurement sensitivity w.r.t. recording station
                if fixed_station == 394
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,9)];
                else
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,6) Hi_obs];
                end
            end
            ri(:,ii+1) = yi - yhati; % measurement residual
            % accumulate observation
            if meas_mode == 1
                Lambda(:,:,ii+1) = Lambda(:,:,ii) + (Hi(1,:,ii+1)*Phi_ti_t0(:,:,ii+1))'*((Ri(1,1))\(Hi(1,:,ii+1)))*Phi_ti_t0(:,:,ii+1);
                N(:,ii+1) = N(:,ii) + (Hi(1,:,ii+1)*Phi_ti_t0(:,:,ii+1))'*((Ri(1,1))\(ri(1,ii+1)));
            elseif meas_mode == 0
                Lambda(:,:,ii+1) = Lambda(:,:,ii) + (Hi(2,:,ii+1)*Phi_ti_t0(:,:,ii+1))'*((Ri(2,2))\(Hi(2,:,ii+1)))*Phi_ti_t0(:,:,ii+1);
                N(:,ii+1) = N(:,ii) + (Hi(2,:,ii+1)*Phi_ti_t0(:,:,ii+1))'*((Ri(2,2))\(ri(2,ii+1)));
            else
                Lambda(:,:,ii+1) = Lambda(:,:,ii) + (Hi(:,:,ii+1)*Phi_ti_t0(:,:,ii+1))'*((Ri)\(Hi(:,:,ii+1)))*Phi_ti_t0(:,:,ii+1);
                N(:,ii+1) = N(:,ii) + (Hi(:,:,ii+1)*Phi_ti_t0(:,:,ii+1))'*((Ri)\(ri(:,ii+1)));
            end

        else    % if no measurements
            Lambda(:,:,ii+1) = Lambda(:,:,ii);
            N(:,ii+1) = N(:,ii);
        end
    end
    
    % solve "normal equations"
    dX_plus(:,jj) = pinv(Lambda(:,:,ii+1))*N(:,ii+1);
    normdX_plus(jj) = norm(dX_plus(:,jj));
    P0_plus(:,:,jj) = pinv(Lambda(:,:,ii+1)); % keep track of covariances of each initial guess
    
    % check to see if tolerance met; how large is the adjustment to the initial state?
    if normdX_plus(jj) < tol
        break
    end

    Xhat0(:,jj+1) = Xhat0(:,jj) + dX_plus(:,jj);        % adjust initial guess
    dX_minus(:,jj+1) = dX_minus(:,jj) - dX_plus(:,jj);  % keep track of adjustments
    jj = jj + 1;
    clear Xhat
end

toc

% Dx_LS
% propagate final estimate
[~,X_wPhi] = ode45(@(t,Xhat0) proj1_ODE(t,Xhat0,H,A,m,rho0,r0,Re,ECEF),time, [Xhat0(:,jj); Phi0_flat], options);

% find difference in state
scale = 1000; % scale residuals (from km)
Dx_LS = X_wPhi - X_wPhi_0;
Dx_LS = Dx_LS(:,1:n);
Dconsts = Dx_LS(1,7:9);         % difference in constants mu, J2, CD
Drs1 = Dx_LS(1,10:12)*scale;    % difference in station 1 ECEF, meters
Drs2 = Dx_LS(1,13:15)*scale;    % difference in station 1 ECEF, meters
Drs3 = Dx_LS(1,16:18)*scale;    % difference in station 1 ECEF, meters
digits = 10;

Dx_formatted = [round(Dx_LS(1,1:6)',digits)*scale; round(Dx_LS(1,7:9)',digits); round(Dx_LS(1,10:end)',digits)*scale]

% covariance over last major iteration
P_plus = zeros(n,n,length(time));
P_plus(:,:,1) = P0;
for kk = 1:length(time)
    P_plus(:,:,kk) = pinv(Lambda(:,:,kk));
end

% RMS values of residuals
% remove NaNs
ri_rho = ri(1,:);
ri_rho = ri_rho(~isnan(ri_rho));
ri_drho = ri(2,:);
ri_drho = ri_drho(~isnan(ri_drho));

RMS_rho = sqrt((sum(ri_rho(1,:).^2,'all'))/length(ri_rho))*scale*100;     % m
RMS_drho = sqrt((sum(ri_drho(1,:).^2,'all'))/length(ri_drho))*scale*1000;  % m/s

RMS = round([RMS_rho; RMS_drho],6)

hours = time/3600;

% plots
figure(1)   % s/c state
marker = 'b-';
ax1 = subplot(6,1,1);
plot(hours,Dx_LS(:,1),marker)
ylabel('\Deltax (km)')
grid on
ax2 = subplot(6,1,2);
plot(hours,Dx_LS(:,2),marker)
ylabel('\Deltay (km)')
grid on
ax3 = subplot(6,1,3);
plot(hours,Dx_LS(:,3),marker)
ylabel('\Deltaz (km)')
grid on
ax4 = subplot(6,1,4);
plot(hours,Dx_LS(:,4)*scale,marker)
ylabel('\Deltaxdot (m/s)')
grid on
ax5 = subplot(6,1,5);
plot(hours,Dx_LS(:,5)*scale,marker)
ylabel('\Deltaydot (m/s)')
grid on
ax6 = subplot(6,1,6);
plot(hours,Dx_LS(:,6)*scale,marker)
ylabel('\Deltazdot (m/s)')
grid on
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
sgtitle('\Deltax_{LS}, Spacecraft State')


hours = time/3600;

% build separate arrays for each station's residuals
% to make plotting by color easier
ri101 = NaN(2,length(ri));
ri337 = NaN(2,length(ri));
ri394 = NaN(2,length(ri));

scale = 1; % scale residuals (from km or m)
if (METER == 0) && (scale == 1000)  % update plot label units
    label = '';
elseif (METER == 0) && (scale == 1)
    label = 'k';
elseif (METER == 1) && (scale == 1)
    label = '';
end

for kk = 1:length(ri)
    if IDs(kk) == 101
        ri101(:,kk) = ri(:,kk)*scale;
    elseif IDs(kk) == 337
        ri337(:,kk) = ri(:,kk)*scale;
    else
        ri394(:,kk) = ri(:,kk)*scale;
    end
end

% measurement residuals (yi - yhati)
figure(44)
subplot(2,1,1)
plot(hours,ri101(1,:),'ro')
hold on
plot(hours,ri337(1,:),'bo')
plot(hours,ri394(1,:),'go')
hold off
grid on
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
ylabel(['Rho Residual (',label,'m)'])
legend('Station 101', 'Station 337','Station 394')
subplot(2,1,2)
plot(hours,ri101(2,:),'ro')
hold on
plot(hours,ri337(2,:),'bo')
plot(hours,ri394(2,:),'go')
hold off
grid on
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
ylabel(['dRho Residual (',label,'m/s)'])
sgtitle(['Post-Fit Measurement Residuals, Batch, Station ', num2str(fixed_station),...
        ' Fixed, ',num2str(jj-1),' Iterations'])

%%
% covariance ellipseoids at the final time
figure(100)
ellipsoid(0,0,0,3*sqrt(P0_plus(1,1,jj-1)),3*sqrt(P0_plus(2,2,jj-1)),3*sqrt(P0_plus(3,3,jj-1)))
pbaspect([1 1 1])
xlabel('r_x (m)')
ylabel('r_y (m)')
zlabel('r_z (m)')
title('Covariance Ellipsoid, S/C Position, Batch')

figure(101)
ellipsoid(0,0,0,3*sqrt(P0_plus(4,4,jj-1)),3*sqrt(P0_plus(5,5,jj-1)),3*sqrt(P0_plus(6,6,jj-1)))
pbaspect([1 1 1])
xlabel('v_x (m/s)')
ylabel('v_y (m/s)')
zlabel('v_z (m/s)')
title('Covariance Ellipsoid, S/C Velocity, Batch')
%%
% covariance plots
LWsig = 1.5;
pos = 40;
vel = 8e-2;

figure(55)
ax1 = subplot(6,1,1);
plot(hours,3*sqrt(squeeze(P_plus(1,1,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(1,1,:)))*1e3,'r:','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('x (m)')
xlim([hours(1) hours(end)])
ylim([-pos pos])

ax2 = subplot(6,1,2);
plot(hours,3*sqrt(squeeze(P_plus(2,2,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(2,2,:)))*1e3,'r:','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('y (m)')
xlim([hours(1) hours(end)])
ylim([-pos pos])

ax3 = subplot(6,1,3);
plot(hours,3*sqrt(squeeze(P_plus(3,3,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(3,3,:)))*1e3,'r:','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('z (m)')
xlim([hours(1) hours(end)])
ylim([-pos pos])

ax4 = subplot(6,1,4);
plot(hours,3*sqrt(squeeze(P_plus(4,4,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(4,4,:)))*1e3,'r:','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('xdot (m/s)')
xlim([hours(1) hours(end)])
ylim([-vel vel])

ax5 = subplot(6,1,5);
plot(hours,3*sqrt(squeeze(P_plus(5,5,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(5,5,:)))*1e3,'r:','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('ydot (m/s)')
xlim([hours(1) hours(end)])
ylim([-vel vel])

ax6 = subplot(6,1,6);
plot(hours,3*sqrt(squeeze(P_plus(6,6,:)))*1e3,'r:','LineWidth',LWsig)
hold on
plot(hours,-3*sqrt(squeeze(P_plus(6,6,:)))*1e3,'r:','LineWidth',LWsig)
hold off
grid on
xlabel('Time (hours)')
ylabel('zdot (m/s)')
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')
xlim([hours(1) hours(end)])
ylim([-vel vel])
sgtitle('Spacecraft State 3\sigma Covariance Bounds, Batch, Range-Rate Only')