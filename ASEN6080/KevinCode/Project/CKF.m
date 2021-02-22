% ASEN 6080
% Project 1
% CKF Estimation

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

% choose # of major iterations
iteration_prompt = 'Enter # of major iterations: ';
maxMajor = input(iteration_prompt);

tol = 1e-7;

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
      R2ecef0; R3ecef0];
  
X0_wPhi = [X0; Phi0_flat];              % attach Phi to state vector
  
P0 = diag([1e6 1e6 1e6 1e6 1e6 1e6 1e20...
           1e6 1e6 1e-10 1e-10 1e-10 1e6...
           1e6 1e6 1e6 1e6 1e6]);       % a priori covariance
       
Ri = diag([1e-4 1e-6]);                 % measurement noise
                                        % sigma rho = 1 cm
                                        % sigma dhro = 1 mm/s
else
% convert to KILOMETERS
% unpack data
IDs = obs_wNaN(:,2)';   % station IDs
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
%x0sc = x0sc + dx0;
    
n = 18;                                 % number of states
Phi0 = eye(n);                          % initial STM
Phi0_flat = reshape(Phi0,n^2,1);        % prepare STM to pass into ode45

X0 = [x0sc; mu; j2; CD; R1ecef;...      % 18x1 initial state vector
      R2ecef0; R3ecef0];
  
X0_wPhi = [X0; Phi0_flat];              % attach Phi to state vector
  
P0 = diag([1 1 1 1 1 1 1e2...
           1e6 1e6 1e-16 1e-16 1e-16...
           1 1 1 1 1 1]);               % a priori covariance
       
Ri = diag([1e-10 1e-12]);               % measurement noise
                                        % sigma rho = 1 cm
                                        % sigma dhro = 1 mm/s
end                                        
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % set integration tolerances
NoMeas = length(time)-1;                            
jj = 1;

% pre-allocate
xNL0 = zeros(n,maxMajor);
Dx0_plus = zeros(n,maxMajor);
normDX_plus = zeros(1,maxMajor);

% initialize
xNL0(:,jj) = X0;

% begin major iterations
while jj <= maxMajor

    % pre-allocate
    Hi = zeros(2,n,length(time));
    Ki = zeros(n,2,length(time));
    xNL = zeros(length(time),n);
    xNL_wPhi = zeros(length(time),n^2+n);
    DX_plus = zeros(n,length(time));
    DX_minus = zeros(n,length(time));
    P_plus = zeros(n,n,length(time));
    P_minus = zeros(n,n,length(time));
    ri = NaN(2,length(time));
    prefits = NaN(2,length(time));
    postfits = NaN(2,length(time));

    % initialize
    DX_plus(:,1) = Dx0_plus(:,jj);    % redundant but keep consistent with block diagram
    P_plus(:,:,1) = P0;
    xNL(1,:) = xNL0(:,jj)';
    xNL_wPhi(1,:) = X0_wPhi';
    X_wPhi = [xNL(1,:)'; Phi0_flat];

    % begin minor iterations
    for ii = 1:NoMeas

        % read next observation
        yi = y_data(:,ii+1);
        yi = yi(~isnan(yi));

        % propagate from i-1 to i
        [~,X_wPhi] = ode45(@(t,X0_wPhi) proj1_ODE(t,X0_wPhi,H,A,m,rho0,r0,Re,ECEF),[0 dt], X_wPhi, options);
        xNL_wPhi(ii+1,:) = X_wPhi(end,:);
        xNL(ii+1,:) = xNL_wPhi(ii+1,1:n);
        X_wPhi = [xNL(ii+1,:)'; Phi0_flat];
        Phi_im1_i = reshape(xNL_wPhi(ii+1,n+1:end),n,n);

        % time update
        DX_minus(:,ii+1) = Phi_im1_i*DX_plus(:,ii);
        P_minus(:,:,ii+1) = Phi_im1_i*P_plus(:,:,ii)*Phi_im1_i'; % plus P.N. if applicable

        if ECEF == 1
            % prepare ECI to ECEF rotation matrix
            tha = omegaE*360/2/pi*time(ii+1);    % theta_GMST

            % transformation from ECI to ECEF
            C_eci2ecef = [cosd(tha) sind(tha) 0;
                         -sind(tha) cosd(tha) 0;
                             0          0     1];

            % transform station positions to ECI
            R1eci = C_eci2ecef'*R1ecef;
            R2eci = C_eci2ecef'*xNL(ii+1,13:15)';
            R3eci = C_eci2ecef'*xNL(ii+1,16:18)';

        else % if ECI
            R1eci = xNL(ii+1,10:12)';
            R2eci = xNL(ii+1,13:15)';
            R3eci = xNL(ii+1,16:18)';
        end     
        
        % process observations
        % determine which station(s) is recording measurements, if any
        IDi = IDs(:,ii+1);      % extract station IDs that are recording
        IDi = IDi(~isnan(IDi)); % remove NaN
        if size(IDi,1) == 1     % if 1 measurement
            if IDi == 101       % if station 1 recording
                yhati = obs_prediction(xNL(ii+1,1:6),R1eci);        % meas pred
                Hi_sc = Htilde_sc_rho_rhod(xNL(ii+1,1:6),R1eci);    % meas sensitivity
                Hi_obs = Htilde_obs_rho_rhod(xNL(ii+1,1:6),R1eci);  % measurement sensitivity w.r.t. recording station
                if fixed_station == 101
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,9)];
                else
                    Hi(:,:,ii+1) = [Hi_sc Hi_obs zeros(2,6)];
                end            
            elseif IDi == 337   % if station 2 recording
                yhati = obs_prediction(xNL(ii+1,1:6),R2eci);        % meas pred
                Hi_sc = Htilde_sc_rho_rhod(xNL(ii+1,1:6),R2eci);    % meas sensitivity
                Hi_obs = Htilde_obs_rho_rhod(xNL(ii+1,1:6),R2eci);  % measurement sensitivity w.r.t. recording station
                if fixed_station == 337
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,9)];
                else
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,3) Hi_obs zeros(2,3)];
                end                
            else                % if station 3 recording
                yhati = obs_prediction(xNL(ii+1,1:6),R3eci);        % meas pred
                Hi_sc = Htilde_sc_rho_rhod(xNL(ii+1,1:6),R3eci);    % meas sensitivity
                Hi_obs = Htilde_obs_rho_rhod(xNL(ii+1,1:6),R3eci);  % measurement sensitivity w.r.t. recording station
                if fixed_station == 394
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,9)];
                else
                    Hi(:,:,ii+1) = [Hi_sc zeros(2,6) Hi_obs];
                end             
            end

            ri(:,ii+1) = yi - yhati; % measurement residual
            Ki(:,:,ii+1) = P_minus(:,:,ii+1)*Hi(:,:,ii+1)'/(Hi(:,:,ii+1)*P_minus(:,:,ii+1)*Hi(:,:,ii+1)' + Ri);   % Kalman gain

            % measurement updates
            DX_plus(:,ii+1) = DX_minus(:,ii+1) + Ki(:,:,ii+1)*(ri(:,ii+1) - Hi(:,:,ii+1)*DX_minus(:,ii+1));
            P_plus(:,:,ii+1) = (eye(n) - Ki(:,:,ii+1)*Hi(:,:,ii+1))*P_minus(:,:,ii+1)*(eye(n) - Ki(:,:,ii+1)*Hi(:,:,ii+1))' + Ki(:,:,ii+1)*Ri*Ki(:,:,ii+1)';

            % residuals
            prefits(:,ii+1) = ri(:,ii+1) - Hi(:,:,ii+1)*DX_minus(:,ii+1);
            postfits(:,ii+1) = ri(:,ii+1) - Hi(:,:,ii+1)*DX_plus(:,ii+1);
            
        else    % if no measurements
            DX_plus(:,ii+1) = DX_minus(:,ii+1);
            P_plus(:,:,ii+1) = P_minus(:,:,ii+1);
        end
    end

    % check for convergence
    normDX_plus(jj) = norm(DX_plus(:,ii+1));
    if norm(DX_plus(:,ii+1)) < tol
        break
    end

    % if no convergence
    % propagate from 0 to i to get STM from 0 to i
    [~,X] = ode45(@(t,X0_wPhi) proj1_ODE(t,X0_wPhi,H,A,m,rho0,r0,Re,ECEF),[0 time(end)], [xNL0(:,jj); Phi0_flat], options);
    xNL_wPhi_i0 = X(end,:);
    Phi_i_0 = reshape(xNL_wPhi_i0(1,n+1:end),n,n);
    
    % invert STM to propagate from i to 0 to get new guess for Deltax_plus
    Dx0_plus(:,jj+1) = pinv(Phi_i_0)*DX_plus(:,ii+1);
    xNL0(:,jj+1) = xNL0(:,jj) + Dx0_plus(:,jj+1);   % new initial NL state guess
    
    jj = jj + 1;
    
end

dx_hat = xNL0(:,jj) - xNL0(:,1);
% disp('Estimated deviation from true initial state: ')
% disp(dx_hat)
toc

% propagate over entire simulation for initial and final guesses
[~,X_wPhi_0] = ode45(@(t,X0_wPhi) proj1_ODE(t,X0_wPhi,H,A,m,rho0,r0,Re,ECEF),time, X0_wPhi, options);
[~,X_wPhi] = ode45(@(t,xNL0) proj1_ODE(t,xNL0,H,A,m,rho0,r0,Re,ECEF),time, [xNL0(:,jj); Phi0_flat], options);

% Dx_CKF
% find difference in state
scale = 1000; % scale residuals (from km)
Dx_CKF = X_wPhi - X_wPhi_0;
Dx_CKF = Dx_CKF(:,1:n);
Dconsts = Dx_CKF(1,7:9);         % difference in constants mu, J2, CD
Drs1 = Dx_CKF(1,10:12)*scale;    % difference in station 1 ECEF, meters
Drs2 = Dx_CKF(1,13:15)*scale;    % difference in station 1 ECEF, meters
Drs3 = Dx_CKF(1,16:18)*scale;    % difference in station 1 ECEF, meters

Dx_formatted = [round(Dx_CKF(1,1:6)',10)*scale; round(Dx_CKF(1,7:9)',10); round(Dx_CKF(1,10:end)',10)*scale]

hours = time/3600;

% plots
% figure(10)   % s/c state
% marker = 'b-';
% ax1 = subplot(6,1,1);
% plot(hours,Dx_CKF(:,1),marker)
% ylabel('\Deltax (km)')
% grid on
% ax2 = subplot(6,1,2);
% plot(hours,Dx_CKF(:,2),marker)
% ylabel('\Deltay (km)')
% grid on
% ax3 = subplot(6,1,3);
% plot(hours,Dx_CKF(:,3),marker)
% ylabel('\Deltaz (km)')
% grid on
% ax4 = subplot(6,1,4);
% plot(hours,Dx_CKF(:,4)*scale,marker)
% ylabel('\Deltaxdot (m/s)')
% grid on
% ax5 = subplot(6,1,5);
% plot(hours,Dx_CKF(:,5)*scale,marker)
% ylabel('\Deltaydot (m/s)')
% grid on
% ax6 = subplot(6,1,6);
% plot(hours,Dx_CKF(:,6)*scale,marker)
% ylabel('\Deltazdot (m/s)')
% grid on
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
% xlim([hours(1) hours(end)])
% xlabel('Time (hours)')
% sgtitle('\Deltax_{CKF}, Spacecraft State')


% build separate arrays for each station's residuals
% to make plotting by color easier
ri101 = NaN(2,length(ri));
ri337 = NaN(2,length(ri));
ri394 = NaN(2,length(ri));
prefits101 = NaN(2,length(ri));
prefits337 = NaN(2,length(ri));
prefits394 = NaN(2,length(ri));
postfits101 = NaN(2,length(ri));
postfits337 = NaN(2,length(ri));
postfits394 = NaN(2,length(ri));

scale = 1000; % scale residuals (from km)
if (METER == 0) && (scale == 1000)
    label = '';
elseif (METER == 0) && (scale == 1)
    label = 'k';
elseif (METER == 1) && (scale == 1)
    label = '';
end

for kk = 1:length(ri)
    if IDs(kk) == 101
        postfits101(:,kk) = postfits(:,kk)*scale;
        prefits101(:,kk) = postfits(:,kk)*scale;
        ri101(:,kk) = ri(:,kk)*scale;
    elseif IDs(kk) == 337
        postfits337(:,kk) = postfits(:,kk)*scale;
        prefits337(:,kk) = postfits(:,kk)*scale;
        ri337(:,kk) = ri(:,kk)*scale;
    else
        postfits394(:,kk) = postfits(:,kk)*scale;
        prefits394(:,kk) = postfits(:,kk)*scale;
        ri394(:,kk) = ri(:,kk)*scale;
    end
end

% measurement residuals (yi - yhati)
% figure(6)
% subplot(2,1,1)
% plot(hours,ri101(1,:),'ro')
% hold on
% plot(hours,ri337(1,:),'bo')
% plot(hours,ri394(1,:),'go')
% hold off
% grid on
% xlim([hours(1) hours(end)])
% xlabel('Time (hours)')
% ylabel(['Rho Residual (',label,'m)'])
% legend('Station 101', 'Station 337','Station 394','location','southeast')
% subplot(2,1,2)
% plot(hours,ri101(2,:),'ro')
% hold on
% plot(hours,ri337(2,:),'bo')
% plot(hours,ri394(2,:),'go')
% hold off
% grid on
% xlim([hours(1) hours(end)])
% %ylim([-1e-4 1e-4])
% xlabel('Time (hours)')
% ylabel(['dRho Residual (',label,'m/s)'])
% sgtitle('Meaurement Residuals, CKF')

% pre-fit measurement residuals
figure(7)
subplot(2,1,1)
plot(hours,prefits101(1,:),'ro')
hold on
plot(hours,prefits337(1,:),'bo')
plot(hours,prefits394(1,:),'go')
hold off
grid on
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
ylabel(['Rho Residual (',label,'m)'])
legend('Station 101', 'Station 337','Station 394','location','southeast')
subplot(2,1,2)
plot(hours,prefits101(2,:),'ro')
hold on
plot(hours,prefits337(2,:),'bo')
plot(hours,prefits394(2,:),'go')
hold off
grid on
xlim([hours(1) hours(end)])
%ylim([-1e-4 1e-4])
xlabel('Time (hours)')
ylabel(['dRho Residual (',label,'m/s)'])
sgtitle('Pre-Fit Meaurement Residuals, CKF')

% post-fit measurement residuals
figure(8)
subplot(2,1,1)
plot(hours,postfits101(1,:),'ro')
hold on
plot(hours,postfits337(1,:),'bo')
plot(hours,postfits394(1,:),'go')
hold off
grid on
xlim([hours(1) hours(end)])
xlabel('Time (hours)')
ylabel(['Rho Residual (',label,'m)'])
legend('Station 101', 'Station 337','Station 394','location','north')
subplot(2,1,2)
plot(hours,postfits101(2,:),'ro')
hold on
plot(hours,postfits337(2,:),'bo')
plot(hours,postfits394(2,:),'go')
hold off
grid on
xlim([hours(1) hours(end)])
%ylim([-1e-4 1e-4])
xlabel('Time (hours)')
ylabel(['dRho Residual (',label,'m/s)'])
sgtitle(['Post-Fit Meaurement Residuals, CKF, ',num2str(maxMajor),' Iteration(s)'])

% RMS values of residuals
% remove NaNs
ri_rho = postfits(1,:);
ri_rho = ri_rho(~isnan(ri_rho));
ri_drho = postfits(2,:);
ri_drho = ri_drho(~isnan(ri_drho));

RMS_rho = sqrt((sum(ri_rho(1,:).^2,'all'))/length(ri_rho))*scale*100;     % cm
RMS_drho = sqrt((sum(ri_drho(1,:).^2,'all'))/length(ri_drho))*scale*1000; % mm/s

RMS = round([RMS_rho; RMS_drho],8)
%%
% compute traces of position and velocity covariances
tr_pos = zeros(1,size(P_plus,3));
tr_vel = zeros(1,size(P_plus,3));
for mm = 1:size(P_plus,3)
    tr_pos(mm) = trace(P_plus(1:3,1:3,mm));
    tr_vel(mm) = trace(P_plus(4:6,4:6,mm));
end

figure(300)
subplot(2,1,1)
semilogy(hours,tr_pos)
grid on
xlabel('Time (hours)')
ylabel('tr(P_{rr}) (m^2)')
subplot(2,1,2)
semilogy(hours,tr_vel)
grid minor
xlabel('Time (hours)')
ylabel('tr(P_{vv}) (m/s)^2')
sgtitle('Traces of Position and Velocity Covariances, CKF')

%%

% covariance ellipseoids at the final time
figure(200)
ellipsoid(0,0,0,3*sqrt(P_plus(1,1,ii)),3*sqrt(P_plus(2,2,ii)),3*sqrt(P_plus(3,3,ii)))
pbaspect([1 1 1])
xlabel('r_x (m)')
ylabel('r_y (m)')
zlabel('r_z (m)')
title('Covariance Ellipsoid, S/C Position, CKF')

figure(201)
ellipsoid(0,0,0,3*sqrt(P_plus(4,4,ii)),3*sqrt(P_plus(5,5,ii)),3*sqrt(P_plus(6,6,ii)))
pbaspect([1 1 1])
xlabel('v_x (m/s)')
ylabel('v_y (m/s)')
zlabel('v_z (m/s)')
title('Covariance Ellipsoid, S/C Velocity, CKF')
%%
LWsig = 1.5;
pos = 50;
vel = 1e-1;

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
sgtitle('Spacecraft State 3\sigma Covariance Bounds, CKF')