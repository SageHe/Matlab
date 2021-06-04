clear all;clc
%% Start of ASEN 6080 Project 1 main script
%{
Author: Sage Herrin
Due Date: 2-16-2021

This main script utilizes both a non-linear batch filet and a CKF to
estimate the position and velocity of a spacecraft, constants including
mu, J2, and drag coefficient Cd, as well as three observation station positions 
%}
% Run three iterations of the Batch filter first
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
meas = load('observations.txt');

rawmeas = meas;

meas = formatmeas(meas);
meas(:,3:4) = meas(:,3:4)/1000; % converting range/range rate measurements into km and km/s

t = meas(:,1);

J2 = 1.082626925638815e-3;
mu_e = 3.986004415e5; %km^3/s^2
CD = 2;

theta0 = 0; % check if theta0 is initialized at 0 (seems so from p. 2, "second past epoch")
stat101pos_ecef = [-5127.5100 -3794.1600 0.0];
stat337pos_ecef = [3860.9100 3238.4900 3898.0940];
stat394pos_ecef = [549.5050 -1380.8720 6182.1970];
statspos_ecef = [stat101pos_ecef;stat337pos_ecef;stat394pos_ecef];

Xs_ECI = stats_ECI(statspos_ecef,0,theta0);

P0 = diag([1e-6 1e-6 1e-6 1e-10 1e-10 1e-10 1e2 1e6 1e6 1e-16 1e-16 1e-16 1 1 1 1 1 1]);
R = 1e-3*diag([1e-10 1e-12]);

x0 = [757.700 5222.607 4851.5000 2.21321 4.67834 -5.37130 mu_e J2 CD Xs_ECI(1,1:3) Xs_ECI(2,1:3) Xs_ECI(3,1:3)]';
n = size(x0,1);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];

dt = t(2) - t(1);
%propagate initial state
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
X_true = x(:,1:18);


%initialize filter
j = 1;
xhat0(:,j) =  x0;
Deltax_minus(:,1) = zeros(n,1);
Deltax_plus(:,1) = ones(n,1);
while j < 4 %norm(Deltax_plus(:,end)) > 1e-7
ri = NaN(numel(t),2);
% xhat(:,1) = xhat0(:,j);
Lambda(:,:,j) = pinv(P0);
N = pinv(P0)*Deltax_minus;%P0\Deltax_minus; %better alternative to inv(P0)*Deltax_minus
%propagate for entire measurement duration for each major iteration
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);
for i = 2:numel(t)
    %read next observation
    yi = meas(i,3:4);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %integrate from t0 to ti
%     Phi = eye(n);
%     Phi_flat = reshape(Phi,n^2,1);
%     Z = [xhat0;Phi_flat];
%     tspan = [0:dt:t(i)];
%     [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    %determine station that produced meas
%     statnum = find(~isnan(yi));
%     if ~isempty(statnum)
%         statnum = statnum(2)/2;
%     end
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if meas(i,2) == 101%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(1,:));
            Hi(:,10:12) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(1,:));
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(1,:));
        elseif meas(i,2) == 337%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(2,:));
            Hi(:,13:15) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(2,:));
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(3,:));
            Hi(:,16:18) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(3,:));
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(3,:));
        end
        %accumulate observation
        
        ri(i,:) = yi - hi(i,:);
        if isnan(hi(i,:))
            ri(i,:) = zeros(1,2);
        end
        Phi_flat_ti_t0 = x(i,n+1:end);
        Phi_ti_t0 = reshape(Phi_flat_ti_t0,n,n);
        Lambda(:,:,i) = Lambda(:,:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*Hi*Phi_ti_t0;
        N(:,i) = N(:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*ri(i,:)';
    else    %no measurements occured during the timestep
        Lambda(:,:,i) = Lambda(:,:,i-1);
        N(:,i) = N(:,i-1);
    end
end
%solve normal equations
Deltax_plus(:,j) = pinv(Lambda(:,:,i))*N(:,i);
P0_plus(:,:,j) = pinv(Lambda(:,:,i));

xhat0(:,j+1) = xhat0(:,j) + Deltax_plus(:,j);
Deltax_minus(:,j+1) = Deltax_minus(:,j) - Deltax_plus(:,j);
j = j + 1;
end
    
for i = 1:numel(t)
    P_plus(:,:,i) = pinv(Lambda(:,:,i)); %might need to implement square root method here due to singular matrix
end

%propagate final state estimate
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);

riplot = NaN(numel(t),6);
for i = 2:numel(t)
    if meas(i,2) == 101
        riplot(i,1:2) = ri(i,:);
    elseif meas(i,2) == 337
        riplot(i,3:4) = ri(i,:);
    elseif meas(i,2) == 394
        riplot(i,5:6) = ri(i,:);
    end
end

figure
subplot(2,1,1)
hold on
plot(t,riplot(:,1),'*')
plot(t,riplot(:,3),'*')
plot(t,riplot(:,5),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Residual (km)')
title('Range Residaul VS Time')
legend('Station 101','Station 337','Station 394')
subplot(2,1,2)
hold on
plot(t,riplot(:,2),'*')
plot(t,riplot(:,4),'*')
plot(t,riplot(:,6),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
title('Range Rate Residuals VS Time')
legend('Station 101','Station 337','Station 394')
%calculate covariance ellipses 
errell = diag(P0_plus(:,:,end));
errell = 3*sqrt(errell);
figure
ellipsoid(0,0,0,errell(1),errell(2),errell(3))
axis equal
xlabel('S/C X Position')
ylabel('S/C Y Position')
zlabel('S/C Z Position')
title('Spacecraft Position Covariance Ellipse, 3\sigma,Batch')

figure
ellipsoid(0,0,0,errell(4),errell(5),errell(6))
axis equal
xlabel('S/C X Velocity')
ylabel('S/C Y Velocity')
zlabel('S/C Z Velocity')
title('Spacecraft Velocity Covariance Ellipse, 3\sigma,Batch')

%calculate presented residual RMS values
rangeRMSres = ri(~isnan(ri(:,1)));
rangeres_RMS = sqrt((1/size(rangeRMSres,1))*sum(rangeRMSres.^2));
rrRMSres = ri(:,2);
rrRMSres = rrRMSres(~isnan(rrRMSres));
rrres_RMS = sqrt((1/size(rrRMSres,1))*sum(rrRMSres.^2));

%plot estimated states and final estimated state errors for batch filter
%position state
figure
subplot(3,1,1)
plot(t,xhat(:,1))
xlabel('Time (s)')
ylabel('r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,2))
xlabel('Time (s)')
ylabel('r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,3))
xlabel('Time (s)')
ylabel('r_z (km)')
grid on
grid minor
sgtitle('S/C Position State VS Time,Batch')
%velocity state
figure
subplot(3,1,1)
plot(t,xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km)')
grid on
grid minor
sgtitle('S/C Velocity State VS Time,Batch')
% constants
figure
subplot(3,1,1)
plot(t,xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D estimates,Batch')
%station positions
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Estimates VS Time,Batch')
%position state errors
figure
subplot(3,1,1)
plot(t,X_true(:,1)-xhat(:,1))
xlabel('Time (s)')
ylabel('\Delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,2)-xhat(:,2))
xlabel('Time (s)')
ylabel('\Delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,3)-xhat(:,3))
xlabel('Time (s)')
ylabel('\Delta r_z (km)')
grid on
grid minor
sgtitle('S/C Position State Error VS Time,Batch')
%velocity state errors
figure
subplot(3,1,1)
plot(t,X_true(:,4)-xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,5)-xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,6)-xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km/s)')
grid on
grid minor
sgtitle('S/C Velocity State Error VS Time,Batch')
% constants state errors
figure
subplot(3,1,1)
plot(t,X_true(:,7)-xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,8)-xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,9)-xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D Estimate Errors,Batch')
%station positions errors
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(X_true(i,10:12)-xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(X_true(i,13:15)-xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(X_true(i,16:18)-xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Error VS Time, Batch')

%% This section of the main script will now utilize a CKF to estimate the desired states
clear;clc

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
meas = load('observations.txt');
rawmeas = meas;

meas = formatmeas(meas);
meas(:,3:4) = meas(:,3:4)/1000; % converting range/range rate measurements into km and km/s

t = meas(:,1);

J2 = 1.082626925638815e-3;
mu_e = 3.986004415e5; %km^3/s^2
CD = 2;

theta0 = 0; % check if theta0 is initialized at 0 (seems so from p. 2, "second past epoch")
stat101pos_ecef = [-5127.5100 -3794.1600 0.0];
stat337pos_ecef = [3860.9100 3238.4900 3898.0940];
stat394pos_ecef = [549.5050 -1380.8720 6182.1970];
statspos_ecef = [stat101pos_ecef;stat337pos_ecef;stat394pos_ecef];

Xs_ECI = stats_ECI(statspos_ecef,0,theta0);

P0 = diag([1e-6 1e-6 1e-6 1e-10 1e-10 1e-10 1e2 1e6 1 1e-16 1e-16 1e-16 1 1 1 1 1 1]);
R = 1e-3*diag([1e-10 1e-12]);

x0 = [757.700 5222.607 4851.5000 2.21321 4.67834 -5.37130 mu_e J2 CD Xs_ECI(1,1:3) Xs_ECI(2,1:3) Xs_ECI(3,1:3)]';
dx = [.1 -.1 .1 0.01 -0.01 0.01 zeros(1,12)]';

x0 = x0 + dx;
n = size(x0,1);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];

dt = t(2) - t(1);
tspan = [0:dt:t(end)];
[~,X_true] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);

j = 1;
x_hat0(1,:) = x0;
x_hat0 = x_hat0';
Deltax_minus = zeros(n,1);
Deltax_plus = Deltax_minus;
Deltax0_plus = ones(n,1);
while j < 3
%Initialize Filter
ti_minus = t(1,1);
x_t_im1 = x_hat0(:,j);
x_hat(:,1) = x_t_im1;
P_plus(:,:,1) = P0; 
ri = NaN(numel(t),2);
ri_pre = ri;
ri_post = ri;
% Deltax_i_m1_plus = zeros(7,1);
% numstats = size(rangemeas,2);
% y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
% y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
% Read next observation
for i = 2:numel(t)
    ti = t(i);
%     statnum = find(~isnan(parsedrange(i,2:4)));
    yi = meas(i,3:4);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
%     n = size(x0,1);
    %Propagate from t_i-1 to t_i
%     dt = ti - ti_minus;
    Phi = eye(n);
    Phi_flat = reshape(Phi,size(Phi,1)^2,1);
    Z = [x_t_im1;Phi_flat];
    tspan = [ti_minus:dt/2:ti];
    [~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
    x_ti = x(end,1:n);
    %Time Update
    Phi = reshape(x(end,n+1:end),n,n);
    Deltax_minus(:,i) = Phi*Deltax_plus(:,i-1);
    P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi';
    %Process Observations
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if meas(i,2) == 101%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:9) = Htilde_sc_proj(x_ti(1,1:6),Xs_ECI(1,:));
            Hi(:,10:12) = Htilde_obs_proj(x_ti(1,1:6),Xs_ECI(1,:));
           [hi(i,:)] = predictmeas(x_ti(1,1:6),Xs_ECI(1,:));
        elseif meas(i,2) == 337%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:9) = Htilde_sc_proj(x_ti(1,1:6),Xs_ECI(2,:));
            Hi(:,13:15) = Htilde_obs_proj(x_ti(1,1:6),Xs_ECI(2,:));
           [hi(i,:)] = predictmeas(x_ti(1,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:9) = Htilde_sc_proj(x_ti(1,1:6),Xs_ECI(3,:));
            Hi(:,16:18) = Htilde_obs_proj(x_ti(1,1:6),Xs_ECI(3,:));
           [hi(i,:)] = predictmeas(x_ti(1,1:6),Xs_ECI(3,:));
        end
        ri(i,:) = yi - hi(i,:);
        Ki = P_minus(:,:,i)*Hi'*inv(Hi*P_minus(:,:,i)*Hi' + R);
        % measurement update
        Deltax_plus(:,i) = Deltax_minus(:,i) + Ki*(ri(i,:)' - Hi*Deltax_minus(:,i));
        P_plus(:,:,i) = (eye(n) - Ki*Hi)*P_minus(:,:,i)*(eye(n) - Ki*Hi)' + Ki*R*Ki';
    else
        hi(i,:) = NaN(1,2);
        Deltax_plus(:,i) = Deltax_minus(:,i);
        P_plus(:,:,i) = P_minus(:,:,i);
    end
%     [hi(i,:)] = predictmeas(x_ti(1:6),Xs_ECI);
%     ri(i,:) = yi - hi(i,:);
    
%     Hi = zeros(2*numstats,n);
%     for k = 1:numstats
%         if ~isnan(ri(i,2*k))
%             Hi(2*k-1:2*k,1:6) = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(k,:));
%         elseif isnan(ri(i,2*k))
%             ri(i,2*k-1:2*k) = [0 0];
%         end
% %         Hi = Htilde_sc_rho_rhod(x(end,1:6),[Rs_ECI Vs_ECI]);
%     end
%     Ki = P_minus*Hi'*inv(Hi*P_minus*Hi' + R);
%     %Measurement update
%     Deltax_plus(:,i) = Deltax_minus(:,i) + Ki*(ri(i,:)' - Hi*Deltax_minus(:,i));
%     P_plus = (eye(7) - Ki*Hi)*P_minus*(eye(7) - Ki*Hi)' + Ki*R*Ki';
    %set up for next minor iteration
    ti_minus = ti;
    x_t_im1 = x_ti';
    x_hat(:,i) = x_ti' + Deltax_plus(:,i);
    %calculate pre-fit and post-fit residuals
    ri_pre(i,:) = yi - hi(i,:) - (Hi*Deltax_minus(:,i))';
    ri_post(i,:) = yi - hi(i,:) - (Hi*Deltax_plus(:,i))';
%     Pi_mi = P_plus;
%     Deltax_i_m1_plus = Deltaxi_plus;
    
    
end
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x_ti';Phi_flat];
tspan = [ti 0];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);

Phi = reshape(x(end,n+1:end),n,n);
Deltax0_plus = (Phi)*Deltax_plus(:,i);
x_hat0(:,j+1) = x_hat0(:,j) + Deltax0_plus;
j = j + 1;
end
ri_pre_plot = NaN(numel(t),6);
ri_post_plot = NaN(numel(t),6);
for i = 2:numel(t)
    if meas(i,2) == 101
        ri_pre_plot(i,1:2) = ri_pre(i,:);
        ri_post_plot(i,1:2) = ri_post(i,:);
    elseif meas(i,2) == 337
        ri_pre_plot(i,3:4) = ri_pre(i,:);
        ri_post_plot(i,3:4) = ri_post(i,:);
    elseif meas(i,2) == 394
        ri_pre_plot(i,5:6) = ri_pre(i,:);
        ri_post_plot(i,5:6) = ri_post(i,:);
    end
end
x_hat = x_hat';
figure
subplot(2,1,1)
hold on
plot(t,ri_pre_plot(:,1),'*')
plot(t,ri_pre_plot(:,3),'*')
plot(t,ri_pre_plot(:,5),'*')
xlabel('Time (s)')
ylabel('Range Residuals (km)')
grid on
grid minor
legend('Station 101','Station 337','Station 394')
subplot(2,1,2)
hold on
plot(t,ri_pre_plot(:,2),'*')
plot(t,ri_pre_plot(:,4),'*')
plot(t,ri_pre_plot(:,6),'*')
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
grid on 
grid minor
legend('Station 101','Station 337','Station 394')
sgtitle('Pre-fit CKF Measurement Residuals VS Time')

figure
subplot(2,1,1)
hold on
plot(t,ri_post_plot(:,1),'*')
plot(t,ri_post_plot(:,3),'*')
plot(t,ri_post_plot(:,5),'*')
xlabel('Time (s)')
ylabel('Range Residuals (km)')
grid on
grid minor
legend('Station 101','Station 337','Station 394')
subplot(2,1,2)
hold on
plot(t,ri_post_plot(:,2),'*')
plot(t,ri_post_plot(:,4),'*')
plot(t,ri_post_plot(:,6),'*')
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
grid on 
grid minor
legend('Station 101','Station 337','Station 394')
sgtitle('Post-fit CKF Measurement Residuals VS Time')

%calculate residual RMS values
rangeRMSres = ri_post(~isnan(ri_post(:,1)));
rangeres_RMS = sqrt((1/size(rangeRMSres,1))*sum(rangeRMSres.^2));
rrRMSres = ri_post(:,2);
rrRMSres = rrRMSres(~isnan(rrRMSres));
rrres_RMS = sqrt((1/size(rrRMSres,1))*sum(rrRMSres.^2));
%calculate trace of sequential filter pos. and vel. covariances
for i = 1:size(P_plus,3)
    pos_trace(i) = trace(P_plus(1:3,1:3,i));
    vel_trace(i) = trace(P_plus(4:6,4:6,i));
end
figure
subplot(2,1,1)
semilogy(t,pos_trace)
xlabel('Time (s)')
ylabel('Position Covariance Trace (m^2)')
grid on
grid minor
subplot(2,1,2)
semilogy(t,vel_trace)
xlabel('Time (s)')
ylabel('Velocity Covariance Trace (m/s)^2')
grid on
grid minor
sgtitle('CKF Position and Velocity Covariance Traces VS Time')

%calculate covariance ellipses 
errell = diag(P_plus(:,:,end));
errell = 3*sqrt(errell);
figure
ellipsoid(0,0,0,errell(1),errell(2),errell(3))
axis equal
xlabel('S/C X Position')
ylabel('S/C Y Position')
zlabel('S/C Z Position')
title('Spacecraft Position Covariance Ellipse, 3\sigma,CKF')

figure
ellipsoid(0,0,0,errell(4),errell(5),errell(6))
axis equal
xlabel('S/C X Velocity')
ylabel('S/C Y Velocity')
zlabel('S/C Z Velocity')
title('Spacecraft Velocity Covariance Ellipse, 3\sigma,CKF')

%plot estimated states and final estimated state errors for CKF 
%position state
figure
subplot(3,1,1)
plot(t,x_hat(:,1))
xlabel('Time (s)')
ylabel('r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,x_hat(:,2))
xlabel('Time (s)')
ylabel('r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,x_hat(:,3))
xlabel('Time (s)')
ylabel('r_z (km)')
grid on
grid minor
sgtitle('S/C Position State VS Time,CKF')
%velocity state
figure
subplot(3,1,1)
plot(t,x_hat(:,4))
xlabel('Time (s)')
ylabel('v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(t,x_hat(:,5))
xlabel('Time (s)')
ylabel('v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(t,x_hat(:,6))
xlabel('Time (s)')
ylabel('v_z (km/s)')
grid on
grid minor
sgtitle('S/C Velocity State VS Time,CKF')
% constants
figure
subplot(3,1,1)
plot(t,x_hat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,x_hat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,x_hat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D estimates,CKF')
%station positions
for i = 1:size(x_hat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(x_hat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(x_hat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(x_hat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Estimates VS Time,CKF')
%position state errors
figure
subplot(3,1,1)
plot(t,X_true(:,1)-x_hat(:,1))
xlabel('Time (s)')
ylabel('\Delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,2)-x_hat(:,2))
xlabel('Time (s)')
ylabel('\Delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,3)-x_hat(:,3))
xlabel('Time (s)')
ylabel('\Delta r_z (km)')
grid on
grid minor
sgtitle('S/C Position State Error VS Time,CKF')
%velocity state errors
figure
subplot(3,1,1)
plot(t,X_true(:,4)-x_hat(:,4))
xlabel('Time (s)')
ylabel('v_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,5)-x_hat(:,5))
xlabel('Time (s)')
ylabel('v_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,6)-x_hat(:,6))
xlabel('Time (s)')
ylabel('v_z (km)')
grid on
grid minor
sgtitle('S/C Velocity State Error VS Time,CKF')
% constants state errors
figure
subplot(3,1,1)
plot(t,X_true(:,7)-x_hat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,8)-x_hat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,9)-x_hat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D Estimate Errors,CKF')
%station positions errors
for i = 1:size(x_hat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(X_true(i,10:12)-x_hat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(X_true(i,13:15)-x_hat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(X_true(i,16:18)-x_hat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Error VS Time, CKF')

%% This section explores the strengths of the range and range rate data by generating solutions for each data type alone using the Batch filter, starting with range measurements
clear all;clc
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
meas = load('observations.txt');

rawmeas = meas;

meas = formatmeas(meas);
meas(:,3:4) = meas(:,3:4)/1000; % converting range/range rate measurements into km and km/s

t = meas(:,1);

J2 = 1.082626925638815e-3;
mu_e = 3.986004415e5; %km^3/s^2
CD = 2;

theta0 = 0; % check if theta0 is initialized at 0 (seems so from p. 2, "second past epoch")
stat101pos_ecef = [-5127.5100 -3794.1600 0.0];
stat337pos_ecef = [3860.9100 3238.4900 3898.0940];
stat394pos_ecef = [549.5050 -1380.8720 6182.1970];
statspos_ecef = [stat101pos_ecef;stat337pos_ecef;stat394pos_ecef];

Xs_ECI = stats_ECI(statspos_ecef,0,theta0);

P0 = diag([1 1 1 1 1 1 1e2 1e6 1e6 1e-16 1e-16 1e-16 1 1 1 1 1 1]);
% R = diag([1e-10 1e-12]);
R = 1e-10;

x0 = [757.700 5222.607 4851.5000 2.21321 4.67834 -5.37130 mu_e J2 CD Xs_ECI(1,1:3) Xs_ECI(2,1:3) Xs_ECI(3,1:3)]';
n = size(x0,1);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];

dt = t(2) - t(1);
%propagate initial state
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
X_true = x(:,1:18);


%initialize filter
j = 1;
xhat0(:,j) =  x0;
Deltax_minus(:,1) = zeros(n,1);
Deltax_plus(:,1) = ones(n,1);
while j < 4 %norm(Deltax_plus(:,end)) > 1e-7
ri = NaN(numel(t),1);
% xhat(:,1) = xhat0(:,j);
Lambda(:,:,j) = pinv(P0);
N = pinv(P0)*Deltax_minus;%P0\Deltax_minus; %better alternative to inv(P0)*Deltax_minus
%propagate for entire measurement duration for each major iteration
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);
for i = 2:numel(t)
    %read next observation
    yi = meas(i,3);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %integrate from t0 to ti
%     Phi = eye(n);
%     Phi_flat = reshape(Phi,n^2,1);
%     Z = [xhat0;Phi_flat];
%     tspan = [0:dt:t(i)];
%     [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    %determine station that produced meas
%     statnum = find(~isnan(yi));
%     if ~isempty(statnum)
%         statnum = statnum(2)/2;
%     end
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if meas(i,2) == 101%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(1,:));
            Hi(:,10:12) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(1,:));
            Hi(2,:) = [];
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(1,:));
        elseif meas(i,2) == 337%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(2,:));
            Hi(:,13:15) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(2,:));
            Hi(2,:) = [];
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(3,:));
            Hi(:,16:18) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(3,:));
            Hi(2,:) = [];
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(3,:));
        end
        %accumulate observation
        
        ri(i,:) = yi - hi(i,1);
        if isnan(hi(i,:))
            ri(i,:) = zeros(1,2);
        end
        Phi_flat_ti_t0 = x(i,n+1:end);
        Phi_ti_t0 = reshape(Phi_flat_ti_t0,n,n);
        Lambda(:,:,i) = Lambda(:,:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*Hi*Phi_ti_t0;
        N(:,i) = N(:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*ri(i,:)';
    else    %no measurements occured during the timestep
        Lambda(:,:,i) = Lambda(:,:,i-1);
        N(:,i) = N(:,i-1);
    end
end
%solve normal equations
Deltax_plus(:,j) = pinv(Lambda(:,:,i))*N(:,i);
P0_plus(:,:,j) = pinv(Lambda(:,:,i));

xhat0(:,j+1) = xhat0(:,j) + Deltax_plus(:,j);
Deltax_minus(:,j+1) = Deltax_minus(:,j) - Deltax_plus(:,j);
j = j + 1;
end
    
for i = 1:numel(t)
    P_plus(:,:,i) = pinv(Lambda(:,:,i)); %might need to implement square root method here due to singular matrix
end

%propagate final state estimate
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);

riplot = NaN(numel(t),6);
for i = 2:numel(t)
    if meas(i,2) == 101
        riplot(i,1:2) = ri(i,:);
    elseif meas(i,2) == 337
        riplot(i,3:4) = ri(i,:);
    elseif meas(i,2) == 394
        riplot(i,5:6) = ri(i,:);
    end
end

figure
hold on
plot(t,riplot(:,1),'*')
plot(t,riplot(:,3),'*')
plot(t,riplot(:,5),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Residual (km)')
title('Range Residaul VS Time')
legend('Station 101','Station 337','Station 394')

figure
hold on
plot(t,riplot(:,2),'*')
plot(t,riplot(:,4),'*')
plot(t,riplot(:,6),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
title('Range Rate Residuals VS Time')
legend('Station 101','Station 337','Station 394')
%calculate covariance ellipses 
errell = diag(P0_plus(:,:,end));
errell = 3*sqrt(errell);
figure
ellipsoid(0,0,0,errell(1),errell(2),errell(3))
axis equal
xlabel('S/C X Position')
ylabel('S/C Y Position')
zlabel('S/C Z Position')
title('Spacecraft Position Covariance Ellipse, 3\sigma,Batch,Range Only')

figure
ellipsoid(0,0,0,errell(4),errell(5),errell(6))
axis equal
xlabel('S/C X Velocity')
ylabel('S/C Y Velocity')
zlabel('S/C Z Velocity')
title('Spacecraft Velocity Covariance Ellipse, 3\sigma,Batch,Range Only')

%calculate presented residual RMS values
rangeRMSres = ri(~isnan(ri(:,1)));
rangeres_RMS = sqrt((1/size(rangeRMSres,1))*sum(rangeRMSres.^2));
% rrRMSres = ri(:,2);
% rrRMSres = rrRMSres(~isnan(rrRMSres));
% rrres_RMS = sqrt((1/size(rrRMSres,1))*sum(rrRMSres.^2));

%plot estimated states and final estimated state errors for batch filter
%position state
figure
subplot(3,1,1)
plot(t,xhat(:,1))
xlabel('Time (s)')
ylabel('r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,2))
xlabel('Time (s)')
ylabel('r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,3))
xlabel('Time (s)')
ylabel('r_z (km)')
grid on
grid minor
sgtitle('S/C Position State VS Time,Batch')
%velocity state
figure
subplot(3,1,1)
plot(t,xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km)')
grid on
grid minor
sgtitle('S/C Velocity State VS Time,Batch')
% constants
figure
subplot(3,1,1)
plot(t,xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D estimates,Batch')
%station positions
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Estimates VS Time,Batch')
%position state errors
figure
subplot(3,1,1)
plot(t,X_true(:,1)-xhat(:,1))
xlabel('Time (s)')
ylabel('\Delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,2)-xhat(:,2))
xlabel('Time (s)')
ylabel('\Delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,3)-xhat(:,3))
xlabel('Time (s)')
ylabel('\Delta r_z (km)')
grid on
grid minor
sgtitle('S/C Position State Error VS Time,Batch')
%velocity state errors
figure
subplot(3,1,1)
plot(t,X_true(:,4)-xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,5)-xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,6)-xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km/s)')
grid on
grid minor
sgtitle('S/C Velocity State Error VS Time,Batch')
% constants state errors
figure
subplot(3,1,1)
plot(t,X_true(:,7)-xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,8)-xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,9)-xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D Estimate Errors,Batch')
%station positions errors
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(X_true(i,10:12)-xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(X_true(i,13:15)-xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(X_true(i,16:18)-xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Error VS Time, Batch')
%% This section explores the strengths of the range and range rate data by generating solutions for each data type alone using the Batch filter, now only using range rate measurements
clear all;clc
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
meas = load('observations.txt');

rawmeas = meas;

meas = formatmeas(meas);
meas(:,3:4) = meas(:,3:4)/1000; % converting range/range rate measurements into km and km/s

t = meas(:,1);

J2 = 1.082626925638815e-3;
mu_e = 3.986004415e5; %km^3/s^2
CD = 2;

theta0 = 0; % check if theta0 is initialized at 0 (seems so from p. 2, "second past epoch")
stat101pos_ecef = [-5127.5100 -3794.1600 0.0];
stat337pos_ecef = [3860.9100 3238.4900 3898.0940];
stat394pos_ecef = [549.5050 -1380.8720 6182.1970];
statspos_ecef = [stat101pos_ecef;stat337pos_ecef;stat394pos_ecef];

Xs_ECI = stats_ECI(statspos_ecef,0,theta0);

P0 = diag([1 1 1 1 1 1 1e2 1e6 1e6 1e-16 1e-16 1e-16 1 1 1 1 1 1]);
% R = diag([1e-10 1e-12]);
R = 1e-12;

x0 = [757.700 5222.607 4851.5000 2.21321 4.67834 -5.37130 mu_e J2 CD Xs_ECI(1,1:3) Xs_ECI(2,1:3) Xs_ECI(3,1:3)]';
n = size(x0,1);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];

dt = t(2) - t(1);
%propagate initial state
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
X_true = x(:,1:18);


%initialize filter
j = 1;
xhat0(:,j) =  x0;
Deltax_minus(:,1) = zeros(n,1);
Deltax_plus(:,1) = ones(n,1);
while j < 4 %norm(Deltax_plus(:,end)) > 1e-7
ri = NaN(numel(t),1);
% xhat(:,1) = xhat0(:,j);
Lambda(:,:,j) = pinv(P0);
N = pinv(P0)*Deltax_minus;%P0\Deltax_minus; %better alternative to inv(P0)*Deltax_minus
%propagate for entire measurement duration for each major iteration
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);
for i = 2:numel(t)
    %read next observation
    yi = meas(i,4);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %integrate from t0 to ti
%     Phi = eye(n);
%     Phi_flat = reshape(Phi,n^2,1);
%     Z = [xhat0;Phi_flat];
%     tspan = [0:dt:t(i)];
%     [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    %determine station that produced meas
%     statnum = find(~isnan(yi));
%     if ~isempty(statnum)
%         statnum = statnum(2)/2;
%     end
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if meas(i,2) == 101%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(1,:));
            Hi(:,10:12) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(1,:));
            Hi(1,:) = [];
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(1,:));
        elseif meas(i,2) == 337%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(2,:));
            Hi(:,13:15) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(2,:));
            Hi(1,:) = [];
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(3,:));
            Hi(:,16:18) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(3,:));
            Hi(1,:) = [];
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(3,:));
        end
        %accumulate observation
        
        ri(i,:) = yi - hi(i,2);
        if isnan(hi(i,:))
            ri(i,:) = zeros(1,2);
        end
        Phi_flat_ti_t0 = x(i,n+1:end);
        Phi_ti_t0 = reshape(Phi_flat_ti_t0,n,n);
        Lambda(:,:,i) = Lambda(:,:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*Hi*Phi_ti_t0;
        N(:,i) = N(:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*ri(i,:)';
    else    %no measurements occured during the timestep
        Lambda(:,:,i) = Lambda(:,:,i-1);
        N(:,i) = N(:,i-1);
    end
end
%solve normal equations
Deltax_plus(:,j) = pinv(Lambda(:,:,i))*N(:,i);
P0_plus(:,:,j) = pinv(Lambda(:,:,i));

xhat0(:,j+1) = xhat0(:,j) + Deltax_plus(:,j);
Deltax_minus(:,j+1) = Deltax_minus(:,j) - Deltax_plus(:,j);
j = j + 1;
end
    
for i = 1:numel(t)
    P_plus(:,:,i) = pinv(Lambda(:,:,i)); %might need to implement square root method here due to singular matrix
end

%propagate final state estimate
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);

riplot = NaN(numel(t),6);
for i = 2:numel(t)
    if meas(i,2) == 101
        riplot(i,1:2) = ri(i,:);
    elseif meas(i,2) == 337
        riplot(i,3:4) = ri(i,:);
    elseif meas(i,2) == 394
        riplot(i,5:6) = ri(i,:);
    end
end

figure
hold on
plot(t,riplot(:,1),'*')
plot(t,riplot(:,3),'*')
plot(t,riplot(:,5),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Residual (km)')
title('Range Residaul VS Time')
legend('Station 101','Station 337','Station 394')

figure
hold on
plot(t,riplot(:,2),'*')
plot(t,riplot(:,4),'*')
plot(t,riplot(:,6),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
title('Range Rate Residuals VS Time')
legend('Station 101','Station 337','Station 394')
%calculate covariance ellipses 
errell = diag(P0_plus(:,:,end));
errell = 3*sqrt(errell);
figure
ellipsoid(0,0,0,errell(1),errell(2),errell(3))
axis equal
xlabel('S/C X Position')
ylabel('S/C Y Position')
zlabel('S/C Z Position')
title('Spacecraft Position Covariance Ellipse, 3\sigma,Batch,Range Rate Only')

figure
ellipsoid(0,0,0,errell(4),errell(5),errell(6))
axis equal
xlabel('S/C X Velocity')
ylabel('S/C Y Velocity')
zlabel('S/C Z Velocity')
title('Spacecraft Velocity Covariance Ellipse, 3\sigma,Batch,Range Rate Only')

%calculate presented residual RMS values
rr_RMSres = ri(~isnan(ri(:,1)));
rr_res_RMS = sqrt((1/size(rr_RMSres,1))*sum(rr_RMSres.^2));
% rrRMSres = ri(:,2);
% rrRMSres = rrRMSres(~isnan(rrRMSres));
% rrres_RMS = sqrt((1/size(rrRMSres,1))*sum(rrRMSres.^2));

%plot estimated states and final estimated state errors for batch filter
%position state
figure
subplot(3,1,1)
plot(t,xhat(:,1))
xlabel('Time (s)')
ylabel('r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,2))
xlabel('Time (s)')
ylabel('r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,3))
xlabel('Time (s)')
ylabel('r_z (km)')
grid on
grid minor
sgtitle('S/C Position State VS Time,Batch')
%velocity state
figure
subplot(3,1,1)
plot(t,xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km)')
grid on
grid minor
sgtitle('S/C Velocity State VS Time,Batch')
% constants
figure
subplot(3,1,1)
plot(t,xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D estimates,Batch')
%station positions
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Estimates VS Time,Batch')
%position state errors
figure
subplot(3,1,1)
plot(t,X_true(:,1)-xhat(:,1))
xlabel('Time (s)')
ylabel('\Delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,2)-xhat(:,2))
xlabel('Time (s)')
ylabel('\Delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,3)-xhat(:,3))
xlabel('Time (s)')
ylabel('\Delta r_z (km)')
grid on
grid minor
sgtitle('S/C Position State Error VS Time,Batch')
%velocity state errors
figure
subplot(3,1,1)
plot(t,X_true(:,4)-xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,5)-xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,6)-xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km/s)')
grid on
grid minor
sgtitle('S/C Velocity State Error VS Time,Batch')
% constants state errors
figure
subplot(3,1,1)
plot(t,X_true(:,7)-xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,8)-xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,9)-xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D Estimate Errors,Batch')
%station positions errors
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(X_true(i,10:12)-xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(X_true(i,13:15)-xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(X_true(i,16:18)-xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Error VS Time, Batch')

%% This section looks at the esitmatation of the states now with no fixed stations, starting with the batch filter
clear all;clc
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
meas = load('observations.txt');

rawmeas = meas;

meas = formatmeas(meas);
meas(:,3:4) = meas(:,3:4)/1000; % converting range/range rate measurements into km and km/s

t = meas(:,1);

J2 = 1.082626925638815e-3;
mu_e = 3.986004415e5; %km^3/s^2
CD = 2;

theta0 = 0; % check if theta0 is initialized at 0 (seems so from p. 2, "second past epoch")
stat101pos_ecef = [-5127.5100 -3794.1600 0.0];
stat337pos_ecef = [3860.9100 3238.4900 3898.0940];
stat394pos_ecef = [549.5050 -1380.8720 6182.1970];
statspos_ecef = [stat101pos_ecef;stat337pos_ecef;stat394pos_ecef];

Xs_ECI = stats_ECI(statspos_ecef,0,theta0);

P0 = diag([1 1 1 1 1 1 1e2 1e6 1e6 1 1 1 1 1 1 1e-16 1e-16 1e-16]);
R = diag([1e-10 1e-12]);

x0 = [757.700 5222.607 4851.5000 2.21321 4.67834 -5.37130 mu_e J2 CD Xs_ECI(1,1:3) Xs_ECI(2,1:3) Xs_ECI(3,1:3)]';
n = size(x0,1);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];

dt = t(2) - t(1);
%propagate initial state
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
X_true = x(:,1:18);


%initialize filter
j = 1;
xhat0(:,j) =  x0;
Deltax_minus(:,1) = zeros(n,1);
Deltax_plus(:,1) = ones(n,1);
while j < 4 %norm(Deltax_plus(:,end)) > 1e-7
ri = NaN(numel(t),2);
% xhat(:,1) = xhat0(:,j);
Lambda(:,:,j) = pinv(P0);
N = pinv(P0)*Deltax_minus;%P0\Deltax_minus; %better alternative to inv(P0)*Deltax_minus
%propagate for entire measurement duration for each major iteration
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);
for i = 2:numel(t)
    %read next observation
    yi = meas(i,3:4);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %integrate from t0 to ti
%     Phi = eye(n);
%     Phi_flat = reshape(Phi,n^2,1);
%     Z = [xhat0;Phi_flat];
%     tspan = [0:dt:t(i)];
%     [~,x] = ode45(@(t,Z) keplerJ2_wPhi_ODE(t,Z),tspan,Z,opts);
    %determine station that produced meas
%     statnum = find(~isnan(yi));
%     if ~isempty(statnum)
%         statnum = statnum(2)/2;
%     end
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if meas(i,2) == 101%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(1,:));
            Hi(:,10:12) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(1,:));
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(1,:));
        elseif meas(i,2) == 337%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(2,:));
            Hi(:,13:15) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(2,:));
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:9) = Htilde_sc_proj(xhat(i,1:6),Xs_ECI(3,:));
            Hi(:,16:18) = Htilde_obs_proj(xhat(i,1:6),Xs_ECI(3,:));
           [hi(i,:)] = predictmeas(xhat(i,1:6),Xs_ECI(3,:));
        end
        %accumulate observation
        
        ri(i,:) = yi - hi(i,:);
        if isnan(hi(i,:))
            ri(i,:) = zeros(1,2);
        end
        Phi_flat_ti_t0 = x(i,n+1:end);
        Phi_ti_t0 = reshape(Phi_flat_ti_t0,n,n);
        Lambda(:,:,i) = Lambda(:,:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*Hi*Phi_ti_t0;
        N(:,i) = N(:,i-1) + (Hi*Phi_ti_t0)'*pinv(R)*ri(i,:)';
    else    %no measurements occured during the timestep
        Lambda(:,:,i) = Lambda(:,:,i-1);
        N(:,i) = N(:,i-1);
    end
end
%solve normal equations
Deltax_plus(:,j) = pinv(Lambda(:,:,i))*N(:,i);
P0_plus(:,:,j) = pinv(Lambda(:,:,i));

xhat0(:,j+1) = xhat0(:,j) + Deltax_plus(:,j);
Deltax_minus(:,j+1) = Deltax_minus(:,j) - Deltax_plus(:,j);
j = j + 1;
end
    
for i = 1:numel(t)
    P_plus(:,:,i) = pinv(Lambda(:,:,i)); %might need to implement square root method here due to singular matrix
end

%propagate final state estimate
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [xhat0(:,j);Phi_flat];
tspan = [0:dt:t(end)];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
xhat = x(:,1:18);

riplot = NaN(numel(t),6);
for i = 2:numel(t)
    if meas(i,2) == 101
        riplot(i,1:2) = ri(i,:);
    elseif meas(i,2) == 337
        riplot(i,3:4) = ri(i,:);
    elseif meas(i,2) == 394
        riplot(i,5:6) = ri(i,:);
    end
end

figure
hold on
plot(t,riplot(:,1),'*')
plot(t,riplot(:,3),'*')
plot(t,riplot(:,5),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Residual (km)')
title('Range Residaul VS Time')
legend('Station 101','Station 337','Station 394')

figure
hold on
plot(t,riplot(:,2),'*')
plot(t,riplot(:,4),'*')
plot(t,riplot(:,6),'*')
grid on
grid minor
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
title('Range Rate Residuals VS Time')
legend('Station 101','Station 337','Station 394')
%calculate covariance ellipses 
errell = diag(P0_plus(:,:,end));
errell = 3*sqrt(errell);
figure
ellipsoid(0,0,0,errell(1),errell(2),errell(3))
axis equal
xlabel('S/C X Position')
ylabel('S/C Y Position')
zlabel('S/C Z Position')
title('Spacecraft Position Covariance Ellipse, 3\sigma,Batch')

figure
ellipsoid(0,0,0,errell(4),errell(5),errell(6))
axis equal
xlabel('S/C X Velocity')
ylabel('S/C Y Velocity')
zlabel('S/C Z Velocity')
title('Spacecraft Velocity Covariance Ellipse, 3\sigma,Batch')

%calculate presented residual RMS values
rangeRMSres = ri(~isnan(ri(:,1)));
rangeres_RMS = sqrt((1/size(rangeRMSres,1))*sum(rangeRMSres.^2));
rrRMSres = ri(:,2);
rrRMSres = rrRMSres(~isnan(rrRMSres));
rrres_RMS = sqrt((1/size(rrRMSres,1))*sum(rrRMSres.^2));

%plot estimated states and final estimated state errors for batch filter
%position state
figure
subplot(3,1,1)
plot(t,xhat(:,1))
xlabel('Time (s)')
ylabel('r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,2))
xlabel('Time (s)')
ylabel('r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,3))
xlabel('Time (s)')
ylabel('r_z (km)')
grid on
grid minor
sgtitle('S/C Position State VS Time,Batch')
%velocity state
figure
subplot(3,1,1)
plot(t,xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km)')
grid on
grid minor
sgtitle('S/C Velocity State VS Time,Batch')
% constants
figure
subplot(3,1,1)
plot(t,xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D estimates,Batch')
%station positions
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Estimates VS Time,Batch')
%position state errors
figure
subplot(3,1,1)
plot(t,X_true(:,1)-xhat(:,1))
xlabel('Time (s)')
ylabel('\Delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,2)-xhat(:,2))
xlabel('Time (s)')
ylabel('\Delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,3)-xhat(:,3))
xlabel('Time (s)')
ylabel('\Delta r_z (km)')
grid on
grid minor
sgtitle('S/C Position State Error VS Time,Batch')
%velocity state errors
figure
subplot(3,1,1)
plot(t,X_true(:,4)-xhat(:,4))
xlabel('Time (s)')
ylabel('v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,5)-xhat(:,5))
xlabel('Time (s)')
ylabel('v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,6)-xhat(:,6))
xlabel('Time (s)')
ylabel('v_z (km/s)')
grid on
grid minor
sgtitle('S/C Velocity State Error VS Time,Batch')
% constants state errors
figure
subplot(3,1,1)
plot(t,X_true(:,7)-xhat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,8)-xhat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,9)-xhat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D Estimate Errors,Batch')
%station positions errors
for i = 1:size(xhat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(X_true(i,10:12)-xhat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(X_true(i,13:15)-xhat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(X_true(i,16:18)-xhat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Error VS Time, Batch')
%% This section looks at the effect of not fixing any stations while using the CKF to estimate the states
clear;clc

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
meas = load('observations.txt');
rawmeas = meas;

meas = formatmeas(meas);
meas(:,3:4) = meas(:,3:4)/1000; % converting range/range rate measurements into km and km/s

t = meas(:,1);

J2 = 1.082626925638815e-3;
mu_e = 3.986004415e5; %km^3/s^2
CD = 2;

theta0 = 0; % check if theta0 is initialized at 0 (seems so from p. 2, "second past epoch")
stat101pos_ecef = [-5127.5100 -3794.1600 0.0];
stat337pos_ecef = [3860.9100 3238.4900 3898.0940];
stat394pos_ecef = [549.5050 -1380.8720 6182.1970];
statspos_ecef = [stat101pos_ecef;stat337pos_ecef;stat394pos_ecef];

Xs_ECI = stats_ECI(statspos_ecef,0,theta0);

P0 = diag([1 1 1 1 1 1 1e2 1e6 1 1 1 1 1 1 1 1e-16 1e-16 1e-16]);
R = diag([1e-10 1e-12]);

x0 = [757.700 5222.607 4851.5000 2.21321 4.67834 -5.37130 mu_e J2 CD Xs_ECI(1,1:3) Xs_ECI(2,1:3) Xs_ECI(3,1:3)]';
n = size(x0,1);
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x0;Phi_flat];

dt = t(2) - t(1);
tspan = [0:dt:t(end)];
[~,X_true] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);

j = 1;
x_hat0(1,:) = x0;
x_hat0 = x_hat0';
Deltax_minus = zeros(n,1);
Deltax_plus = Deltax_minus;
Deltax0_plus = ones(n,1);
while j < 2
%Initialize Filter
ti_minus = t(1,1);
x_t_im1 = x_hat0(:,j);
x_hat(:,1) = x_t_im1;
P_plus(:,:,1) = P0;
ri = NaN(numel(t),2);
ri_pre = ri;
ri_post = ri;
% Deltax_i_m1_plus = zeros(7,1);
% numstats = size(rangemeas,2);
% y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
% y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
% Read next observation
for i = 2:numel(t)
    ti = t(i);
%     statnum = find(~isnan(parsedrange(i,2:4)));
    yi = meas(i,3:4);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
%     n = size(x0,1);
    %Propagate from t_i-1 to t_i
%     dt = ti - ti_minus;
    Phi = eye(n);
    Phi_flat = reshape(Phi,size(Phi,1)^2,1);
    Z = [x_t_im1;Phi_flat];
    tspan = [ti_minus:dt/2:ti];
    [~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);
    x_ti = x(end,1:n);
    %Time Update
    Phi = reshape(x(end,n+1:end),n,n);
    Deltax_minus(:,i) = Phi*Deltax_plus(:,i-1);
    P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi';
    %Process Observations
    Hi = zeros(2,n);
    if ~all(isnan(yi))
        if meas(i,2) == 101%~isnan(yi(1)) || ~isnan(yi(2))
            Hi(:,1:9) = Htilde_sc_proj(x_ti(1,1:6),Xs_ECI(1,:));
            Hi(:,10:12) = Htilde_obs_proj(x_ti(1,1:6),Xs_ECI(1,:));
           [hi(i,:)] = predictmeas(x_ti(1,1:6),Xs_ECI(1,:));
        elseif meas(i,2) == 337%~isnan(yi(3)) || ~isnan(yi(4))
            Hi(:,1:9) = Htilde_sc_proj(x_ti(1,1:6),Xs_ECI(2,:));
            Hi(:,13:15) = Htilde_obs_proj(x_ti(1,1:6),Xs_ECI(2,:));
           [hi(i,:)] = predictmeas(x_ti(1,1:6),Xs_ECI(2,:));
        else 
            Hi(:,1:9) = Htilde_sc_proj(x_ti(1,1:6),Xs_ECI(3,:));
            Hi(:,16:18) = Htilde_obs_proj(x_ti(1,1:6),Xs_ECI(3,:));
           [hi(i,:)] = predictmeas(x_ti(1,1:6),Xs_ECI(3,:));
        end
        ri(i,:) = yi - hi(i,:);
        Ki = P_minus(:,:,i)*Hi'*inv(Hi*P_minus(:,:,i)*Hi' + R);
        % measurement update
        Deltax_plus(:,i) = Deltax_minus(:,i) + Ki*(ri(i,:)' - Hi*Deltax_minus(:,i));
        P_plus(:,:,i) = (eye(n) - Ki*Hi)*P_minus(:,:,i)*(eye(n) - Ki*Hi)' + Ki*R*Ki';
    else
        hi(i,:) = NaN(1,2);
        Deltax_plus(:,i) = Deltax_minus(:,i);
        P_plus(:,:,i) = P_minus(:,:,i);
    end
%     [hi(i,:)] = predictmeas(x_ti(1:6),Xs_ECI);
%     ri(i,:) = yi - hi(i,:);
    
%     Hi = zeros(2*numstats,n);
%     for k = 1:numstats
%         if ~isnan(ri(i,2*k))
%             Hi(2*k-1:2*k,1:6) = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(k,:));
%         elseif isnan(ri(i,2*k))
%             ri(i,2*k-1:2*k) = [0 0];
%         end
% %         Hi = Htilde_sc_rho_rhod(x(end,1:6),[Rs_ECI Vs_ECI]);
%     end
%     Ki = P_minus*Hi'*inv(Hi*P_minus*Hi' + R);
%     %Measurement update
%     Deltax_plus(:,i) = Deltax_minus(:,i) + Ki*(ri(i,:)' - Hi*Deltax_minus(:,i));
%     P_plus = (eye(7) - Ki*Hi)*P_minus*(eye(7) - Ki*Hi)' + Ki*R*Ki';
    %set up for next minor iteration
    ti_minus = ti;
    x_t_im1 = x_ti';
    x_hat(:,i) = x_ti' + Deltax_plus(:,i);
    %calculate pre-fit and post-fit residuals
    ri_pre(i,:) = yi - hi(i,:) - (Hi*Deltax_minus(:,i))';
    ri_post(i,:) = yi - hi(i,:) - (Hi*Deltax_plus(:,i))';
%     Pi_mi = P_plus;
%     Deltax_i_m1_plus = Deltaxi_plus;
    
    
end
Phi = eye(n);
Phi_flat = reshape(Phi,n^2,1);
Z = [x_ti';Phi_flat];
tspan = [ti 0];
[~,x] = ode45(@(t,Z) keplerJ2CD_wPhi_ODE(t,Z),tspan,Z,opts);

Phi = reshape(x(end,n+1:end),n,n);
Deltax0_plus = (Phi)*Deltax_plus(:,i);
x_hat0(:,j+1) = x_hat0(:,j) + Deltax0_plus;
j = j + 1;
end
ri_pre_plot = NaN(numel(t),6);
ri_post_plot = NaN(numel(t),6);
for i = 2:numel(t)
    if meas(i,2) == 101
        ri_pre_plot(i,1:2) = ri_pre(i,:);
        ri_post_plot(i,1:2) = ri_post(i,:);
    elseif meas(i,2) == 337
        ri_pre_plot(i,3:4) = ri_pre(i,:);
        ri_post_plot(i,3:4) = ri_post(i,:);
    elseif meas(i,2) == 394
        ri_pre_plot(i,5:6) = ri_pre(i,:);
        ri_post_plot(i,5:6) = ri_post(i,:);
    end
end
x_hat = x_hat';
figure
subplot(2,1,1)
hold on
plot(t,ri_pre_plot(:,1),'*')
plot(t,ri_pre_plot(:,3),'*')
plot(t,ri_pre_plot(:,5),'*')
xlabel('Time (s)')
ylabel('Range Residuals (km)')
grid on
grid minor
legend('Station 101','Station 337','Station 394')
subplot(2,1,2)
hold on
plot(t,ri_pre_plot(:,2),'*')
plot(t,ri_pre_plot(:,4),'*')
plot(t,ri_pre_plot(:,6),'*')
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
grid on 
grid minor
legend('Station 101','Station 337','Station 394')
sgtitle('Pre-fit CKF Measurement Residuals VS Time')

figure
subplot(2,1,1)
hold on
plot(t,ri_post_plot(:,1),'*')
plot(t,ri_post_plot(:,3),'*')
plot(t,ri_post_plot(:,5),'*')
xlabel('Time (s)')
ylabel('Range Residuals (km)')
grid on
grid minor
legend('Station 101','Station 337','Station 394')
subplot(2,1,2)
hold on
plot(t,ri_post_plot(:,2),'*')
plot(t,ri_post_plot(:,4),'*')
plot(t,ri_post_plot(:,6),'*')
xlabel('Time (s)')
ylabel('Range Rate Residuals (km/s)')
grid on 
grid minor
legend('Station 101','Station 337','Station 394')
sgtitle('Post-fit CKF Measurement Residuals VS Time')

%calculate residual RMS values
rangeRMSres = ri_post(~isnan(ri_post(:,1)));
rangeres_RMS = sqrt((1/size(rangeRMSres,1))*sum(rangeRMSres.^2));
rrRMSres = ri_post(:,2);
rrRMSres = rrRMSres(~isnan(rrRMSres));
rrres_RMS = sqrt((1/size(rrRMSres,1))*sum(rrRMSres.^2));
%calculate trace of sequential filter pos. and vel. covariances
for i = 1:size(P_plus,3)
    pos_trace(i) = trace(P_plus(1:3,1:3,i));
    vel_trace(i) = trace(P_plus(4:6,4:6,i));
end
figure
subplot(2,1,1)
semilogy(t,pos_trace)
xlabel('Time (s)')
ylabel('Position Covariance Trace (m^2)')
grid on
grid minor
subplot(2,1,2)
semilogy(t,vel_trace)
xlabel('Time (s)')
ylabel('Velocity Covariance Trace (m/s)^2')
grid on
grid minor
sgtitle('CKF Position and Velocity Covariance Traces VS Time')

%calculate covariance ellipses 
errell = diag(P_plus(:,:,end));
errell = 3*sqrt(errell);
figure
ellipsoid(0,0,0,errell(1),errell(2),errell(3))
axis equal
xlabel('S/C X Position')
ylabel('S/C Y Position')
zlabel('S/C Z Position')
title('Spacecraft Position Covariance Ellipse, 3\sigma,CKF')

figure
ellipsoid(0,0,0,errell(4),errell(5),errell(6))
axis equal
xlabel('S/C X Velocity')
ylabel('S/C Y Velocity')
zlabel('S/C Z Velocity')
title('Spacecraft Velocity Covariance Ellipse, 3\sigma,CKF')

%plot estimated states and final estimated state errors for CKF 
%position state
figure
subplot(3,1,1)
plot(t,x_hat(:,1))
xlabel('Time (s)')
ylabel('r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,x_hat(:,2))
xlabel('Time (s)')
ylabel('r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,x_hat(:,3))
xlabel('Time (s)')
ylabel('r_z (km)')
grid on
grid minor
sgtitle('S/C Position State VS Time,CKF')
%velocity state
figure
subplot(3,1,1)
plot(t,x_hat(:,4))
xlabel('Time (s)')
ylabel('v_x (km/s)')
grid on
grid minor
subplot(3,1,2)
plot(t,x_hat(:,5))
xlabel('Time (s)')
ylabel('v_y (km/s)')
grid on
grid minor
subplot(3,1,3)
plot(t,x_hat(:,6))
xlabel('Time (s)')
ylabel('v_z (km/s)')
grid on
grid minor
sgtitle('S/C Velocity State VS Time,CKF')
% constants
figure
subplot(3,1,1)
plot(t,x_hat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,x_hat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,x_hat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D estimates,CKF')
%station positions
for i = 1:size(x_hat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(x_hat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(x_hat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(x_hat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Estimates VS Time,CKF')
%position state errors
figure
subplot(3,1,1)
plot(t,X_true(:,1)-x_hat(:,1))
xlabel('Time (s)')
ylabel('\Delta r_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,2)-x_hat(:,2))
xlabel('Time (s)')
ylabel('\Delta r_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,3)-x_hat(:,3))
xlabel('Time (s)')
ylabel('\Delta r_z (km)')
grid on
grid minor
sgtitle('S/C Position State Error VS Time,CKF')
%velocity state errors
figure
subplot(3,1,1)
plot(t,X_true(:,4)-x_hat(:,4))
xlabel('Time (s)')
ylabel('v_x (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,5)-x_hat(:,5))
xlabel('Time (s)')
ylabel('v_y (km)')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,6)-x_hat(:,6))
xlabel('Time (s)')
ylabel('v_z (km)')
grid on
grid minor
sgtitle('S/C Velocity State Error VS Time,CKF')
% constants state errors
figure
subplot(3,1,1)
plot(t,X_true(:,7)-x_hat(:,7))
xlabel('Time (s)')
ylabel('\mu_{Earth} (km)')
grid on
grid minor
subplot(3,1,2)
plot(t,X_true(:,8)-x_hat(:,8))
xlabel('Time (s)')
ylabel('J_2')
grid on
grid minor
subplot(3,1,3)
plot(t,X_true(:,9)-x_hat(:,9))
xlabel('Time (s)')
ylabel('C_D')
grid on
grid minor
sgtitle('\mu_{Earth},J_2, and C_D Estimate Errors,CKF')
%station positions errors
for i = 1:size(x_hat,1)
    Xspos_ECEF(i,1:3) = stats_ECEF(X_true(i,10:12)-x_hat(i,10:12),t(i),theta0);
    Xspos_ECEF(i,4:6) = stats_ECEF(X_true(i,13:15)-x_hat(i,13:15),t(i),theta0);
    Xspos_ECEF(i,7:9) = stats_ECEF(X_true(i,16:18)-x_hat(i,16:18),t(i),theta0);
end
figure
subplot(3,1,1)
hold on
plot(t,Xspos_ECEF(:,1))
plot(t,Xspos_ECEF(:,2))
plot(t,Xspos_ECEF(:,3))
xlabel('Time (s)')
ylabel('Station 101 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,2)
hold on
plot(t,Xspos_ECEF(:,4))
plot(t,Xspos_ECEF(:,5))
plot(t,Xspos_ECEF(:,6))
xlabel('Time (s)')
ylabel('Station 337 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
subplot(3,1,3)
hold on
plot(t,Xspos_ECEF(:,7))
plot(t,Xspos_ECEF(:,8))
plot(t,Xspos_ECEF(:,9))
xlabel('Time (s)')
ylabel('Station 394 (km)')
grid on
grid minor
legend('X Position','Y Position','Z Position')
sgtitle('Station State Error VS Time, CKF')
