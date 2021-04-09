clear all;clc
tic
%load in measurements from HW2-3
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('J3_Meas_Noised')
%clean up measurement data 
rangemeasnoisy(14319:end,:) = []; %Clean up end of measurement file that has no measurements
rrmeasnoisy(14319:end,:) = [];
t(14319:end) = [];
truth_state(14319:end,:) = [];

mu = 3.986004415e5;
J2 = 1082.63e-6;
J3 = -2.5323e-6;
%Calculate station positions in ECEF 
Re = 6378; %earth radius in km
theta0 = deg2rad(122);
stat1pos_ecef = [Re*cosd(-35.398333)*cosd(148.981944);Re*cosd(-35.398333)*sind(148.981944);Re*sind(-35.398333)]; 
stat2pos_ecef = [Re*cosd(40.427222)*cosd(355.749444);Re*cosd(40.427222)*sind(355.749444);Re*sind(40.427222)];
stat3pos_ecef = [Re*cosd(35.247164)*cosd(243.205);Re*cosd(35.247164)*sind(243.205);Re*sind(35.247164)];
statspos_ecef = [stat1pos_ecef';stat2pos_ecef';stat3pos_ecef'];

sigma_rho = 1e-3;
sigma_rhodot = 1e-6;
sigma_r = 1;
sigma_v = 1e-3;
%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2 1e-12 1e-12]);
R_meas = diag([sigma_rho^2 sigma_rhodot^2]);

%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
C = 0;
x0 = [r;v];
n = size(x0,1);
dt = t(2) - t(1);

% dx = zeros(7,1);
dx = [0.5 -0.5 0.5 0.5e-3 -0.5e-3 0.5e-3]';
% dx =  [-0.0892392099513218; 0.0988127534016745; -0.00483070366315786; -4.61709120432531e-05;-4.17141733572296e-05;0.000146192215370242];
% dx = [5 -5 5 0.5e-6 -0.5e-6 0.5e-6 0 0]';
% dx = [0.8 0.6 -0.2 1e-6 1e-7 -6e-7]';
x_hat0(1,:) = x0 + dx;
dx_minus = zeros(n,1);
dxc_minus = dx_minus;
dx_plus = dx_minus;
Deltax0_plus = ones(n,1);
dc = J3;
j = 1;
x_hat0(j,:) = x0 + dx;
x_hat0 = x_hat0';
while norm(Deltax0_plus) > 1e-5
%Initialize filter
x_t_im1 = x_hat0(:,j);
x_hat(:,1) = x_t_im1;
x_hat_plus(:,1) = x_t_im1;
Pxx_plus(:,:,1) = 1e4*eye(n);
Px_plus(:,:,1) = Pxx_plus(:,:,1);
Pxc_plus(:,1) = zeros(n,1);
Pcc_minus = J3^2;
Sxc_plus(:,1) = zeros(n,1);

ri = NaN(numel(t),6);

y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
for i = 2:numel(t)
    %read next observation
    yi = y(i,:);
    Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    inds = find(~isnan(yi));
%     statnum = inds(2)/2;
    %integrate reference traj.
    Phi(:,:,i) = eye(n);
    Phi_flat = reshape(Phi(:,:,i),n^2,1);
    Theta(:,i) = zeros(n,1);
    Theta_flat = Theta(:,i);
    Z = [x_t_im1;Phi_flat;Theta_flat];
    tspan = [0:dt/2:10];
    [~,x] = ode45(@(t,Z) kepler_wPhiandTheta_ode(t,Z,n),tspan,Z,opts);
    x_ti = x(end,1:n);
    Phi(:,:,i) = reshape(x(end,n+1:(n^2 + n)),n,n);
    Theta(:,i) = x(end,(n^2 + n + 1):end);
    %compute time update
    Px_minus(:,:,i) = Phi(:,:,i)*Px_plus(:,:,i-1)*Phi(:,:,i)';
    Sxc_minus(:,i) = Phi(:,:,i)*Sxc_plus(:,i-1) + Theta(:,i);
    Pxx_minus(:,:,i) = Px_minus(:,:,i) + Sxc_minus(:,i)*Pcc_minus*Sxc_minus(:,i)';
    Pxc_minus(:,i) = Sxc_minus(:,i)*Pcc_minus;
    %report dx_minus and dxc_minus
    dx_minus(:,i) = Phi(:,:,i)*dx_plus(:,i-1);
    dxc_minus(:,i) = dx_minus(:,i) + Sxc_minus(:,i)*dc;
    if ~all(isnan(yi))
        %measurement update
        [hi(i,:)] = predictmeas(x_ti(1:6),Xs_ECI);
        ri(i,:) = yi - hi(i,:);
        statnum = inds(2)/2; 
        Hxi = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(statnum,:));
%         Hxi = [Hxi zeros(2,1)];
        Hci = zeros(2,1);
        K(:,:,i) = Px_minus(:,:,i)*Hxi'*pinv(Hxi*Px_minus(:,:,i)*Hxi' + R_meas);
        Px_plus(:,:,i) = (eye(n) - K(:,:,i)*Hxi)*Px_minus(:,:,i)*(eye(n) - K(:,:,i)*Hxi)' + K(:,:,i)*R_meas*K(:,:,i)';
        Sxc_plus(:,i) = (eye(n) - K(:,:,i)*Hxi)*Sxc_minus(:,i) - K(:,:,i)*Hci;
        Pxx_plus(:,:,i) = Px_plus(:,:,i) + Sxc_plus(:,i)*Pcc_minus*Sxc_plus(:,i)';
        Pxc_plus(:,i) = Sxc_plus(:,i)*Pcc_minus;
        %report dx_plus and dxc_plus
        dx_plus(:,i) = dx_minus(:,i) + K(:,:,i)*(ri(i,inds)' - Hxi*dx_minus(:,i));
        dxc_plus(:,i) = dx_plus(:,i) + Sxc_plus(:,i)*dc;
    else
        dx_plus(:,i) = dx_minus(:,i);
        dxc_plus(:,i) = dxc_minus(:,i);
        Px_plus(:,:,i) = Px_minus(:,:,i);
        Sxc_plus(:,i) = Sxc_minus(:,i);
        Pxx_plus(:,:,i) = Pxx_minus(:,:,i);
        Pxc_plus(:,i) = Pxc_minus(:,i);
    end
    
    x_t_im1 = x_ti';
    x_hat(:,i) = x_ti' + dx_plus(:,i);
    
    Px_cov(i,:) = [2*sqrt(Px_plus(1,1,i)) 2*sqrt(Px_plus(2,2,i)) 2*sqrt(Px_plus(3,3,i)) 2*sqrt(Px_plus(4,4,i)) 2*sqrt(Px_plus(5,5,i)) 2*sqrt(Px_plus(6,6,i))];
    Pxx_cov(i,:) = [2*sqrt(Pxx_plus(1,1,i)) 2*sqrt(Pxx_plus(2,2,i)) 2*sqrt(Pxx_plus(3,3,i)) 2*sqrt(Pxx_plus(4,4,i)) 2*sqrt(Pxx_plus(5,5,i)) 2*sqrt(Pxx_plus(6,6,i))];
    posnorm(i,:) = norm(x_hat(1:3,i)' - truth_state(i,1:3));
    velnorm(i,:) = norm(x_hat(4:6,i)' - truth_state(i,4:6));

end
Phiprop = eye(n);
Phi_flat = reshape(Phiprop,n^2,1);
Theta_flat = zeros(n,1);
Z = [x_ti';Phi_flat;Theta_flat];
tspan = [t(i) 0];
[~,x] = ode45(@(t,Z) kepler_wPhiandTheta_ode(t,Z,n),tspan,Z,opts);

Phiprop = reshape(x(end,n+1:(n^2 + n)),n,n);
Thetaprop = x(end,(n^2+n+1):end)';
Deltax0_plus = (Phiprop)*dxc_plus(:,i)
norm(Deltax0_plus)
x_hat0(:,j+1) = x_hat0(:,j) + Deltax0_plus;
j = j + 1;
end
% plot post-fit state error with 2 sigma cov. envelopes and include rms
% errors
%% 
state_error = x_hat - truth_state(:,1:6)';
% Px_cov = [2*sqrt
%position state errors
figure
subplot(3,1,1)
hold on
plot(t,state_error(1,:))
plot(t,Px_cov(:,1),'--r')
plot(t,Pxx_cov(:,1),'--g')
plot(t,-Px_cov(:,1),'--r')
plot(t,-Pxx_cov(:,1),'--g')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta r_x (km)')
legend('State Error','Px^+','Pxx^+','FontSize',12)
ylim([-1 1])
subplot(3,1,2)
hold on
plot(t,state_error(2,:))
plot(t,Px_cov(:,2),'--r')
plot(t,Pxx_cov(:,2),'--g')
plot(t,-Px_cov(:,2),'--r')
plot(t,-Pxx_cov(:,2),'--g')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta r_y (km)')
legend('State Error','Px^+','Pxx^+','FontSize',12)
ylim([-1 1])
subplot(3,1,3)
hold on
plot(t,state_error(3,:))
plot(t,Px_cov(:,3),'--r')
plot(t,Pxx_cov(:,2),'--g')
plot(t,-Px_cov(:,3),'--r')
plot(t,-Pxx_cov(:,3),'--g')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta r_z (km)')
legend('State Error','Px^+','Pxx^+','FontSize',12)
ylim([-1 1])
sgtitle('CCA Sequential Filter Position State Errors')
%velocity state errors
figure
subplot(3,1,1)
hold on
plot(t,state_error(4,:))
plot(t,Px_cov(:,4),'--r')
plot(t,Pxx_cov(:,4),'--g')
plot(t,-Px_cov(:,4),'--r')
plot(t,-Pxx_cov(:,4),'--g')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta v_x (km/s)')
ylim([-4e-4 4e-4])
legend('State Error','Px^+','Pxx^+','FontSize',12)
subplot(3,1,2)
hold on
plot(t,state_error(5,:))
plot(t,Px_cov(:,5),'--r')
plot(t,Pxx_cov(:,5),'--g')
plot(t,-Px_cov(:,5),'--r')
plot(t,-Pxx_cov(:,5),'--g')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta v_y (km/s)')
ylim([-4e-4 4e-4])
legend('State Error','Px^+','Pxx^+','FontSize',12)
subplot(3,1,3)
hold on
plot(t,state_error(6,:))
plot(t,Px_cov(:,6),'--r')
plot(t,Pxx_cov(:,6),'--g')
plot(t,-Px_cov(:,6),'--r')
plot(t,-Pxx_cov(:,6),'--g')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta v_z (km/s)')
ylim([-4e-4 4e-4])
legend('State Error','Px^+','Pxx^+','FontSize',12)
sgtitle('CCA Sequential Filter Velocity State Errors')
% velnorm(3) = [];
posRMS = sqrt(((1/(size(x_hat,2))))*sum(posnorm.^2));
velRMS = sqrt(((1/(size(x_hat,2))))*sum(velnorm.^2));

% velnorm(3:4) = [];
% posRMS = sqrt(((1/(size(x_hat,2)-200)))*sum(posnorm(201:end).^2));
% velRMS = sqrt(((1/(size(x_hat,2)-2)))*sum(velnorm.^2));
% 

% part d -- Use Psi to map Pc from tf to t0 and back again (see assignment)
Pc_plus_tf = [Pxx_plus(:,:,i) Pxc_plus(:,i);Pxc_plus(:,i)' Pcc_minus];
Psi_tf_t0 = [Phiprop Thetaprop;zeros(1,n) 1];
Pc_plus_t0 = Psi_tf_t0*Pc_plus_tf*Psi_tf_t0';
Pxx_plus_new(:,:,1) = Pc_plus_t0(1:6,1:6,1);
% integrate new state est. forwards in time from epoch
clear x
x_t_im1 = x_hat0(:,j);
Phi(:,:,i) = eye(n);
Phi_flat = reshape(Phi(:,:,i),n^2,1);
Theta(:,i) = zeros(n,1);
Theta_flat = Theta(:,i);
Z = [x_t_im1;Phi_flat;Theta_flat];
% tspan = [0 t(end)];
tspan = linspace(0,t(end),numel(t));
[t_int,x] = ode45(@(t,Z) kepler_wPhiandTheta_ode(t,Z,n),tspan,Z,opts);
Phi_t0_tf = reshape(x(end,n+1:(n^2 + n)),n,n);
covbound(1,:) = [2*sqrt(Pxx_plus_new(1,1,1)) 2*sqrt(Pxx_plus_new(2,2,1)) 2*sqrt(Pxx_plus_new(3,3,1)) 2*sqrt(Pxx_plus_new(4,4,1)) 2*sqrt(Pxx_plus_new(5,5,1)) 2*sqrt(Pxx_plus_new(6,6,1))]; 
clear posnorm velnorm
for i = 2:size(x,1)
    Psi(:,:,i) = [reshape(x(i,n+1:(n^2 + n)),n,n) x(i,(n^2 + n + 1):end)'; zeros(1,n) 1];
    Pc_plus(:,:,i) = Psi(:,:,i)*Pc_plus_t0*Psi(:,:,i)';
    Pxx_plus_new(:,:,i) = Pc_plus(1:6,1:6,i);
    covbound(i,:) = [2*sqrt(Pxx_plus_new(1,1,i)) 2*sqrt(Pxx_plus_new(2,2,i)) 2*sqrt(Pxx_plus_new(3,3,i)) 2*sqrt(Pxx_plus_new(4,4,i)) 2*sqrt(Pxx_plus_new(5,5,i)) 2*sqrt(Pxx_plus_new(6,6,i))]; 
    posnorm(i,:) = norm(x(i,1:3) - truth_state(i,1:3));
    velnorm(i,:) = norm(x(i,4:6) - truth_state(i,4:6));
    
end
clear state_error

posRMS2 = sqrt(((1/(size(x,2))))*sum(posnorm.^2));
velRMS2 = sqrt(((1/(size(x,2))))*sum(velnorm.^2));


state_error = x(:,1:6) - truth_state(:,1:6);
state_error = state_error';
figure
subplot(3,1,1)
hold on
plot(t,state_error(1,:))
plot(t,covbound(:,1),'--r')
plot(t,-covbound(:,1),'--r')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta r_x (km)')
legend('State Error','Pxx^+','FontSize',12)
ylim([-1 1])
subplot(3,1,2)
hold on
plot(t,state_error(2,:))
plot(t,covbound(:,2),'--r')
plot(t,-covbound(:,2),'--r')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta r_y (km)')
legend('State Error','Pxx^+','FontSize',12)
ylim([-1 1])
subplot(3,1,3)
hold on
plot(t,state_error(3,:))
plot(t,covbound(:,3),'--r')
plot(t,-covbound(:,3),'--r')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta r_z (km)')
legend('State Error','Pxx^+','FontSize',12)
ylim([-1 1])
sgtitle('Re-Integrated Position State Errors')
%velocity state errors
figure
subplot(3,1,1)
hold on
plot(t,state_error(4,:))
plot(t,covbound(:,4),'--r')
plot(t,-covbound(:,4),'--r')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta v_x (km/s)')
ylim([-4e-4 4e-4])
legend('State Error','Pxx^+','FontSize',12)
subplot(3,1,2)
hold on
plot(t,state_error(5,:))
plot(t,covbound(:,5),'--r')
plot(t,-covbound(:,5),'--r')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta v_y (km/s)')
ylim([-4e-4 4e-4])
legend('State Error','Pxx^+','FontSize',12)
subplot(3,1,3)
hold on
plot(t,state_error(6,:))
plot(t,covbound(:,6),'--r')
plot(t,-covbound(:,6),'--r')
grid on
grid minor
xlabel('Time (s)')
ylabel('\delta v_z (km/s)')
ylim([-4e-4 4e-4])
legend('State Error','Pxx^+','FontSize',12)
sgtitle('Re-Integrated Velocity State Errors')

% Theta_t0_tf = x(end,(n^2 + n + 1):end)';
% %built new Psi matrix to map Pc_0 from t0 to tf
% Psi_t0_tf = [Phi_t0_tf Theta_t0_tf;zeros(1,n) 1];
% clear Pc_plus_tf
% Pc_plus_tf = Psi_t0_tf*Pc_plus_t0*Psi_t0_tf';
