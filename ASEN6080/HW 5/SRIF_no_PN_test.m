clear all;clc
%load in measurements from HW1
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('HW2_Nom_Measurements')
load('y_noisy');
y_data = y_noisy(1:6,:);
IDs = y_noisy(7:9,:);
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

rangemeasnoisy = rangemeas + rho_noise;
rrmeasnoisy = rrmeas + rhodot_noise;

%Define a priori covariance and uncertainty matrices
P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 sigma_v^2 sigma_v^2 sigma_v^2]);
R_meas = diag([sigma_rho^2 sigma_rhodot^2]);
R0 = chol(pinv(P0));

Vi = chol(R_meas);
Vi_inv = inv(Vi);
%calculate initial true unperturbed state
[r,v] = calcposvel(10000,0.001,40,80,40,0);
x0 = [r;v];
n = size(x0,1);
dt = t(2) - t(1);

Phi = eye(size(x0,1));
Phi_flat = reshape(Phi,size(x0,1)^2,1);
Z = [x0;Phi_flat];
tspan = [0:dt:t(end)];
[~,X_true] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

% dx = zeros(7,1);
% dx = [.5 .5 .5 1e-7 1e-7 1e-7 0]';
% dx = [-.8949639 -0.266738025 -0.667648229 0.0012933487 -0.00051752255 0.00015225489 0]';
% dx = [0.5 -0.5 0.5 0.5e-3 -0.5e-3 0.5e-3]';
dx = [0.8 0.6 -0.2 1e-6 1e-7 -6e-7]';
j = 1;
x_hat0(1,:) = x0 + dx;
x_hat0 = x_hat0';
Deltab0 = zeros(n,1);
Deltab_plus = Deltab0;
Ri_plus(:,:,1) = R0;
while j < 2 %norm(Deltax0_plus) > 1e-8
%Initialize Filter
ti_minus = t(1);
x_t_im1 = x_hat0(:,j);
x_hat(:,1) = x_t_im1;

numstats = size(rangemeas,2);
% y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
% y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
    for i = 2:numel(t)
        %Read next observation
        ti = t(i);
%         yi = y(i,:);
        IDi = IDs(:,i);
        IDi = IDi(~isnan(IDi));
        yi = y_data(:,i);
        yi = yi(~isnan(yi));
        Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
        %Propagate NL trajectory
        Phi(:,:,i) = eye(n);
        Phi_flat = reshape(Phi(:,:,i),size(Phi,1)^2,1);
        Z = [x_t_im1;Phi_flat];
        tspan = [0:dt/2:10];
        [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);
        x_ti = x(end,1:n);
        Phi(:,:,i) = reshape(x(end,n+1:end),n,n);
        %Time update
        Deltab_minus(:,i) = Deltab_plus(:,i-1);
        Ri_minus(:,:,i) = Ri_plus(:,:,i-1)*pinv(Phi(:,:,i));        %should be same as Ri_plus(:,:,i-1)*inv(Phi(:,:,i))
        %Calculate H and r for pre-whitening
        [hi(i,:)] = predictmeas(x_ti(1:6),Xs_ECI);
%         ri(i,:) = yi - hi(i,:);
        ri(i,:) = yi - hi(i,2*IDi-1:2*IDi)'; 

        Hi = zeros(2*numstats,n);
%         for k = 1:numstats
%             if ~isnan(ri(i,2*k))
%                 Hi(2*k-1:2*k,1:6) = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(k,:));
%             elseif isnan(ri(i,2*k))
%                 ri(i,2*k-1:2*k) = [0 0];
%             end
%         end

        Hi = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(IDi,:));

        
%         for z = size(Hi,2):-1:1
%             if Hi(z,:) == 0
%                 Hi(z,:) = [];
%             end
%         end
        if ~isempty(IDi)
            ri_tilde(i,:) = Vi_inv*nonzeros(ri(i,:));
            Htilde = Vi_inv*Hi;
            
            Ri_loop = Ri_minus(:,:,i);
            Deltab_loop = Deltab_minus(:,i);
        for k = 1:sum(~isnan(yi))
            Aminus = [Ri_loop Deltab_loop;Htilde(k,:) ri_tilde(i,k)];
            [Q,R_decomp] = qr(Aminus);
            Ri_loop = R_decomp(1:6,1:6);
            Deltab_loop = R_decomp(1:6,end);
            e(i,k) = R_decomp(end,end);
        end
        Ri_plus(:,:,i) = Ri_loop;
        Deltab_plus(:,i) = Deltab_loop;
        else
            Ri_plus(:,:,i) = Ri_minus(:,:,i);    
            Deltab_plus(:,i) = Deltab_minus(:,i);
%         Htilde = zeros(2,n);
%         ri_tilde(i,:) = zeros(1,2);
        end
        %Meas. Update for each meas. k
        Deltax_plus(:,i) = pinv(Ri_plus(:,:,i))*Deltab_plus(:,i);                 %getting matrices singular to working precision here, might need to use chol. decomp.?
        P_plus(:,:,i) = pinv(Ri_plus(:,:,i))*(pinv(Ri_plus(:,:,i)))';
        x_t_im1 = x_ti';
        x_hat(:,i) = x_ti' + Deltax_plus(:,i);
    end
end