%% Homework 2 - Orbit estimation using CKF, EKF, and Batch filtering -- CKF 
clear all;close all;clc
% sigma = logspace(-14,-5,10);
sigma = 4.29193426012878e-06;
% sigma = 8e-7;
for ii = 1:numel(sigma)

    %load in measurements from HW1

    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    % load('HW1_Nom_Measurements_V2')
    load('J3_Measurements')
    % load('HW3_test_measurements')
    rangemeas(14319:end,:) = []; %Clean up end of measurement file that has no measurements
    rrmeas(14319:end,:) = [];
    t(14319:end) = [];
    truth_state(14319:end,:) = [];
    % load('HW1_Noisy_Measurements')
    J2 = 1082.63e-6;
    mu = 3.986004415e5;
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
    R = 4*diag([sigma_rho sigma_rhodot sigma_rho sigma_rhodot sigma_rho sigma_rhodot]);
    %calculate initial true unperturbed state
    [r,v] = calcposvel(10000,0.001,40,80,40,0);
    Period = 2*pi*sqrt((10000^3)/mu);
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
    dx = [0.5 -0.5 0.5 0.5e-3 -0.5e-3 0.5e-3]';
    j = 1;
    x_hat0(1,:) = x0 + dx;
    x_hat0 = x_hat0';
    Deltax_minus = zeros(n,1);
    Deltax_plus = Deltax_minus;
    Deltax0_plus = ones(n,1);
    while norm(Deltax0_plus) > 1e-7
    %Initialize Filter
    ti_minus = t(1,1);
    x_t_im1 = x_hat0(:,j);
    x_hat(:,1) = x_t_im1;
    P_plus = P0;
    % Deltax_i_m1_plus = zeros(7,1);
    numstats = size(rangemeas,2);
    % y = [rangemeas(:,1) rrmeas(:,1) rangemeas(:,2) rrmeas(:,2) rangemeas(:,3) rrmeas(:,3)]; %save noisy measurements below for later
    y = [rangemeasnoisy(:,1) rrmeasnoisy(:,1) rangemeasnoisy(:,2) rrmeasnoisy(:,2) rangemeasnoisy(:,3) rrmeasnoisy(:,3)];
    % Read next observation
    for i = 2:numel(t)
        ti = t(i);
    %     statnum = find(~isnan(parsedrange(i,2:4)));
        yi = y(i,:);
        Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
    %     n = size(x0,1);
        %Propagate from t_i-1 to t_i
    %     dt = ti - ti_minus;
        Phi = eye(n);
        Phi_flat = reshape(Phi,size(Phi,1)^2,1);
        Z = [x_t_im1;Phi_flat];
        tspan = [ti_minus:dt/2:ti];
        [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);
        x_ti = x(end,1:n);
        %Time Update
        Phi = reshape(x(end,n+1:end),n,n);
        Deltax_minus(:,i) = Phi*Deltax_plus(:,i-1);
        %calculate discrete-time P.N
        Qsigma = sigma(ii)*ones(1,3);
        Q_DT = calc_DT_PN(dt,Qsigma);

        P_minus = Phi*P_plus*Phi' + Q_DT;
        %Process Observations
        [hi(i,:)] = predictmeas(x_ti(1:6),Xs_ECI);
        ri(i,:) = yi - hi(i,:);

        Hi = zeros(2*numstats,n);
        for k = 1:numstats
            if ~isnan(ri(i,2*k))
                Hi(2*k-1:2*k,1:6) = Htilde_sc_rho_rhod(x_ti(1:6),Xs_ECI(k,:));
            elseif isnan(ri(i,2*k))
                ri(i,2*k-1:2*k) = [0 0];
            end
    %         Hi = Htilde_sc_rho_rhod(x(end,1:6),[Rs_ECI Vs_ECI]);
        end
        Ki = P_minus*Hi'*inv(Hi*P_minus*Hi' + R);
        %Measurement update
        Deltax_plus(:,i) = Deltax_minus(:,i) + Ki*(ri(i,:)' - Hi*Deltax_minus(:,i));
        P_plus = (eye(n) - Ki*Hi)*P_minus*(eye(n) - Ki*Hi)' + Ki*R*Ki';
        %set up for next minor iteration
        covbound(i,:) = [3*sqrt(P_plus(1,1)) 3*sqrt(P_plus(2,2)) 3*sqrt(P_plus(3,3)) 3*sqrt(P_plus(4,4)) 3*sqrt(P_plus(5,5)) 3*sqrt(P_plus(6,6))]; 
        
        ti_minus = ti;
        x_t_im1 = x_ti';
        x_hat(:,i) = x_ti' + Deltax_plus(:,i);
    %     Pi_mi = P_plus;
    %     Deltax_i_m1_plus = Deltaxi_plus;
    %     inds(i,:) = find(~isnan(y(i,:)));
        ri_pre(i,:) = yi - hi(i,:) - (Hi*Deltax_minus(:,i))';
        ri_post(i,:) = yi - hi(i,:) - (Hi*Deltax_plus(:,i))';
        if ~all(isnan(yi))
            ri_post_plot(i,:) = ri_post(i,~isnan(ri_post(i,:)));
        else
            ri_post_plot(i,:) = NaN(1,2);
        end
        posnorm(i,:) = norm(x_hat(1:3,i)' - truth_state(i,1:3));
        velnorm(i,:) = norm(x_hat(4:6,i)' - truth_state(i,4:6));
    end
    Phi = eye(n);
    Phi_flat = reshape(Phi,n^2,1);
    Z = [x_ti';Phi_flat];
    tspan = [ti 0];
    [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

    Phi = reshape(x(end,n+1:end),n,n);
    Deltax0_plus = (Phi)*Deltax_plus(:,i)
    x_hat0(:,j+1) = x_hat0(:,j) + Deltax0_plus;
    j = j + 1;
    end

    state_error = x_hat - truth_state(:,1:6)';
    % state_error(7,:) = [];
    state_error = state_error';
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,1))
    plot(t,covbound(:,1),'--r')
    plot(t,-covbound(:,1),'--r')
    grid on
    grid minor
    xlabel('Time (s)')
    ylabel('\delta r_x')
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,2))
    plot(t,covbound(:,2),'--r')
    plot(t,-covbound(:,2),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_y')
    grid on
    grid minor
    % xlabel('Time (s)')
    % ylabel('y state error')
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,3))
    plot(t,covbound(:,3),'--r')
    plot(t,-covbound(:,3),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_z')
    grid on
    grid minor
    sgtitle('Position State Error with 3\sigma Covariance Bound, CKF')
    
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,4))
    plot(t,covbound(:,4),'--r')
    plot(t,-covbound(:,4),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_x')
    grid on
    grid minor
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,5))
    plot(t,covbound(:,5),'--r')
    plot(t,-covbound(:,5),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_y')
    grid on
    grid minor
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,6))
    plot(t,covbound(:,6),'--r')
    plot(t,-covbound(:,6),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_z')
    grid on
    grid minor
    sgtitle('Velocity State Error with 3\sigma Covariance Bound, CKF')

    figure
    subplot(2,1,1)
    plot(t,ri_post_plot(:,1),'*')
    grid on
    xlabel('Time (s)')
    ylabel('Range')
    subplot(2,1,2)
    plot(t,ri_post_plot(:,2),'*')
    xlabel('Time (s)')
    ylabel('Range Rate')
%     xlim([0 10e4])
    ylim([-1e-5 1e-5])
    grid on
    sgtitle('CKF Post-fit Residuals, \sigma = 4.2919e-6')
    % for i = 1:size(ri_post,1)
    %     ri_post
    rangeRMS = ri_post_plot(~isnan(ri_post_plot(:,1)));
    rangeRMS = sqrt((1/size(rangeRMS,1))*sum(rangeRMS.^2));
    rangeRMSout(ii) = rangeRMS;
    rrRMS = ri_post_plot(:,2);
    rrRMS = rrRMS(~isnan(rrRMS));
    rrRMS = sqrt((1/size(rrRMS,1))*sum(rrRMS.^2));
    rrRMSout(ii) = rrRMS;

    posRMS(ii) = sqrt((1/size(x_hat,1))*sum(posnorm.^2));
    velRMS(ii) = sqrt((1/size(x_hat,1))*sum(velnorm.^2));
    
    clear x_hat0
end

% figure
% semilogx(sigma,rangeRMSout*1e3)
% hold on
% semilogx(sigma,rrRMSout*1e6)
% xlabel('Tuned Parameter \sigma')
% ylabel('RMS value (m or mm/s)')
% legend('$\rho$','$\dot{\rho}$','Interpreter','latex')
% title('Measurement Residual RMS VS \sigma')
% grid on
% grid minor
% figure
% semilogx(sigma,posRMS*1e3)
% hold on
% semilogx(sigma,velRMS*1e6)
% xlabel('Tuned Parameter \sigma')
% ylabel('RMS value (m or mm/s)')
% legend('r','v')
% title('3D State Error RMS VS \sigma')
% grid on
% grid minor
