%% Homework 4 - CKF without SNC implemented 
clear all;clc
% sigma = logspace(-14,-5,10);
% sigma = 4.29193426012878e-06;
sigma = 0;
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
    while j < 2%norm(Deltax0_plus) > 1e-7
    %Initialize Filter
    ti_minus = t(1,1);
    x_t_im1 = x_hat0(:,j);
    x_nom = x_hat0(:,j);
    x_hat(:,1) = x_t_im1;
    P_plus(:,:,1) = P0;
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
        Phi(:,:,i) = eye(n);
        Phi_flat = reshape(Phi(:,:,i),size(Phi,1)^2,1);
        Z = [x_t_im1;Phi_flat];
        tspan = [ti_minus:dt/2:ti];
        [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);
        x_nom(:,i) = x(end,1:6);
        x_ti = x(end,1:n);
        %Time Update
        Phi(:,:,i) = reshape(x(end,n+1:end),n,n);
        Deltax_minus(:,i) = Phi(:,:,i)*Deltax_plus(:,i-1);
        %calculate discrete-time P.N
        Qsigma = sigma(ii)*ones(1,3);
        Q_DT = calc_DT_PN(dt,Qsigma);

        P_minus(:,:,i) = Phi(:,:,i)*P_plus(:,:,i-1)*Phi(:,:,i)' + Q_DT;
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
        Ki = P_minus(:,:,i)*Hi'*inv(Hi*P_minus(:,:,i)*Hi' + R);
        %Measurement update
        Deltax_plus(:,i) = Deltax_minus(:,i) + Ki*(ri(i,:)' - Hi*Deltax_minus(:,i));
        P_plus(:,:,i) = (eye(n) - Ki*Hi)*P_minus(:,:,i)*(eye(n) - Ki*Hi)' + Ki*R*Ki';
        %set up for next minor iteration
        covbound(i,:) = [3*sqrt(P_plus(1,1,i)) 3*sqrt(P_plus(2,2,i)) 3*sqrt(P_plus(3,3,i)) 3*sqrt(P_plus(4,4,i)) 3*sqrt(P_plus(5,5,i)) 3*sqrt(P_plus(6,6,i))]; 
        
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
    Phi_prop = eye(n);
    Phi_flat = reshape(Phi_prop,n^2,1);
    Z = [x_ti';Phi_flat];
    tspan = [ti 0];
    [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

    Phi_prop = reshape(x(end,n+1:end),n,n);
    Deltax0_plus = (Phi_prop)*Deltax_plus(:,i)
    x_hat0(:,j+1) = x_hat0(:,j) + Deltax0_plus;
    j = j + 1;
    end
end

%     state_error = x_hat - truth_state(:,1:6)';
    state_error = truth_state(:,1:6)' - x_hat;
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
    ylabel('\delta r_x (km)')
    ylim([-.5 .2])
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,2))
    plot(t,covbound(:,2),'--r')
    plot(t,-covbound(:,2),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_y (km)')
    ylim([-.5 .2])
    grid on
    grid minor
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,3))
    plot(t,covbound(:,3),'--r')
    plot(t,-covbound(:,3),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_z (km)')
    grid on
    grid minor
    ylim([-.5 .2])
    sgtitle('Position State Error with 3\sigma Covariance Bound, CKF')
    
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,4))
    plot(t,covbound(:,4),'--r')
    plot(t,-covbound(:,4),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_x (km/s)')
    ylim([-.2e-3 .2e-3])
    grid on
    grid minor
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,5))
    plot(t,covbound(:,5),'--r')
    plot(t,-covbound(:,5),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_y (km/s)')
    ylim([-.2e-3 .2e-3])
    grid on
    grid minor
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,6))
    plot(t,covbound(:,6),'--r')
    plot(t,-covbound(:,6),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_z (km/s)')
    ylim([-.2e-3 .2e-3])
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
    
    %Implement smoothing algorithm
    xhat_i_l = zeros(size(x_hat,1),size(x_hat,2));
    xhat_i_l(:,end) = Deltax_plus(:,end);
    P_i_l = zeros(size(P_plus,1),size(P_plus,2),size(P_plus,3));
    P_i_l(:,:,end) = P_plus(:,:,end);
    x_hat_smoothed = zeros(size(x_hat,1),size(x_hat,2));
    x_hat_smoothed(:,end) = x_hat(:,end);
    for i = numel(t)-1:-1:1
        %Pull CKF values from solution
        P_i_i = P_plus(:,:,i);
        P_ip1_i = P_minus(:,:,i+1); %possibly needs to be P_minus(:,:,i+1)
        PhiS = Phi(:,:,i+1);
        xhat_i_i = Deltax_plus(:,i);
        %smooth at time ti
        Si = P_i_i*PhiS'*inv(P_ip1_i);
        xhat_i_l(:,i) = xhat_i_i + Si*(xhat_i_l(:,i+1) - PhiS*xhat_i_i);
        P_i_l(:,:,i) = P_i_i + Si*(P_i_l(:,:,i+1) - P_ip1_i)*Si';
        covbound(i,:) = [3*sqrt(P_i_l(1,1,i)) 3*sqrt(P_i_l(2,2,i)) 3*sqrt(P_i_l(3,3,i)) 3*sqrt(P_i_l(4,4,i)) 3*sqrt(P_i_l(5,5,i)) 3*sqrt(P_i_l(6,6,i))]; 
        x_hat_smoothed(:,i) = x_nom(:,i) + xhat_i_l(:,i);
    end
    
    clear state_error
    state_error = (truth_state(:,1:6)' - x_hat_smoothed)';
    
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,1))
    plot(t,covbound(:,1),'--r')
    plot(t,-covbound(:,1),'--r')
    grid on
    grid minor
    xlabel('Time (s)')
    ylabel('\delta r_x (km)')
    ylim([-.3 .3])
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,2))
    plot(t,covbound(:,2),'--r')
    plot(t,-covbound(:,2),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_y (km)')
    ylim([-.4 .4])
    grid on
    grid minor
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,3))
    plot(t,covbound(:,3),'--r')
    plot(t,-covbound(:,3),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_z (km)')
    ylim([-.3 .3])
    grid on
    grid minor
    sgtitle('Position State Error with 3\sigma Covariance Bound, CKF Smoothed')
    
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,4))
    plot(t,covbound(:,4),'--r')
    plot(t,-covbound(:,4),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_x (km/s)')
    ylim([-.2e-3 .2e-3])
    grid on
    grid minor
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,5))
    plot(t,covbound(:,5),'--r')
    plot(t,-covbound(:,5),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_y (km/s)')
    ylim([-.2e-3 .2e-3])
    grid on
    grid minor
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,6))
    plot(t,covbound(:,6),'--r')
    plot(t,-covbound(:,6),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_z (km/s)')
    ylim([-.2e-3 .2e-3])
    grid on
    grid minor
    sgtitle('Velocity State Error with 3\sigma Covariance Bound, CKF Smoothed')

    figure
    error_difference = ((truth_state(:,1)-x_hat(1,:)')-(truth_state(:,1)-x_hat_smoothed(1,:)'));
    plot(t,((truth_state(:,1)-x_hat(1,:)')-(truth_state(:,1)-x_hat_smoothed(1,:)')))
    
    
    clear x_hat0

