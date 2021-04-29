%% Homework 2 - Orbit estimation using CKF, EKF, and Batch filtering -- CKF 
clear all;clc
tic
% sigma = logspace(-14,-5,10);
sigma = 0; % not sure if SNC will help/converge as well for this, but should be needed to prevent filter saturation
% sigma = 8e-7;
for ii = 1:numel(sigma)

    %load in measurements for part 2
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    load('2a_meas')
    load('Project2_Prob2_truth_traj_50days');
    data = formatmeas(refmat_data);
    t = data(:,1);
    %Calculate station positions in ECEF 
    Re = 6378.1363; %earth radius in km
    theta0 = deg2rad(0);
    h1 = 0.691750; h2 = 0.834539; h3 = 1.07114904;

    stat1pos_ecef = [(Re+h1)*cosd(-35.398333)*cosd(148.981944);(Re+h1)*cosd(-35.398333)*sind(148.981944);(Re+h1)*sind(-35.398333)]; 
    stat2pos_ecef = [(Re+h2)*cosd(40.427222)*cosd(-355.749444);(Re+h2)*cosd(40.427222)*sind(-355.749444);(Re+h2)*sind(40.427222)];
    stat3pos_ecef = [(Re+h3)*cosd(35.247164)*cosd(243.205);(Re+h3)*cosd(35.247164)*sind(243.205);(Re+h3)*sind(35.247164)];
    statspos_ecef = [stat1pos_ecef';stat2pos_ecef';stat3pos_ecef'];
    %Add Gaussian noise with specified standard deviations
    sigma_rho = 5e-3;
    sigma_rhodot = 0.5e-6;
    
    sigma_r = 100;
    sigma_v = 0.1;
    sigma_cr = 0.1;
    %Define a priori covariance and uncertainty matrices
    P0 = diag([sigma_r sigma_r sigma_r sigma_v sigma_v sigma_v sigma_cr]);
    R = diag([sigma_rho^2 sigma_rhodot^2 sigma_rho^2 sigma_rhodot^2 sigma_rho^2 sigma_rhodot^2]);
    %initialize a priori state estimate, leave as pos., vel., and Cd for
    %now
    x0 = [-274096790 -92859240 -40199490 32.67 -8.94 -3.88 1.2]';
    n = size(x0,1);
%     dt = t(2) - t(1);

    Phi = eye(size(x0,1));
    Phi_flat = reshape(Phi,size(x0,1)^2,1);
    Z = [x0;Phi_flat];
%     tspan = [0:dt:t(end)];
%     [~,X_true] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);

    % dx = zeros(7,1);
    % dx = [.5 .5 .5 1e-7 1e-7 1e-7 0]';
    % dx = [-.8949639 -0.266738025 -0.667648229 0.0012933487 -0.00051752255 0.00015225489 0]';
%     dx = [0.5 -0.5 0.5 0.5e-3 -0.5e-3 0.5e-3]';
    j = 1;
    x_hat0(1,:) = x0;
    x_hat0 = x_hat0';
    Deltax_minus = zeros(n,1);
    Deltax_plus = Deltax_minus;
    Deltax0_plus = 10*ones(n,1);
    while norm(Deltax0_plus(1:3)) > 2e-2
    %Initialize Filter
    ti_minus = t(1,1);
    x_t_im1 = x_hat0(:,j);
    x_hat(:,1) = x_t_im1;
    P_plus = P0;
    
    y = data(:,2:7);
    numstats = 3;
    % Read next observation
    for i = 2:numel(t) % 50 days:ind. 6495    100 days:ind. 12818    150days: ind. 18683
        ti = t(i);
        dt = ti - ti_minus;
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
        [~,x] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);
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
        covbound(i,:) = [3*sqrt(P_plus(1,1)) 3*sqrt(P_plus(2,2)) 3*sqrt(P_plus(3,3)) 3*sqrt(P_plus(4,4)) 3*sqrt(P_plus(5,5)) 3*sqrt(P_plus(6,6)) 3*sqrt(P_plus(7,7))]; 
        
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
        if rem(i,1000) == 0
            i
        end
%         posnorm(i,:) = norm(x_hat(1:3,i)' - truth_state(i,1:3));
%         posnorm(i,:) = norm(x_hat(1:3,i)' - X_true(i,1:3));
%         velnorm(i,:) = norm(x_hat(4:6,i)' - truth_state(i,4:6));
%         velnorm(i,:) = norm(x_hat(4:6,i)' - X_true(i,4:6));
    end
    toc
    Phi = eye(n);
    Phi_flat = reshape(Phi,n^2,1);
    Z = [x_ti';Phi_flat];
    tspan = [ti 0];
    [~,x] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);

    Phi = reshape(x(end,n+1:end),n,n);
    Deltax0_plus = (Phi)*Deltax_plus(:,i)
    x_hat0(:,j+1) = x_hat0(:,j) + Deltax0_plus;
    j
    j = j + 1;
    end
    ind = find(ri_pre(:,1) == min(ri_pre(:,1)));
    ri_pre(ind,:) = NaN;
%%    
%     tspan = t;
%     X0 = x_hat(:,1);
%     Z = [X0;Phi_flat];
%     
%     [t_check,x_check] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);
    mu_E = 398600.432896939;

%     figure
%     subplot(2,1,1)
%     plot(t,ri_pre(:,1),'*')
%     hold on
%     plot(t,ri_pre(:,3),'*')
%     plot(t,ri_pre(:,5),'*')
%     xlabel('Time (s)')
%     ylabel('Pre-fit Meas. Res. (km)')
%     ylim([-.04 .04])
%     grid on
%     grid minor
%     subplot(2,1,2)
%     plot(t,ri_pre(:,2),'*')
%     hold on
%     plot(t,ri_pre(:,4),'*')
%     plot(t,ri_pre(:,6),'*')
%     xlabel('Time (s)')
%     ylabel('Pre-fit Meas. Res (km/s)')
%     grid on
%     grid minor
%     sgtitle('Pre-Fit Measurement Residuals')
%     
%     figure
%     subplot(2,1,1)
%     plot(t,ri_post_plot(:,1),'*')
%     xlabel('Time (s)')
%     ylabel('Post-fit Meas. Res (km)')
%     ylim([-.04 .04])
%     grid on
%     grid minor
%     subplot(2,1,2)
%     plot(t,ri_post_plot(:,2),'*')
%     xlabel('Time (s)')
%     ylabel('Post-fit Meas. Res. (km/s)')
%     grid on
%     grid minor
%     sgtitle('Post-Fit Measurement Residuals')
%     
%     figure
%     subplot(3,1,1)
%     plot(t,covbound(:,1))
%     hold on
%     plot(t,-covbound(:,1))
%     xlabel('Time (s)')
%     ylabel('km')
%     ylim([-40 40])
%     grid on
%     grid minor
%     subplot(3,1,2)
%     plot(t,covbound(:,2))
%     hold on
%     plot(t,-covbound(:,2))
%     xlabel('Time (s)')
%     ylabel('km')
%     ylim([-40 40])
%     grid on
%     grid minor
%     subplot(3,1,3)
%     plot(t,covbound(:,3))
%     hold on
%     plot(t,-covbound(:,3))
%     xlabel('Time (s)')
%     ylabel('km')
%     ylim([-40 40])
%     grid on
%     grid minor
%     sgtitle('3\sigma Covariance Envelopes, Position States')
% 
%     figure
%     subplot(3,1,1)
%     plot(t,covbound(:,4))
%     hold on
%     plot(t,-covbound(:,4))
%     xlabel('Time (s)')
%     ylabel('km/s')
%     ylim([-.03 .03])
%     grid on
%     grid minor
%     subplot(3,1,2)
%     plot(t,covbound(:,4))
%     hold on
%     plot(t,-covbound(:,4))
%     xlabel('Time (s)')
%     ylabel('km/s')
%     ylim([-.03 .03])
%     grid on
%     grid minor
%     subplot(3,1,3)
%     plot(t,covbound(:,4))
%     hold on
%     plot(t,-covbound(:,4))
%     xlabel('Time (s)')
%     ylabel('km/s')
%     ylim([-.03 .03])
%     grid on
%     grid minor
%     sgtitle('3\sigma Covariance Envelopes, Velocity States')
%     
%     figure
%     plot(t,covbound(:,7))
%     hold on
%     plot(t,-covbound(:,7))
%     grid on
%     grid minor
%     xlabel('Time (s)')
%     ylabel('C_r Covariance')
%     title('3\sigma bound of C_r')

    
    opts = odeset('RelTol',10e-12,'AbsTol',10e-12,'Events',@eventfunc1);
    tspan = [0 (300*24*3600)];
    X0 = x_hat(:,1);
    Z = [X0;Phi_flat];
    [t_3soi,x_3soi,te,ye,ie] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);

    rinf = ye(1:3);
    vinf = ye(4:6);

    [BN,BdotR1,BdotR2,BdotT1,BdotT2,LTOF] = Bplane_calcs(rinf,vinf);

    opts = odeset('RelTol',10e-14,'AbsTol',10e-14,'Events',@(t,Z) eventfunc2(t,Z,BN));
    tspan = [t(end) (300*24*3600)];
    X0 = x_hat(:,end);
    Z = [X0;Phi_flat];
    [t_BPI,x_BPI,te_BPI,ye_BPI,ie_BPI] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);

    Phi_t0_tf = reshape(x_BPI(end,8:end),n,n);

    P_BPI = Phi_t0_tf*P_plus*Phi_t0_tf';
    rcov = P_BPI(1:3,1:3);
    B_P_BPI = BN*rcov*BN';
    B_P_BPI = B_P_BPI(2:3,2:3);

    ell_t = linspace(0,2*pi,200);
    for i = 1:numel(ell_t)
    [V,D] = eig(2*sqrtm(B_P_BPI));
    ell(:,i) = V*[sqrt(D(1,1))*cos(t(i));sqrt(D(2,2))*sin(t(i))];
    end
    
    figure
    plot(ell(1,:),ell(2,:),'*')
   
%    opts = odeset('RelTol',10e-12,'AbsTol',10e-12);
% 
%     tspan = [0 (te + LTOF)];
%     X0 = x_hat(:,1);
%     Z = [X0;Phi_flat];
%     [t_BPI,x_BPI] = ode45(@(t,Z) TBG_SRP_ode(t,Z,n),tspan,Z,opts);

    
%     state_error = x_hat - truth_state(:,1:6)';
    state_error = x_hat - X_true(:,1:6)';
    state_error = state_error';
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,1))
%     plot(t,covbound(:,1),'--r')
%     plot(t,-covbound(:,1),'--r')
    grid on
    grid minor
    xlabel('Time (s)')
    ylabel('\delta r_x')
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,2))
%     plot(t,covbound(:,2),'--r')
%     plot(t,-covbound(:,2),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_y')
    grid on
    grid minor
    xlabel('Time (s)')
    ylabel('y state error')
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,3))
%     plot(t,covbound(:,3),'--r')
%     plot(t,-covbound(:,3),'--r')
    xlabel('Time (s)')
    ylabel('\delta r_z')
    grid on
    grid minor
    sgtitle('Position State Error with 3\sigma Covariance Bound, CKF')
    
    figure
    subplot(3,1,1)
    hold on
    plot(t,state_error(:,4))
%     plot(t,covbound(:,4),'--r')
%     plot(t,-covbound(:,4),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_x')
    grid on
    grid minor
    subplot(3,1,2)
    hold on
    plot(t,state_error(:,5))
%     plot(t,covbound(:,5),'--r')
%     plot(t,-covbound(:,5),'--r')
    xlabel('Time (s)')
    ylabel('\delta v_y')
    grid on
    grid minor
    subplot(3,1,3)
    hold on
    plot(t,state_error(:,6))
%     plot(t,covbound(:,6),'--r')
%     plot(t,-covbound(:,6),'--r')
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