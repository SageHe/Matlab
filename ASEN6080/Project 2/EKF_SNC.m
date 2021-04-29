%% Project 2 - estimating state for part 2 using EKF 
clear all;clc
tic
% sigma = logspace(-14,-5,15);
% sigma = 4e-8;
sigma = 0;
for ii = 1:numel(sigma) 
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    %load in measurements for part 2
    load('2a_meas')
    load('Project2_Prob2_truth_traj_50days');
    
    data = formatmeas(refmat_data);
    t = data(:,1);
    %Calculate station positions in ECEF -- STATION 2 LONGITUDE IS NEGATIVE
    %ONLY FOR COMPARISON TO EXAMPLE TRUTH DATA!!!
    
    Re = 6378.1363; %earth radius in km
    theta0 = deg2rad(0);
    stat1pos_ecef = [(Re+h1)*cosd(-35.398333)*cosd(148.981944);(Re+h1)*cosd(-35.398333)*sind(148.981944);(Re+h1)*sind(-35.398333)]; 
    stat2pos_ecef = [(Re+h2)*cosd(40.427222)*cosd(-355.749444);(Re+h2)*cosd(40.427222)*sind(-355.749444);(Re+h2)*sind(40.427222)];
    stat3pos_ecef = [(Re+h3)*cosd(35.247164)*cosd(243.205);(Re+h3)*cosd(35.247164)*sind(243.205);(Re+h3)*sind(35.247164)];
    statspos_ecef = [stat1pos_ecef';stat2pos_ecef';stat3pos_ecef'];
    %define uncertainty on range and range rate (5m and 0.5 mm/s,
    %respectively)
    sigma_r = 100;
    sigma_v = 0.1;
    sigma_cr = 0.1;
    
    sigma_rho = 5e-3;
    sigma_rhodot = 0.5e-6;
    %Define a priori covariance and uncertainty matrices
    P0 = diag([sigma_r sigma_r sigma_r sigma_v sigma_v sigma_v sigma_cr]);
    R = diag([2.5e-6 2.5e-13]);
    % Define initial apriori state
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
%     dx = 50*[0.5 -0.5 0.5 0.5e-7 -0.5e-7 0.5e-7]';
%     dx = [5 -5 5 0.5e-6 -0.5e-6 0.5e-6]';
    % j = 1;
    % initialize filter
    x_hat0(1,:) = x0;
    x_hat0 = x_hat0';
    x_hat(1,:) = x_hat0;
    ti_m1 = t(1);
    P_plus(:,:,1) = P0;
    ri = NaN(numel(t),2);
    % read next observation
    for i = 2:numel(t)
        ti = t(i);
        yi = y(i,:);
        Xs_ECI = stats_ECI(statspos_ecef,t(i),theta0);
        % propagate from t_i-1 to t_i
        Phi = eye(n);
        Phi_flat = reshape(Phi,n^2,1);
        Z = [x_hat(i-1,:)';Phi_flat];
        tspan = [0 10];
        [~,x] = ode45(@(t,Z) keplerJ2OOS_wPhi_ODE(t,Z),tspan,Z,opts);
        x_hat(i,:) = x(end,1:n); % technically x_hat_minus I think
        Phi = reshape(x(end,n+1:end),n,n);
        % time update
        %calculate discrete-time P.N
        Qsigma = sigma(ii)*ones(1,3);
        Q_DT = calc_DT_PN(dt,Qsigma);

        P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi' + Q_DT;
        % compute obs sensitiviteis and Kalman gain
        statnum = find(~isnan(yi));
        if ~isempty(statnum)
            statnum = statnum(2)/2;
        end
        Hi = zeros(2,n);
        if ~all(isnan(yi))
            if statnum == 1%~isnan(yi(1)) || ~isnan(yi(2))
                Hi(:,1:6) = Htilde_sc_rho_rhod(x_hat(i,1:6),Xs_ECI(1,:));
               [hi(i,:)] = predictmeas(x_hat(i,1:6),Xs_ECI(1,:));
            elseif statnum == 2%~isnan(yi(3)) || ~isnan(yi(4))
                Hi(:,1:6) = Htilde_sc_rho_rhod(x_hat(i,1:6),Xs_ECI(2,:));
               [hi(i,:)] = predictmeas(x_hat(i,1:6),Xs_ECI(2,:));
            elseif statnum == 3 
                Hi(:,1:6) = Htilde_sc_rho_rhod(x_hat(i,1:6),Xs_ECI(3,:));
               [hi(i,:)] = predictmeas(x_hat(i,1:6),Xs_ECI(3,:));
            end
            inds = find(~isnan(y(i,:)));
            yi = yi(1,inds);

            ri(i,:) = yi - hi(i,:);
            Ki = P_minus(:,:,i)*Hi'*pinv(Hi*P_minus(:,:,i)*Hi' + R);
            % measurement update
            x_hat(i,:) = x_hat(i,:) + (Ki*ri(i,:)')';
            P_plus(:,:,i) = (eye(n) - Ki*Hi)*P_minus(:,:,i)*(eye(n) - Ki*Hi)' + Ki*R*Ki';

            if statnum == 1
                hi(i,:) = predictmeas(x_hat(i,1:6),Xs_ECI(1,:));
            elseif statnum == 2
                hi(i,:) = predictmeas(x_hat(i,1:6),Xs_ECI(2,:));
            elseif statnum == 3
                hi(i,:) = predictmeas(x_hat(i,1:6),Xs_ECI(3,:));
            end
            ri(i,:) = yi - hi(i,:);
        else
            P_plus(:,:,i) = P_minus(:,:,i);
        end
        covbounds(i,:) = [3*sqrt(P_plus(1,1,i)) 3*sqrt(P_plus(2,2,i)) 3*sqrt(P_plus(3,3,i)) 3*sqrt(P_plus(4,4,i)) 3*sqrt(P_plus(5,5,i)) 3*sqrt(P_plus(6,6,i))]; 
    end

    figure
    subplot(3,1,1)
    plot(t,x_hat(:,1))
    subplot(3,1,2)
    plot(t,x_hat(:,2))
    subplot(3,1,3)
    plot(t,x_hat(:,3))

    figure
    subplot(3,1,1)
    hold on
    plot(t,x_hat(:,1) - truth_state(:,1))
    plot(t,covbounds(:,1),'--r')
    plot(t,-covbounds(:,1),'--r')
%     ylim([-.5 .5])
    grid on
    xlabel('Tims (s)')
    ylabel('\delta r_x')
    subplot(3,1,2)
    hold on
    plot(t,x_hat(:,2) - truth_state(:,2))
    plot(t,covbounds(:,2),'--r')
    plot(t,-covbounds(:,2),'--r')
%     ylim([-.5 .5])
    grid on
    xlabel('Tims (s)')
    ylabel('\delta r_y')
    subplot(3,1,3)
    hold on
    plot(t,x_hat(:,3) - truth_state(:,3))
    plot(t,covbounds(:,3),'--r')
    plot(t,-covbounds(:,3),'--r')
%     ylim([-.5 .5])
    grid on
    xlabel('Tims (s)')
    ylabel('\delta r_z')
    sgtitle('Position State Error with 3\sigma Covariance Bounds, EKF')
    
    figure
    subplot(3,1,1)
    hold on
    plot(t,x_hat(:,4) - truth_state(:,4))
    plot(t,covbounds(:,4),'--r')
    plot(t,-covbounds(:,4),'--r')
    ylim([-1e-3 1e-3])
    grid on
    xlabel('Tims (s)')
    ylabel('\delta v_x')
    subplot(3,1,2)
    hold on
    plot(t,x_hat(:,5) - truth_state(:,5))
    plot(t,covbounds(:,5),'--r')
    plot(t,-covbounds(:,5),'--r')
    ylim([-1e-3 1e-3])
    grid on
    xlabel('Tims (s)')
    ylabel('\delta v_x')
    subplot(3,1,3)
    hold on
    plot(t,x_hat(:,6) - truth_state(:,6))
    plot(t,covbounds(:,6),'--r')
    plot(t,-covbounds(:,6),'--r')
    ylim([-1e-3 1e-3])
    grid on
    xlabel('Tims (s)')
    ylabel('\delta v_x')
    sgtitle('Velocity State Error with 3\sigma Covariance Bounds, EKF')

    figure
    subplot(2,1,1)
    plot(t,ri(:,1),'*')
%     ylim([-.002 .002])
    ylim([-0.01 0.01])
    xlabel('Time (s)')
    ylabel('Range')
    grid on
    subplot(2,1,2)
    plot(t,ri(:,2),'*')
%     ylim([-.5e-5 .5e-5])
    ylim([-1e-5 1e-5])
    xlabel('Time (s)')
    ylabel('Range Rate')
    grid on
    sgtitle('Post-fit Measurement Residuals VS Time, EKF')

    %calculate post-fit measurement residual RMS values for measurements
    for i = 1:size(ri,1)
        if isnan(ri(i,:))
            ri(i,:) = [0 0];
        end
        posnorm(i) = norm(x_hat(i,1:3)-truth_state(i,1:3));
        velnorm(i) = norm(x_hat(i,4:6)-truth_state(i,4:6));
    end
    rangeRMS(ii) = sqrt((1/size(ri,1))*sum(ri(:,1).^2));
    rrRMS(ii) = sqrt((1/size(ri,1))*sum(ri(:,2).^2));
    posRMS(ii) = sqrt((1/size(x_hat,1))*sum(posnorm.^2));
    velRMS(ii) = sqrt((1/size(x_hat,1))*sum(velnorm.^2));
    
    clear x_hat0
end
figure
semilogx(sigma,rangeRMS*1e3)
hold on
semilogx(sigma,rrRMS*1e6)
semilogx(sigma,posRMS*1e3)
semilogx(sigma,velRMS*1e6)
xlabel('Tuned Parameter \sigma')
ylabel('RMS value (m or mm/s)')
legend('$\rho$','$\dot{\rho}$','r','v','Interpreter','latex')
title('Measurement and 3D State Residual RMS Values VS \sigma')
grid on
grid minor
toc