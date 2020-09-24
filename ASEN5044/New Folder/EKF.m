%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theodore Trozinski
% Maya West
% ASEN 5044 Estimation; Final Project Extended Kalman Filter
% Created: April 17, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xp,sigma,NEES,NIS,Pp] = EKF(x0,yIn,P0,Q,R,Gamma,dt,mu,xtruth,test)
% Instructions: set k=0, reality, k=1 is start
% therefore, xp(0) = xp(1) = initial guess
% Initialize
n = length(yIn);
xp = zeros(length(x0),n);
xp(:,1) = x0;
Pp{1} = P0;
Omega = dt*Gamma;

for k = 1:(n-1)
    % Setup Predict Step
    % x_(k+1) comes from NL ODE solver
    opts = odeset('MaxStep',1e-2,'AbsTol',1e-12,'RelTol',1e-12);
    tvec = [dt*(k-1) dt*(k)]; %(dt*k):0.1:(dt*(k+1)); % time vector for ODE
    [~,x_NL] = ode45(@(t,x) OrbitNL(t,x,mu),tvec,xp(:,k),opts); % Predict
    
    % Need linearized F(k) about xp(:,k)
    % Solve for A
    X = xp(1,k);
    Y = xp(3,k);
    A = [0 1 0 0;
      mu*(2*X^2-Y^2)/((Y^2+X^2).^(5/2)) 0 3*mu*X*Y/((X^2+Y^2).^(5/2)) 0;
      0 0 0 1;
      3*mu*X*Y/((X^2+Y^2).^(5/2)) 0 mu*(2*Y^2-X^2)/((Y^2+X^2).^(5/2)) 0];
    F{k} = eye(4) + dt*A; % F, dyn. linearization about xp(k)
    
    % PREDICT
    x_(:,k+1) = x_NL(end,:)';
    P_{k+1} = F{k}*Pp{k}*F{k}' + Omega*Q*Omega';
    
    % UPDATE PREDICTIONS
    % y_{k+1} is real measurement of predicted,
    % H{k+1} is Linearized H about x_(:,k+1)
    t = (k-1)*dt; % check (k+1)*10, maybe k*10?
    [Hin{k+1},~,~,meas{k+1}] = assembleH(t,x_(:,k+1)); % equivalent to h(y_{k+1})
    meas{k+1};
    yIn{k+1};
    [y_{k+1},y{k+1},H{k+1}] = datcheck(meas{k+1},yIn{k+1},Hin{k+1});
    ey{k+1} = y{k+1} - y_{k+1}; % actual data minus predicted Innovation
    
    if isempty(y_{k+1}) == 1 && isempty(y{k+1}) == 1 && isempty(H{k+1}) == 1
        % Break loop, nothing to update. xp = x-
        xp(:,k+1) = x_(:,k+1);
        Pp{k+1} = P_{k+1};
        Sk = zeros(3);
        % If statement to determine number of stations in view
    elseif numel(ey{k+1}) == 3
        % Approx KF Gains from meas. linearization
        K{k+1} = P_{k+1}*H{k+1}' * inv(H{k+1}*P_{k+1}*H{k+1}' + R);
        % For NIS purposes...
        Sk = H{k+1}*P_{k+1}*H{k+1}'+R;
        % Updated total state estimate
        xp(:,k+1) = x_(:,k+1) + K{k+1}*ey{k+1};
        Pp{k+1} = (eye(4)-K{k+1}*H{k+1}) * P_{k+1};
    % if more than one station is in view:
    elseif numel(ey{k+1}) == 6
        % do it once
        K{k+1} = P_{k+1}*H{k+1}(1:3,:)' * inv(H{k+1}(1:3,:)*P_{k+1}*H{k+1}(1:3,:)' + R);
        xp(:,k+1) = x_(:,k+1) + K{k+1}*ey{k+1}(1:3,1);
        Pp{k+1} = (eye(4)-K{k+1}*H{k+1}(1:3,:)) * P_{k+1};
        % For NIS purposes...
        Sk = H{k+1}(1:3,:)*P_{k+1}*H{k+1}(1:3,:)'+R; % it might be possible that we'll need to do it again for this case
        % do it again, using second measurement now all minuses are pluses
        K{k+1} = Pp{k+1}*H{k+1}(4:6,:)' * inv(H{k+1}(4:6,:)*Pp{k+1}*H{k+1}(4:6,:)' + R);
        xp(:,k+1) = xp(:,k+1) + K{k+1}*ey{k+1}(1:3,2);
        Pp{k+1} = (eye(4)-K{k+1}*H{k+1}(4:6,:)) * Pp{k+1};
    elseif numel(y{k+1}) == 0
        xp(:,k+1) = x_(:,k+1);
        Pp{k+1} = Pp{k};
        fprintf('WARNING')
    end
    sigma(:,k+1) = diag(Pp{k+1}.^(1/2));
    NEES(k) = 0;
    NIS(k) = 0;
    if test == 1
    % Calculate NEES and NIS
        exk = xtruth(:,k+1)-xp(:,k+1);
        eyk = ey{k+1}(1:3,1);
    
        NEES(k) = exk'*inv(Pp{k+1})*exk;
        NIS(k) = eyk'*inv(Sk)*eyk;
    end
end
end