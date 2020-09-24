%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theodore Trozinski
% Maya West
% ASEN 5044 Estimation; Final Project Linearized Kalman Filter
% Created: April 6, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx, sigma, NEES, NIS, Pout] = KF(yIn,ystarIn,dx0,Fcells,P0,HIn,n,Q,R,Gamma,xNL,xstar,ytruth,test)
% Kalman filter for use with no forcing function, Linear
N = length(dx0);
dx_ = zeros(N,n);
dx_(:,1) = dx0;
dx = zeros(N,n);
dx(:,1) = dx0;
P = P0;
Pout{1} = P0;
dt = 10;
sigma = zeros(N,n);
Omega = dt*Gamma;
% Loop over each iteration:
for k = 1:(n-1)
% k
    % Update Mats
    F = Fcells{k};
    
    [ystar{k+1},y{k+1},H{k+1}] = datcheck(ystarIn{k+1},yIn{k+1},HIn{k+1});
    
    dy = y{k+1}-ystar{k+1};
    
    % Correct k+1
    if numel(dy) == 3
        % Predict k+1
        dx_(:,k+1) = F*dx(:,k);
        P_ = F*P*(F') + Omega*Q*Omega';
        K = P_*H{k+1}' * inv(H{k+1}*P_*H{k+1}'+R);
        
        % For NIS purposes...
        Sk = H{k+1}*P_*H{k+1}'+R;
        
        % Update k+1
        dx(:,k+1) = dx_(:,k+1) + K * (dy-H{k+1}*dx_(:,k+1));
        P = (eye(N)-K*H{k+1})*P_;

    % for cases with 2 ground stations tracking:
    elseif numel(dy) == 6
        dy = cat(1,dy(1:3,1),dy(1:3,2));
        % Predict
        dx_(:,k+1) = F*dx(:,k);
        P_ = F*P*(F') + Gamma*Q*Gamma';
        K = P_*H{k+1}(1:3,:)' * inv(H{k+1}(1:3,:)*P_*H{k+1}(1:3,:)'+R);
        % For NIS purposes...
        Sk = H{k+1}(1:3,:)*P_*H{k+1}(1:3,:)'+R; % it might be possible that we'll need to do it again for this case
        % Update the First
        dx(:,k+1) = dx_(:,k+1) + K * (dy(1:3)-H{k+1}(1:3,:)*dx_(:,k+1));
        P = (eye(N)-K*H{k+1}(1:3,:))*P_;
        % Predict a Second Time
        dx_(:,k+1) = dx(:,k+1);
        P_ = P;
        K = P_*H{k+1}(4:6,:)' * inv(H{k+1}(4:6,:)*P_*H{k+1}(4:6,:)'+R);
        % For NIS purposes...
        %Sk = blkdiag(Sk,H{k+1}(4:6,:)*P_*H{k+1}(4:6,:)'+R);
        % Update Again
        dx(:,k+1) = dx_(:,k+1) + K * (dy(4:6)-H{k+1}(4:6,:)*dx_(:,k+1));
        P = (eye(N)-K*H{k+1}(4:6,:))*P_;
    end
    Pout{k+1} = P;
    % Store Covariances:
    sigma(:,k+1) = diag(P).^(1/2);
    
    
    if test == 1
        % Calculate NEES and NIS
        exk = xNL(:,k+1) - xstar(:,k+1) - dx(:,k+1);
        eyk = dy(1:3);%yIn{k+1}(1:3,1) - (ystarIn{k+1}(1:3,1) + dy(1:3,1));
    
        NEES(k) = exk'*inv(P)*exk;
        NIS(k) = eyk'*inv(Sk)*eyk;
    else 
        NEES = 0;
        NIS = 0;
    end
    

end



end