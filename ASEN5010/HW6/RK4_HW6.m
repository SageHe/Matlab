function [X,u,err] = RK4_HW6(sigma_BN0,Bw_BN0,tstart,tend,dt,gains)
% global rLMO omegaLMO iLMO thetaLMO thetadotLMO rGMO omegaGMO iGMO thetaGMO thetadotGMO
sigma_BN = [];
Bw_BN = [];
I = [100 0 0;0 75 0;0 0 80];
t = tstart;
L = [0.5 -0.3 0.2]';
f = @(sigma,w,u) [.25*((1-sigma'*sigma)*eye(3) + 2*tilde(sigma) + 2*sigma*sigma')*w... 
    inv(I)*-tilde(w)*I*w + inv(I)*u];
sigma = sigma_BN0;
w = Bw_BN0;
X = [sigma w];
n = 1;
sigma_BR = sigma_BN0;
Bw_BR = Bw_BN0;
while t < tend
    sigma_RN = [0.2*sin(0.05*t) 0.3*cos(0.05*t) -0.3*sin(0.05*t)]';
    sigma_RN_p1 = [0.2*sin(0.05*(t+dt)) 0.3*cos(0.05*(t+dt)) -0.3*sin(0.05*(t+dt))]';
    sigmadot_RN = [0.2*cos(0.05*t)*0.05 -0.3*sin(0.05*t)*0.05 -0.3*cos(0.05*t)*0.05]';
    sigmadot_RN_p1 = [0.2*cos(0.05*(t+dt))*0.05 -0.3*sin(0.05*(t+dt))*0.05 -0.3*cos(0.05*(t+dt))*0.05]';
    w_RN = calc_wRN(sigma_RN,sigmadot_RN);
    w_RN_p1 = calc_wRN(sigma_RN_p1,sigmadot_RN_p1);
    omegadot = (w_RN_p1 - w_RN)/dt;
    RN = MRP2C(sigma_RN);
    BN = MRP2C(sigma);
    omegadot = BN*omegadot;
    w_RN = BN*w_RN;
%     RN_tp1 = MRP2C(sigma_RN_p1);
%     RN_tp2 = MRP2C(sigma_RN_p2);
%     RN_dot1 = (RN_tp1 - RN)/dt;
%     RN_dot2 = (RN_tp2 -RN_tp1)/dt;
%     omegatilde = -RN_dot1*RN';
%     Nw_RN1 = RN'*[mean([omegatilde(3,2),-omegatilde(2,3)]);
%     mean([-omegatilde(3,1),omegatilde(1,3)]);
%     mean([omegatilde(2,1),-omegatilde(1,2)])];
%     omegatilde = -RN_dot2*RN_tp1';
%     Nw_RN2 = RN_tp1'*[mean([omegatilde(3,2),-omegatilde(2,3)]);
%     mean([-omegatilde(3,1),omegatilde(1,3)]);
%     mean([omegatilde(2,1),-omegatilde(1,2)])];
%     omegadot = (Nw_RN2-Nw_RN1)/dt;
%     [RLMO,~] = velandpos(rLMO,omegaLMO,iLMO,thetaLMO+thetadotLMO*t);
%     [RGMO,~] = velandpos(rGMO,omegaGMO,iGMO,thetaGMO+thetadotGMO*t);
%     if RLMO(2) > 0
%         RN = [-1 0 0;0 0 1;0 1 0];
%         Nw_RN = [0 0 0]';
%     elseif acosd(dot(RLMO,RGMO)/(norm(RLMO)*norm(RGMO))) < 35
%         RN = calcRcN(t);
%         Nw_RN = calcNwRcN(t);
%     else
%     RN = calcRnN(t);
%     Nw_RN = calcNwRnN(t);
%     end
    [sigma_BR,Bw_BR] = att_err(sigma,w,RN,w_RN);
      u = -5*sigma -diag([22.3607 19.3649 20])*w + tilde(w)*I*w;
%     u = -gains(2)*sigma_BR -gains(1)*Bw_BR + I*(omegadot - cross(w,w_RN)) + tilde(w)*I*w;
%     u = -gains(2)*sigma_BR -gains(1)*Bw_BR;
    k1 = dt*f(sigma,w,u);
    k2 = dt*f((sigma+k1(:,1)/2),(w+k1(:,2)/2),u);
    k3 = dt*f((sigma+k2(:,1)/2),(w+k2(:,2)/2),u);
    k4 = dt*f((sigma+k3(:,1)),(w+k3(:,2)),u);
    X(:,:,n+1) = X(:,:,n) + (1/6)*(k1+2*k2+2*k3+k4);
    sigma = X(:,1,n+1);
    if norm(X(:,1,n+1)) > 1
        sigma = X(:,1,n+1);
        sigma = -sigma/(sigma'*sigma);
        X(:,1,n+1) = sigma;
    end
    w = X(:,2,n+1);
    t = t + dt;
    n = n + 1;
    err(:,:,n) = [sigma_BR Bw_BR];
end
% sigma_BN = X(:,1,end);
% Bw_BN = X(:,2,end);
u = u;