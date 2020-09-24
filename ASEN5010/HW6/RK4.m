function [X,u] = RK4(sigma_BN0,Bw_BN0,tstart,tend,dt,gains)
global rLMO omegaLMO iLMO thetaLMO thetadotLMO rGMO omegaGMO iGMO thetaGMO thetadotGMO
sigma_BN = [];
Bw_BN = [];
I = [10 0 0;0 5 0;0 0 7.5];
t = tstart;
f = @(sigma,w,u) [.25*((1-sigma'*sigma)*eye(3) + 2*tilde(sigma) + 2*sigma*sigma')*w... 
    inv(I)*-tilde(w)*I*w + inv(I)*u];
sigma = sigma_BN0;
w = Bw_BN0;
X = [sigma w];
n = 1;
while t < tend
    [RLMO,~] = velandpos(rLMO,omegaLMO,iLMO,thetaLMO+thetadotLMO*t);
    [RGMO,~] = velandpos(rGMO,omegaGMO,iGMO,thetaGMO+thetadotGMO*t);
    if RLMO(2) > 0
        RN = [-1 0 0;0 0 1;0 1 0];
        Nw_RN = [0 0 0]';
    elseif acosd(dot(RLMO,RGMO)/(norm(RLMO)*norm(RGMO))) < 35
        RN = calcRcN(t);
        Nw_RN = calcNwRcN(t);
    else
    RN = calcRnN(t);
    Nw_RN = calcNwRnN(t);
    end
    [sigma_BR,Bw_BR] = att_err(sigma,w,RN,Nw_RN);
    u = -gains(2)*sigma_BR-gains(1)*Bw_BR;
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
end
% sigma_BN = X(:,1,end);
% Bw_BN = X(:,2,end);
u = u;