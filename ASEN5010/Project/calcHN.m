function [HN] = calcHN(t)
global rLMO omegaLMO iLMO thetaLMO thetadotLMO
% thetadot = sqrt(42828.3/rLMO^3);
theta = thetaLMO + thetadotLMO*t;
[rLM rdotLM] = velandpos(rLMO,omegaLMO,iLMO,theta);
Bt1 = [rLMO 0 0]';
Bt1 = Bt1/norm(Bt1);
vel = [0 rLMO*thetadotLMO 0]';
Bt2 = cross(Bt1,vel)/norm(cross(Bt1,vel));
Bt3 = cross(Bt2,Bt1);

Nt1 = rLM;
Nt1 = Nt1/norm(Nt1);
Nt2 = cross(rLM,rdotLM)/norm(cross(rLM,rdotLM));
Nt3 = cross(Nt2,Nt1);

HT = [Bt1 Bt2 Bt3];
NT = [Nt1 Nt2 Nt3];

HN = HT*NT';
end