function [RcN] = calcRcN(t)
global rLMO omegaLMO iLMO thetaLMO thetadotLMO rGMO omegaGMO iGMO thetaGMO thetadotGMO
[rLM,rdotLM] = velandpos(rLMO,omegaLMO,iLMO,thetaLMO+thetadotLMO*t);
[rGM,rdotGM] = velandpos(rGMO,omegaGMO,iGMO,thetaGMO+thetadotGMO*t);

deltar = rGM - rLM;
Nt1 = -deltar/norm(deltar);
nhat_3 = [0 0 1]';
Nt2 = cross(deltar,nhat_3)/norm(cross(deltar,nhat_3));
Nt3 = cross(Nt1,Nt2);

RcN = [Nt1 Nt2 Nt3]';
end