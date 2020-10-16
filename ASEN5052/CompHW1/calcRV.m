function [R,V] = calcRV(E,P,e,ehat,ehatperp,mu)
f = atan2(tan(E/2),sqrt((1-e)/(1+e)))*2;
R = (P/(1+e*cos(f)))*(cos(f)*ehat + sin(f)*ehatperp);
V = sqrt(mu/P)*(-sin(f)*ehat + (e+cos(f))*ehatperp);
end
