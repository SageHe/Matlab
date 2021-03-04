function [a,e,i,w,W,P,ehat,ehatperp,f] = compOE(r0,v0,mu)
% mu = 4e5; %km^3/s^2
h = cross(r0,v0); %compute ang. mom.
hhat = h/norm(h);
P = norm(h)^2/mu;
i = acos(dot(hhat,[0 0 1])); %compute inclination
n_Omega = (cross([0 0 1],h))/(norm(cross([0 0 1],h)));
n_Omega_hat = n_Omega/norm(n_Omega);
W = atan2((dot(n_Omega,[0 1 0])),(dot(n_Omega,[1 0 0]))); %compute Omega
e = ((1/mu)*(cross(v0,h))) - (r0/norm(r0)); %compute omega
ehat = e/norm(e);
n_Omega_perp = cross(hhat,n_Omega_hat);
n_Omgegaperp_hat = n_Omega_perp/norm(n_Omega_perp);
w = atan2((dot(ehat,n_Omega_perp)),(dot(ehat,n_Omega)));
a = P/(1 - norm(e)^2);
ehatperp = cross(hhat,ehat);
n_Omgegaperp_hat = n_Omega_perp/norm(n_Omega_perp);
f = atan2((dot(r0,ehatperp)),(dot(r0,ehat)));
% E = 2*atan2(sqrt(1-norm(e))*tan(f/2),sqrt((1+norm(e))));
% n = sqrt(mu/a^3);
% tau = t0 - (1/n)*(E-norm(e)*sin(E));
end