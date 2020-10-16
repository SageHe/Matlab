clear all;close all;clc
%% Problem 1 
% Write program that takes initial S/C position and velocity at initial
% time and predict future/past position of S/C at time t
% Part a -- Calculate orbital elements from arbitrary position vector,
% velocity vecor, and time
mu = 4e5; %Kg^3/s^2
r0 = [6e3 6e3 6e3];
v0 = [-5 5 0];
t0 = 0;
h = cross(r0,v0);
[a,e,i,w,W,tau,P,ehat,ehatperp] = compOE(r0,v0,t0,mu);
n = sqrt(mu/a^3);
T = 2*pi*sqrt(a^3/mu);
t = [0:60:2*T];
for i = 1:numel(t)
    M = n*(t(i) - tau);    
    E = solvekep(M,norm(e));
    [r(i,:),v(i,:)] = calcRV(E,P,norm(e),ehat,ehatperp,mu);
end
figure
plot3(r(:,1),r(:,2),r(:,3))
figure
plot3(v(:,1),v(:,2),v(:,3))
% part ii
clear all
mu = 4e5; %Kg^3/s^2
e = [0 0.25 0.5 0.75 0.99];
rp = 10000;
i = 135;
Omega = 45;
omega = -90;
tau = 0;
nhat_Omega = cosd(Omega)*[1 0 0] + sind(Omega)*[0 1 0];
nhat_Omega_perp = -cosd(i)*sind(Omega)*[1 0 0] + cosd(i)*cosd(Omega)*[0 1 0] + sind(i)*[0 0 1];
ehat = cosd(omega)*nhat_Omega + sind(omega)*nhat_Omega_perp;
ehat_perp = -sind(omega)*nhat_Omega + cosd(omega)*nhat_Omega_perp;
figure
hold on
for i = 1:5
    P = rp*(1+e(i));
    a = P/(1 - e(i)^2);
    T = 2*pi*sqrt(a^3/mu);
    t = [0:60:2*T];
    n = sqrt(mu/a^3);
    for j = 1:numel(t)
        M = n*(t(j) - tau);
        E = solvekep(M,e(i));
        [R(j,:),~] = calcRV(E,P,e(i),ehat,ehat_perp,mu);
    end
    plot3(R(:,1),R(:,2),R(:,3))
end