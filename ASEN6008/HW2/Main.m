clear all;close all;clc
%% Problem 1 -- using Meeus algorithm and lambert solver to solve for v0 and vf for a given TOF to compute C3 and v_inf
mu = 1.32712440018e11;
AU = 1.49597870700e8;
Launch_JD = 2458239;
Arrival_JD = [2458389:1:2458639];
for i = 1:numel(Arrival_JD)
    EphemDepart = Meeus(Launch_JD);
    EphemArrive = Meeus(Arrival_JD(i));

    [r0,ve0] = calcposvel(EphemDepart.Earth.a,EphemDepart.Earth.e,EphemDepart.Earth.i,EphemDepart.Earth.Omega,EphemDepart.Earth.w,EphemDepart.Earth.nu,mu);
    [rf,vmf] = calcposvel(EphemArrive.Mars.a,EphemArrive.Mars.e,EphemArrive.Mars.i,EphemArrive.Mars.Omega,EphemArrive.Mars.w,EphemArrive.Mars.nu,mu);

    [v01,vf1] = solvelambert(r0,rf,(Arrival_JD(i)-Launch_JD)*86400,1);
    [v02,vf2] = solvelambert(r0,rf,(Arrival_JD(i)-Launch_JD)*86400,-1);
    
    vinf_mars1(i) = norm(vf1' - vmf');
    C31(i) = (norm(v01 - ve0))^2;
    vinf_mars2(i) = norm(vf2' - vmf');
    C32(i) = (norm(v02 - ve0))^2;
end
TOF = Arrival_JD - Launch_JD;
figure
hold on
vinf_min_ind = find(vinf_mars1 == min(vinf_mars1));
plot(TOF(vinf_min_ind),vinf_mars1(vinf_min_ind),'*')
plot(TOF(1:90),vinf_mars1(1:90))
xlabel('Days Since May 1, 2018')
ylabel('V_{\infty} at mars (km/s)')
grid on
grid minor
title('V_{\infty} vs TOF, Type 1')
ylim([2 6])
legend('Minimum V_{\infty}, Type 1')
figure
hold on
C3_min_ind = find(C31 == min(C31));
plot(TOF(C3_min_ind),C31(C3_min_ind),'*')
plot(TOF(1:85),C31(1:85))
xlabel('Days Since May 1, 2018')
ylabel('C3 (km^2/s^2)')
grid on
grid minor
title('C3 vs TOF, Type 1')
legend('Minimum C3 Value, Type 1')

figure
hold on
vinf_min_ind = find(vinf_mars2 == min(vinf_mars2));
plot(TOF(vinf_min_ind),vinf_mars2(vinf_min_ind),'*');
plot(TOF(100:120),vinf_mars2(100:120))
xlabel('Days Since May 1, 2018')
ylabel('V_{\infty} at mars (km/s)')
grid on
grid minor
title('V_{\infty} vs TOF, Type 2')
legend('Minimum V_{\infty}, Type 2')
figure
hold on
C3_min_ind = find(C32 == min(C32));
plot(TOF(C3_min_ind),C32(C3_min_ind),'*');
plot(TOF(100:120),C32(100:120))
xlabel('Days Since May 1, 2018')
ylabel('C3 (km^2/s^2)')
grid on
grid minor
title('C3 vs TOF, Type 2')
legend('Minimum C3 Value, Type 2')

%% Question 2
mu = 1.32712440018e11;
AU = 1.49597870700e8;

rp_d = 2.17*AU;
ra_d = 2.57*AU;
e_d = (ra_d - rp_d)/(ra_d + rp_d);
a_d = (rp_d + ra_d)/2;

P_c = 1682*86400;
e_c = 0.0758;
a_c = (((P_c/(2*pi))^2)*mu)^(1/3);
rp_c = a_c*(1 - e_c);
nu = a_d*(1 - e_d^2)/rp_c;
nu = nu - 1;
nu = nu/e_d;
nu = acos(nu);
E = 2*atan(sqrt((1-e_d)/(1+e_d))*tan(.5*nu));

dt = (E - e_d*sin(E))*sqrt(a_d^3/mu);
dt = dt/86400/365;

