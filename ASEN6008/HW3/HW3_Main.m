clear all;close all;clc
%%Problem 1 -- give position of Mars and Jupiter given Launch date of JD =
%%2454085.5 and TOF of 830 days to Jupiter;
mu = 1.32712440018e11;

Depart = 2454085.5;
EphemD = Meeus(Depart);
EphemA = Meeus(Depart + 830);
% pos. of Mars and Jupiter at departure date 
[rMarsD,vMarsD] = calcposvel(EphemD.Mars.a,EphemD.Mars.e,EphemD.Mars.i,EphemD.Mars.Omega,EphemD.Mars.w,EphemD.Mars.nu,mu);
[rJupiterD,vJupiterD] = calcposvel(EphemD.Jupiter.a,EphemD.Jupiter.e,EphemD.Jupiter.i,EphemD.Jupiter.Omega,EphemD.Jupiter.w,EphemD.Jupiter.nu,mu);
% pos. of Mars and Jupiter at arrival date
[rMarsA,vMarsA] = calcposvel(EphemA.Mars.a,EphemA.Mars.e,EphemA.Mars.i,EphemA.Mars.Omega,EphemA.Mars.w,EphemA.Mars.nu,mu);
[rJupiterA,vJupiterA] = calcposvel(EphemA.Jupiter.a,EphemA.Jupiter.e,EphemA.Jupiter.i,EphemA.Jupiter.Omega,EphemA.Jupiter.w,EphemA.Jupiter.nu,mu);
%%Problem 2 -- generate psi vs TOF plot for Mars-Jupiter transfer
psi0 = linspace(-4*pi,4*pi^2,100);
for i = 1:numel(psi0)
    [v0,vf,dt0(i)] = solvelambert_MR(rMarsD,rJupiterA,psi0(i),0,0);
end

psi1 = linspace(4*pi^2,16*pi^2,1000);
for i = 1:numel(psi1)
    [v0,vf,dt1(i)] = solvelambert_MR(rMarsD,rJupiterA,psi1(i),0,0);
end

hyplinet = zeros(1,numel(dt0));
revlinet = 4*pi^2*ones(1,numel(dt0));

figure
hold on
plot(psi0,dt0/86400)
plot(psi1,dt1/86400);
plot(hyplinet,dt0/86400,'--k')
plot(revlinet,dt0/86400,'--k')
% xlim([-20 150])
ylim([0 10000])
grid on
grid minor
xlabel('\Psi(rad^2)')
ylabel('TOF (Days)')
title('\Psi VS TOF, Mars to Jupiter Transfer')
legend('0 rev','1 rev')