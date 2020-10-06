% ASEN 5052, HW 3, 
% Created by: Sage Herrin
clear all;close all;clc
%% Problem 2, part e -- plot ratio of deltaV_BP/deltaV_H for R = 1 -> 100 and ident. the value of R, denoted
% as R*, at which ratio equal 1
R = [1:100];
Ha = sqrt((2.*R)./(1+R));
Hb = sqrt(1./R);
Hc = sqrt((2./(R.*(1+R))));
deltaV_H = 2*pi*(Ha - 1 + Hb - Hc);
deltaV_BP = 2*pi.*(sqrt(2) - 1)*(1 + (1./sqrt(R)));
dv_rat = deltaV_BP./deltaV_H;
figure
plot(R,dv_rat)
xlabel('R')
ylabel('\Delta V_{BP}/\Delta V_H')
title('\Delta V ratio vs R Value')
Rs = interp1(dv_rat(2:end),R(2:end),1);
%% Problem 2, part f -- plot ratio of deltaV_BE/deltaV_H at a few values at and around the value of R*
%varying L in interval [R,100].
clear deltaV_BP R deltaV_H
R = linspace(Rs-4,Rs+5,10);
figure
hold on
for i = 1:numel(R)
    L = [R(i):100];
    Ha = sqrt((2.*R(i))./(1+R(i)));
    Hb = sqrt(1./R(i));
    Hc = sqrt((2./(R(i).*(1+R(i)))));
    deltaV_H = 2*pi*(Ha - 1 + Hb - Hc);
    BEa = sqrt((2.*L)./(1+L));
    BEb = sqrt((2.*R(i))./(L.*(R(i)+L)));
    BEc = sqrt((2./(L.*(1+L))));
    BEd = sqrt((2.*L)./(R(i).*(L+R(i))));
    BEe = sqrt(1./R(i));
    deltaV_BE = 2*pi*(BEa - 1 + BEb - BEc + BEd - BEe);
    dvrat = deltaV_BE./deltaV_H;
    plot(L,dvrat)
end
xlabel('L')
ylabel('\Delta V_{BE}/\Delta V_H')
title('\Delta V ratio vs L value')
%% Problem 2, part g 
% Compare cost of deltaV_BE and deltaV_BP by plotting their values across
%different ranges of R and L.
clear R L deltaV_BE deltaV_BP
R = linspace(1,20,100);
L = linspace(1,20,100);
for i = 1:numel(R)
    BEa = sqrt((2.*L)./(1+L));
    BEb = sqrt((2.*R(i))./(L.*(R(i)+L)));
    BEc = sqrt((2./(L.*(1+L))));
    BEd = sqrt((2.*L)./(R(i).*(L+R(i))));
    BEe = sqrt(1./R(i));
    deltaV_BE(i,:) = 2*pi*(BEa - 1 + BEb - BEc + BEd - BEe);
    deltaV_BP(i) = 2*pi.*(sqrt(2) - 1)*(1 + (1./sqrt(R(i))));
    vdiff(i,:) = deltaV_BE(i,:) - deltaV_BP(i);
end
figure
colormap hsv
contour(R,L,vdiff)
xlabel('R')
ylabel('L')
title('Countour of \Delta V_{BE}-\Delta V_{BP}')
%% Problem 3
% part d -- compute and compare the one-impulse cost with the 3-impulse
% cost for values of rho = 2,5,inf., starting with i = 60 degrees
mu = 4e5;
Re = 6400;
alt = 600;
rhoinf = @(r) 2*sqrt(mu/r)*(sqrt(2) - 1);
rho = [2 5];
dv360 = 2*sqrt(mu/(Re+alt)).*(sqrt((2.*rho)./(rho+1)) - 1 + sind(60/2)*sqrt(2./((1+rho).*rho)));
dv360 = [dv360 rhoinf(Re+alt)];
dv60 = 2*sqrt(mu/(Re+alt))*sind(60/2);
% repeat, but for a plane change of 45 degrees
mu = 4e5;
Re = 6400;
alt = 600;
rho = [2 5];
dv345 = 2*sqrt(mu/(Re+alt)).*(sqrt((2.*rho)./(rho+1)) - 1 + sind(45/2)*sqrt(2./((1+rho).*rho)));
dv345 = [dv345 rhoinf(Re+alt)];
dv45 = 2*sqrt(mu/(Re+alt))*sind(45/2);
% repeat, but for a plane change of 30 degrees
mu = 4e5;
Re = 6400;
alt = 600;
rho = [2 5];
dv330 = 2*sqrt(mu/(Re+alt)).*(sqrt((2.*rho)./(rho+1)) - 1 + sind(30/2)*sqrt(2./((1+rho).*rho)));
dv330 = [dv330 rhoinf(Re+alt)];
dv30 = 2*sqrt(mu/(Re+alt))*sind(30/2);
%% Problem 4, part c 
%plot ratio of deltaV_II/deltaV_I as func of V_inf/V_lc for several values
%of l on [0 1]
l = [0:0.1:1];
vinfvlcrat = linspace(0.5,5,11);
figure
hold on

for i = 1:numel(l)
    dvI(i,:) = (sqrt((vinfvlcrat).^2 + 2) - 1);
    dvII(i,:) = 1 - sqrt((2*l(i))./(1+l(i))) + sqrt((vinfvlcrat).^2 + (2./l(i))) - sqrt(2./(l(i)*(1+l(i))));
    rat(i,:) = dvII(i,:)./dvI(i,:);
    txt = ['l = ',num2str(l(i))];
    color = jet(length(l));
    plot(vinfvlcrat,rat(i,:),'DisplayName',txt,'color',color(i,:))
end
hold off
legend show
xlabel('V_{\infty}/V_{lc}')
ylabel('\Delta V_{II}/\Delta V_{I}')
title('Velocity ratio vs \Delta V ratio')