clear all;close all;clc
Alpha = [-5:12 14 16];
C_L = [-0.32438 -0.21503 -0.10081 0.010503 0.12155 0.24163 0.34336 0.45256 0.56037 0.66625 0.76942 0.86923 0.96386 1.0441 1.0743 1.0807 1.0379 1.034 1.0156 0.97946];
C_D = [0.044251 0.033783 0.028627 0.025864 0.024643 0.025099 0.025635 0.02766 0.030677 0.034855 0.040403 0.04759 0.057108 0.070132 0.090921 0.11193 0.13254 0.15645 0.20959 0.25668];

% plot(Alpha,C_D)
% title('C_D')
% xlabel('Alpha')
% ylabel('C_D')
% figure
% plot(Alpha,C_L)
% title('C_L')
% xlabel('Alpha')
% ylabel('C_L')
% LOD = (C_L./C_D);
% figure
% plot(Alpha,LOD)
% title('Lift VS Drag')
% xlabel('Alpha')
% ylabel('Lift/Drag')

P = polyfit(Alpha(1:13),C_L(1:13),1);
x = [-5:7];
y = P(1).*x + P(2);

% figure
% hold on
% plot(Alpha,C_L)
% plot(x,y)

e = 0.601;
AR = 16.5;
C_D0 = C_D(5);

DP_C_D = C_D0 + ((C_L - C_L(5)).^2)./(pi*e*AR);

DP_C_D_lin = C_D0 + ((y - y(5)).^2./(pi*e*AR));

figure
hold on
plot(C_L(1:13),C_D(1:13))
plot(C_L(1:13),DP_C_D_lin)
title('Drag Polar')
xlabel('C_L')
ylabel('C_D')

[T,a,P,rho] = atmosisa(1829);

mass = 6.4; %Kg

weight = mass*9.81; %Newtons;

cruise_speed = 18; %m/s

S = 0.63; %m^2

cruise_C_L = (mass*9.81)/(.5*rho*(cruise_speed)^2*S);
cruise_C_D = C_D0 + (cruise_C_L - C_L(5))^2 / (pi*e*AR);

plot(cruise_C_L,cruise_C_D,'+')

n = .5; %Efficiency
P = 26.1*50*2; %Max power
P_a = n*P; %Power available


V = [12:55]; % m/s

C_l = weight./(.5*rho.*(V.^2)*S);

C_d = C_D0 + (C_l.^2)./(pi*e*AR);

P_r = sqrt((2*weight^3.*C_d.^2)./(rho*S.*C_l.^3));
P_a = P_a*ones(1,length(P_r));

figure
hold on
plot(V,P_r)
plot(V,P_a)
title('Power available VS Power Required')
xlabel('Velocity')
ylabel('Power')
legend('Power Required','Power Available')

Excess_power = (P_a - P_r);
ROC = Excess_power / weight;
max_EP = max(Excess_power);
max_ROC = max_EP / weight;

figure
plot(V,ROC)
title('Rate of Climb')
xlabel('Velocity')
ylabel('Rate of Climb')

mars_weight = 6.4*3.71; %Newtons
mars_rho = 0.02; %Kg/m^3

T_a = P_a ./ V;
TW = T_a ./ weight;
WS = weight / S;
% %% Design Lab 
% %Explore performance on mars based on paramter variations in TTwister such
% %as power avaialble, span, and wing area
% b = linspace(0.8*max(TW),1.2*max(TW));
% c = linspace(0.8*WS,1.2*WS);
% [B,C] = meshgrid(b,c);
% D = sqrt((B*C + C*sqrt(B^2 - (4*C_D0/pi*e*AR))/(rho*C_D0)));
% 
% figure
% contour(B,C,D,'ShowText','on')
% xlabel('TW')
% ylabel('WS')







