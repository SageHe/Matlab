% ASEN 2004
% Lab 1
% Aerodynamic Performance of the Tempest/TTwistor UAS

% Purpose: Analysis of TTwistor UAS.x

% Assumptions: Level unbalanced flight.

% Inputs: no.

% Outputs: OUTPUTS CALCULAtED PARAMETERS IN COMAND WINDOW! Also graphs the
% parameters for wing area inputed in variable S.

% Created: 01/19/2018
% Modified: 03/07/2018

%NOTES:
% - Values are outputed in the command window for easy access.

clear all
close all
clc;
%% INITIAL CONDITIONS:
JumpSize = 60000;

%Earth flight conditions:
Density = 1.02400617867;
g = 9.81;
S = 0.63;  %Earth wing area.

%Each battery weights 0.984, calculate weight of airplane:
W = (6.4)*g; % BASELINE WEIGTH

%Each battery has 50 amps, there are two batteries:
i = 100;  %[Amps] 

%Airplane conditions:
AR = 16.5;
ntot = 0.5;

%Battery conditions
Volts = 26.1;
Rt = 1;   %[hour]
C = 18;   %[Amps*hour]
n = 1.0;

%Minimum coefficients at alpha = 1 degree, where we have the minimum drag.
CLmin = 0.1216;
CDmin = 0.02464;
CD0 = 0.02511;

% All data from the lab instructions.
alpha = [ -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 14 16];
cl = [-0.32438 -0.21503 -0.10081 0.010503 0.12155 0.24163 0.34336 0.45256 0.56037 0.66625 0.76942 0.86923 0.96386 1.0441 1.0743 1.0807 1.0379 1.034 1.0156 0.97946 ];
cd = [0.044251 0.033783 0.028627 0.025864 0.024643 0.025099 0.025635 0.02766 0.030677 0.034855 0.040403 0.04759 0.057108 0.070132 0.090921 0.11193 0.13254 0.15645 0.20959 0.25668 ];

%Cut data of the linear part of Cl vs Alpha.
alpha2 = [ -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 ];
cl2 = [-0.32438 -0.21503 -0.10081 0.010503 0.12155 0.24163 0.34336 0.45256 0.56037 0.66625 0.76942 0.86923];
cd2 = [0.044251 0.033783 0.028627 0.025864 0.024643 0.025099 0.025635 0.02766 0.030677 0.034855 0.040403 0.04759];

%% INITIAL CALCULATIONS
%Get Oswald efficiency factor:
Oswald = 1.78*(1-0.045*AR.^(0.68))-0.64;

%Get k:
k = 1/(pi*AR*Oswald);

%Get the Drag coefficient for plot.:
CL = linspace(-0.4,1.2,5801);
%CD = CD0 + k*CL .^2;
CD = CDmin + k.*(CL - CLmin).^2;

%Fit linear portion of Cl:
coeff = polyfit(cl2,cd2,2);

Cd_fit = coeff(1).*CL.^2+coeff(2).*CL+coeff(3);

%Bsterry: get t.
t = Rt/(i^n)*(C/Rt)^n;

%Get the maximum power provided by the batteries.
Pmax = (Volts * C /Rt*(Rt/t)^(1/n))*0.5;

%Get the t for each element in the array
tg = Rt/(i^n)*(C/Rt)^n;

%Get array of discharges rates to see how power varies
ns = linspace(1,1.3); 

%% INITIALIZE:

%Initialize 1 dimensional arrays.
Stall_Speed = zeros(1,length(S));
vE = zeros(1,length(S));
vR = zeros(1,length(S));
Emax = zeros(1,length(S));
Rmax = zeros(1,length(S));
RCmax = zeros(1,length(S));
ws = zeros(1,length(S));
Vmax = zeros(1,length(S));
twmaxspeed = zeros(1,length(S));
 
%Initialize 2 dimensional arrays for plotting. 
Clspec = zeros(length(S),JumpSize);
Cdspec = zeros(length(S),JumpSize);
LDspec = zeros(length(S),JumpSize);
Power = zeros(length(S),JumpSize);
Excess = zeros(length(S),JumpSize);
RC = zeros(length(S),JumpSize);
E = zeros(length(S),JumpSize);
R = zeros(length(S),JumpSize);
tw = zeros(length(S),JumpSize);
v = zeros(length(S),JumpSize);

%% CALCULATIONS:
%Begin for loop. Each loop will calculate all parameters for each wing area
%in the initial array S.
for ii=1:length(S)
    
%Calculate the stall speed:
Stall_Speed(ii) = sqrt(2*W/(Density*S(ii)*max(1.3518)));

%Create array of velocities:
v(ii,:)  = linspace(Stall_Speed(ii),Stall_Speed(ii)+50,JumpSize);

%Get the required Cl and Cd from drag polar provided fit:
Clspec(ii,:) = 2*W./(Density.*v(ii,:).^2.*S(ii));
Cdspec(ii,:)  = CDmin + k.*(Clspec(ii,:) - CLmin).^2;  
LDspec(ii,:) = Clspec(ii,:)./Cdspec(ii,:);

%Get the required power based on speed:
Power(ii,:) = W.*v(ii,:)./(LDspec(ii,:));

%Calculate excess power and rate of climb:
Excess(ii,:) = Pmax-Power(ii,:);
RC(ii,:) = Excess(ii,:) ./ W;
RCmax(ii) = max(RC(ii,:));

%Calculate endurance for each velocity:
E(ii,:) = Rt.^(1-n)*(ntot.*Volts.*C./(Density.*v(ii,:).^3.*S(ii).*CD0./2+2.*W.^2.*k./(Density.*v(ii,:).*S(ii)))).^n;
%Calculate range for each velocity:
R(ii,:) = E(ii,:).*v(ii,:)*3.6;

%Calculate velocity for maximum endurance:
vE(ii) = sqrt((2.*W./(Density.*S(ii))).*sqrt(k./(3*CDmin)));
%Calculate velocity for maximum range:
vR(ii) = sqrt((2.*W./(Density.*S(ii))).*sqrt(k./CDmin));

%Calculate the max endurance:
Emax(ii) = Rt.^(1-n)*(ntot.*Volts.*C./((2./sqrt(Density.*S(ii))).*(CD0.^(1/4)).*((2.*W*sqrt(k./3)).^(3/2)))).^n;
%= max(E(ii,:));
%Calculate the max range:
Rmax(ii) = (Rt.^(1-n)*((ntot.*Volts.*C)./((1./sqrt(Density.*S(ii))).*(CD0.^(1/4)).*((2.*W*sqrt(k)).^(3/2)))).^n).*3.6.*vR(ii);
%= max(R(ii,:));

%Calculate Weight / surface area :
ws(ii) = W/S(ii);

%Find the maximum velocity
 u = find(Power(ii,:)<Pmax+.1 &  Power(ii,:)>Pmax-.1);
Vmax(ii) =  v(ii,u(1));

%Get thrust over weight:
tw(ii,:) = Pmax ./ v(ii,:);

%Get thrust over weight at the maximum speed:
twmaxspeed(ii) = tw(ii,u(1)) ;

end

%Calculate endurance for different discharge coefficients:
Rs = zeros(1,100);
for jj = 1:100
Rs(jj) = (Rt.^(1-ns(jj))*((ntot.*Volts.*C)./((1./sqrt(Density.*S(ii))).*(CD0^(1/4))*((2.*W*sqrt(k))^(3/2))))^ns(jj))*3.6*vR(ii);
end

%Calculate drag polar for ploting:
coefs=polyfit(alpha(1:13),cl(1:13),1);
line=polyval(coefs,alpha(1:13));

lowest=find(cd==min(cd));
C_dp=cd(lowest)+k*(line-line(lowest)).^2;

vv=18;
C_l_op=W*g/(.5*Density*vv^2*S);
C_d_op=cd(lowest)+k*(C_l_op-line(lowest)).^2;
C_dp2=.026+(cl.^2/(pi*Oswald*AR));

C_d_0=cd(lowest)+k*(0-line(lowest)).^2;

%% Print the information: 

fprintf('The wing areas are:')
S
fprintf('The wing loadings are:')
ws
fprintf('The t/w at max speed are:')
twmaxspeed
fprintf('The stall speeds are:')
Stall_Speed
fprintf('The maximum speeds are:')
Vmax
fprintf('The maximum endurances are:')
Emax
fprintf('The speed for maximum endurances are:')
vE
fprintf('The maximum ranges are:')
Rmax
fprintf('The speed for maximum ranges are:')
vR
fprintf('The maximum rates of climb are:')
RCmax

%% Plot:

%Create legend:
Leg=cell(length(S),1);
for pl1 = 1:length(S)
Leg{pl1} =strcat('Wing area =  ', num2str(S(pl1)),'  [m^2]');
end

%LIFT COEFFICIENT VS ALPHA
figure
plot(alpha,cl,'-x','LineWidth',1.1)
title('Coefficient of Lift')
xlabel('Angle of attack (deg)')
ylabel('C_L')
hold on
plot(alpha(1:13),line,'r','LineWidth',1.1)
grid on
legend('Experimental Data', 'Line of Best Fit')

figure
plot(alpha,cl./cd,'x-','LineWidth',1.1)
title('Lift to Drag Ratio')
xlabel('Angle of attack (deg)')
ylabel('C_l/C_d')
grid on

figure
plot(cl(1:12),cd(1:12),'x-','LineWidth',1.1)
title('Drag Polar')
xlabel('C_L')
ylabel('C_D')
xlim([0 1])
hold on
plot(cl(1:12),C_dp(1:12),'--','LineWidth',1.1)
plot(cl(1:12),C_dp2(1:12),'--','LineWidth',1.1)
plot(C_l_op,C_d_op,'go','LineWidth',1.2)
plot(0,C_d_0,'ro','LineWidth',1.2)
legend('Experimental','Theoretical Using C_D min','Theoretical Using C_D at 0 Lift','Operational Point','Zero-Lift','Location','SouthOutside')
grid on


% DRAGG COEFFICIENT VS ALPHA
figure 
plot(alpha, cd,'-x')
title('Drag coefficient vs Angle of attack')
xlabel('Angle of attack [Degrees]');
grid on
 set(findall(gca, 'Type', 'Line'),'LineWidth',1.2);
ylabel('Cd');


%RATE OF CLIMB VS VELOCITY
figure
hold on
grid on
for pl1 = 1:length(S)
plot(v(pl1,:),RC(pl1,:))
end
legend(Leg)
title('Rate of climb vs Velocity')
xlabel('Velocity [m/s]');
ylabel('Rate of climb [m/s]');
 set(findall(gca, 'Type', 'Line'),'LineWidth',1.2);
hold off

% %RANGE VS DISCHARGE PARAMETER
figure
hold on
grid on
plot(ns,Rs)
title('Maximum Range vs Discharge parameter')
xlabel('Discharge parameter');
ylabel('Range [Km]');
 set(findall(gca, 'Type', 'Line'),'LineWidth',1.2);
hold off

%ENDURANCE VS VELOCITY
figure
hold on
grid on
for pl1 = 1:length(S)
plot(v(pl1,:),E(pl1,:))
end
legend(Leg)
title('Endurance vs Velocity')
xlabel('Velocity [m/s]');
ylabel('Endurance [h]');
 set(findall(gca, 'Type', 'Line'),'LineWidth',1.2);
hold off
% 
%RANGE VS VELOCITY
figure
hold on
grid on
for pl1 = 1:length(S)
plot(v(pl1,:),R(pl1,:))
end
legend(Leg)
title('Range vs Velocity')
xlabel('Velocity [m/s]');
ylabel('Range [km]');
 set(findall(gca, 'Type', 'Line'),'LineWidth',1.2);
hold off


%POWERS VS VELOCITY
figure
hold on
grid on
for pl1 = 1:length(S)
plot(v(pl1,:),Power(pl1,:))
end
plot(v(1,:),Pmax*ones(1,length(v)))
Leg{length(S)+1} =strcat('Power available [w]');
legend(Leg)
title('Power required and available vs Velocity')
xlabel('Velocity [m/s]');
ylabel('Power [w]');
 set(findall(gca, 'Type', 'Line'),'LineWidth',1.2);
hold off




