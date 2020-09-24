clear all;close all; clc
set(0,'defaulttextInterpreter','latex');
%Read in text files of various materials 
Alum_280mA_data = load('Aluminum_30V_280mA_10DegC_1.txt'); %Aluminum data at 30 volts and 280 mA, 10C
Alum_300mA_data = load('Aluminum_30V_300mA_10DegC.txt'); %Aluminum data at 30 volts and 300 mA, 10C
Brass_190mA_data = load('Brass_20V_190mA_10DegC.txt'); %Brass at 20 volts and 190 mA, 10C
Brass_240mA_data = load('Brass_25V_240mA_8DegC.txt'); %Bass at 25 volts and 240 mA, 8C
Steel_300mA_data = load('Steel_30V_300mA_10DegC.txt'); %Steel at 30 volts and 300 mA, 10C

%% Alum_280_mA data set
T_0_alum_280 = mean(Alum_280mA_data(1,[2:end]));
%Thermocouples spaced by .5 inches -> total length of .0889 m
x = linspace(0,.0889,numel(Alum_280mA_data(1,[2:9])));
T_alum_280_steady = Alum_280mA_data(end,[2:end]);
P_A_280 = polyfit(x,T_alum_280_steady,1);
%Solve for x_0 distance using T_0 as avg value from zero time in data file
x_0_alum_280 = (T_0_alum_280 - P_A_280(2))/P_A_280(1);
% P_A_280(2) = P_A_280(2) - (x_0_alum_280);
Q_A_280 = 30*.280;
k_A_280 = 130;
diam = (2.54/100); %Diameter in meters
A = pi*(diam/2)^2;
H_an_A_280 = Q_A_280/(k_A_280*A); %Analytical slope, aluminum 
%% Alum_300mA data
T_0_alum_300 = mean(Alum_300mA_data(1,[2:end]));
x = linspace(0,.0889,numel(Alum_300mA_data(1,[2:end])));
T_alum_300_steady = Alum_300mA_data(end,[2:end]);
P_A_300 = polyfit(x,T_alum_300_steady,1);
x_0_alum_300 = (T_0_alum_300 - P_A_300(2))/P_A_300(1);
Q_A_300 = 30*.3;
k_A_300 = 130;
H_an_A_300 = Q_A_300/(k_A_300*A);
%% Brass_190mA data
T_0_brass_190 = mean(Brass_190mA_data(1,[2:end]));
x = linspace(0,.0889,numel(Brass_190mA_data(1,[2:end])));
T_brass_190_steady = Brass_190mA_data(end,[2:end]);
P_B_190 = polyfit(x,T_brass_190_steady,1);
x_0_brass_190 = (T_0_brass_190 - P_B_190(2))/P_B_190(1);
% P_B_190(2) = P_B_190(2) - x_0_brass_190;
Q_B_190 = 20*.190;
k_B_190 = 115;
H_an_B_190 = Q_B_190/(k_B_190*A);
%% Brass_240mA data
T_0_brass_240 = mean(Brass_240mA_data(1,[2:end]));
x = linspace(0,.0889,numel(Brass_240mA_data(1,[2:end])));
T_brass_240_steady = Brass_240mA_data(end,[2:end]);
P_B_240 = polyfit(x,T_brass_240_steady,1);
x_0_brass_240 = (T_0_brass_240 - P_B_240(2))/P_B_240(1);
Q_B_240 = 25*.240;
k_B_240 = 115;
H_an_B_240 = Q_B_240/(k_B_240*A);
%% Steel_300mA data
T_0_steel_300 = mean(Steel_300mA_data(1,[2:end]));
x = linspace(0,.0889,numel(Steel_300mA_data(1,[2:end])));
T_steel_300_steady = Steel_300mA_data(end,[2:end]);
P_S_300 = polyfit(x,T_steel_300_steady,1);
x_0_steel_300 = (T_0_steel_300 - P_S_300(2))/P_S_300(1);
Q_S_300 = 30*.3;
k_S_300 = 16.2;
H_an_S_300 = Q_S_300/(k_S_300*A);
%Plotting experimental and analytical slopes and experimental temp.
%distribution 
%Aluminum 280 mA plot
x = linspace(0,.0899,8);
x_A_280 = [x_0_alum_280 x];
T_alum_280_steady = [T_0_alum_280 T_alum_280_steady];
y_A_280_poly = polyval(P_A_280,x_A_280);
y_A_280 = H_an_A_280*x_A_280 + T_0_alum_280 + H_an_A_280*abs(x_0_alum_280);
x_A_280 = x_A_280 - x_0_alum_280;
figure(1)
grid on
grid minor
hold on
plot(x_A_280,y_A_280)
set(gca,'TickLabelInterpreter','latex');
plot(x_A_280,T_alum_280_steady,'o')
set(gca,'TickLabelInterpreter','latex');
plot(x_A_280,y_A_280_poly)
set(gca,'TickLabelInterpreter','latex');
title('Analytical VS Experimental SS Temp Distribution, Al, 280mA','fontsize',14)
xlabel('Distance From x_{0} [m]','Interpreter','tex')
ylabel('Temperature [C]')
legend('analytical temp distribution','experimental temp distribution','Location','northwest')
%Aluminum 300 mA plot
x_A_300 = [x_0_alum_300 x];
T_alum_300_steady = [T_0_alum_300 T_alum_300_steady];
y_A_300_poly = polyval(P_A_300,x_A_300);
y_A_300 = H_an_A_300*x_A_300 + T_0_alum_300 + H_an_A_300*abs(x_0_alum_300);
x_A_300 = x_A_300 - x_0_alum_300;
figure(2)
grid on
grid minor
hold on
plot(x_A_300,y_A_300)
set(gca,'TickLabelInterpreter','latex');
plot(x_A_300,T_alum_300_steady,'o')
set(gca,'TickLabelInterpreter','latex');
plot(x_A_300,y_A_300_poly);
set(gca,'TickLabelInterpreter','latex');
title('Analytical VS Experimental SS Temp Distribution, Al, 300mA','fontsize',14)
xlabel('Distance From x_{0} [m]','Interpreter','tex')
ylabel('Temperature [C]')
legend('analytical temp distribution','experimental temp distribution','Location','northwest')
% Brass 190 mA plot
x_b_190 = [x_0_brass_190 x];
T_brass_190_steady = [T_0_brass_190 T_brass_190_steady];
y_B_190_poly = polyval(P_B_190,x_b_190);
y_B_190 = H_an_B_190*x_b_190 + T_0_brass_190 + H_an_B_190*abs(x_0_brass_190);
x_b_190 = x_b_190 - x_0_brass_190;
figure(3)
grid on
grid minor
hold on
plot(x_b_190,y_B_190)
set(gca,'TickLabelInterpreter','latex');
plot(x_b_190,T_brass_190_steady,'o')
set(gca,'TickLabelInterpreter','latex');
plot(x_b_190,y_B_190_poly);
set(gca,'TickLabelInterpreter','latex');
title('Analytical VS Experimental SS Temp Distribution, Brass, 190mA','fontsize',14)
xlabel('Distance From x_{0} [m]','Interpreter','tex')
ylabel('Temperature [C]')
legend('analytical temp distribution','experimental temp distribution','Location','northwest')
%Brass 240 mA plot
x_b_240 = [x_0_brass_240 x];
T_brass_240_steady = [T_0_brass_240 T_brass_240_steady];
y_B_240_poly = polyval(P_B_240,x_b_240);
y_B_240 = H_an_B_240*x_b_240 + T_0_brass_240 + H_an_B_240*abs(x_0_brass_240);
x_b_240 = x_b_240 - x_0_brass_240;
figure(4)
hold on
grid on
grid minor
plot(x_b_240,y_B_240)
set(gca,'TickLabelInterpreter','latex');
plot(x_b_240,T_brass_240_steady,'o')
set(gca,'TickLabelInterpreter','latex');
plot(x_b_240,y_B_240_poly);
set(gca,'TickLabelInterpreter','latex');
title('Analytical VS Experimental SS Temp Distribution, Brass, 240mA','fontsize',14)
xlabel('Distance From x_{0} [m]','Interpreter','tex')
ylabel('Temperature [C]')
legend('analytical temp distribution','experimental temp distribution','Location','northwest')
%Steel 300 mA plot
x_s_300 = [x_0_steel_300 x];
T_steel_300_steady = [T_0_steel_300 T_steel_300_steady];
y_S_300_poly = polyval(P_S_300,x_s_300);
y_S_300 = H_an_S_300*x_s_300 + T_0_steel_300 + H_an_S_300*abs(x_0_steel_300);
x_s_300 = x_s_300 - x_0_steel_300;
figure(5)
grid on
grid minor
hold on
plot(x_s_300,y_S_300)
set(gca,'TickLabelInterpreter','latex');
plot(x_s_300,T_steel_300_steady,'o')
set(gca,'TickLabelInterpreter','latex');
plot(x_s_300,y_S_300_poly)
set(gca,'TickLabelInterpreter','latex');
title('Analytical VS Experimental SS Temp Distribution, Steel, 300mA','fontsize',14)
xlabel('Distance From x_{0} [m]','Interpreter','tex')
ylabel('Temperature [C]')
legend('analytical temp distribution','experimental temp distribution','Location','northwest')
%% Question 3 
%Alum 280mA T_0 distribution 
figure(6)
plot(x,Alum_280mA_data(1,[2:end]))
grid on
grid minor
set(gca,'TickLabelInterpreter','latex');
title('T_{0} Distribution, Al, 280mA','Interpreter','tex','FontSize',14)
xlabel('Distance From First Thermocouple [m]')
ylabel('Temperature [C]')
ylim ([11 13])
%Alum 300mA T_0 distribution
figure(7)
plot(x,Alum_300mA_data(1,[2:end]))
grid on
grid minor
set(gca,'TickLabelInterpreter','latex');
title('T_{0} Distribution, Al, 300mA','Interpreter','tex','FontSize',14)
xlabel('Distance From First Thermocouple [m]')
ylabel('Temperature [C]')
ylim ([11 13])
%Brass 190mA T_0 disribution 
figure(8)
plot(x,Brass_190mA_data(1,[2:end]))
grid on
grid minor
set(gca,'TickLabelInterpreter','latex');
title('T_{0} Distribution, Brass, 190mA','Interpreter','tex','FontSize',14)
xlabel('Distance From First Thermocouple [m]')
ylabel('Temperature [C]')
ylim ([11 12])
%Brass 240mA T_0 distribution
figure(9)
plot(x,Brass_240mA_data(1,[2:end]))
grid on
grid minor
set(gca,'TickLabelInterpreter','latex');
title('T_{0} Distribution, Brass, 240mA','Interpreter','tex','FontSize',14)
xlabel('Distance From First Thermocouple [m]')
ylabel('Temperature [C]')
ylim ([9 11])
%Steel 300mA T_0 distribution
figure(10)
plot(x,Steel_300mA_data(1,[2:end]))
grid on
grid minor
set(gca,'TickLabelInterpreter','latex');
title('T_{0} Distribution, Steel, 300mA','Interpreter','tex','FontSize',14)
xlabel('Distance From First Thermocouple [m]')
ylabel('Temperature [C]')
% ylim ([11 13])