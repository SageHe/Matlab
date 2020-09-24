clear all;close all;clc
V0 = 0;
ST = linspace(0,120,1000);
eta_th = 0.4;
h_PR = 18400;
g_c = 32.174;
TSFC = (ST.*g_c + 2*V0)/(2*eta_th*h_PR);
TSFC = TSFC/3600;

figure
grid on
grid minor 
plot(ST,TSFC)

V0 = 500;
TSFC = (ST.*g_c + 2*V0)/(2*eta_th*h_PR);
TSFC = TSFC/3600;

figure(2)
grid on
grid minor
plot(ST,TSFC)