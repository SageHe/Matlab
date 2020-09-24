clear all;close  all;clc
%define zeta, w_n, and w_d to be used 
zeta = 0.5411;
w_n = 3.383;
w_d = w_n*sqrt(1 - zeta^2);

t = linspace(0,5,1000);

x = (1 - exp(-zeta.*w_n.*t).*(cos(w_d.*t) + (zeta./sqrt(1 - zeta^2)).*sin(w_d.*t)))*.5;

line = linspace(.5,.5,1000);

figure
hold on
grid on
grid minor 
plot(t,x);
plot(t,line,'--')
title('Modeled Response')
xlabel('Time (s)')
ylabel('Radians')
legend('Expected Response','Step Input')