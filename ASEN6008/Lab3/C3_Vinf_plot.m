clear all;close all;clc
%Lab 3 - plot C3 and V_inf for trip durations between 150 and 250 days in
%10 day increments
TD = [150:10:250];

C3 = [20.486084 18.485843 17.145671 16.310663 15.892231 15.865466 16.275925 17.291178 19.327284 23.450093 32.973770];
V_inf = [5.299445 4.662140 4.125164 3.680125 3.322242 3.051535 2.872987 2.800007 2.860964 3.118345 3.734438];

figure
plot(TD,C3)
grid on
grid minor
xlabel('Transfer Duration (Days)')
ylabel('C3 $\frac{km^2}{s^2}$','Interpreter','latex')
title('C3 at Earth VS Transfer Duration','Interpreter','latex')

figure
plot(TD,V_inf)
grid on
grid minor
xlabel('Transfer Duration (Days)')
ylabel('$V_{\infty}$','Interpreter','latex')
title('$V_{\infty}$ at Mars VS Transfer Duration','Interpreter','latex')