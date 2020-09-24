tic
%% Function that reads in NACA airfoil paramters, chord length, panel
%number, and angle of attack, returning Cp plot,c_l, and plot of airfoil
clear all;close all;clc
[c_l,x1,x2,x3,alpha] = Vortex_Panel(0,0,15,1,1000,5);
%% Question 1 - plot Cp and output c_l for airfoil given input paramters
fprintf('The coefficient of lift for a NACA %d%d%d airfoil is %.4f at %d degrees angle of attack\n',x1, x2, x3, c_l, alpha)
%% Question 2 - Apply function to NACA 0015 at various angles of attack and perform quantitative error analysis
%Flow is computed at angle of attack of 5 degrees instead of 0 to see
%comparison of different resolutions better since c_l for a symmetric
%airfoil at 0 degree angle of attack is always 0. The nominal panel number
%is set to 1000 in the first iteration of the vortex panel.
N = [10:10:400];
C_L = [];
for i = 1:length(N)
    [C_L] = Vortex_Panel_alt(0,0,15,1,N(i),5);
    C_l(i) = C_L;
end
A = ones(1,length(N))*.95*c_l;
B = ones(1,length(N))*.99*c_l;
C = ones(1,length(N))*.999*c_l;
figure
hold on
plot(N,C_l)
plot(N,A)
plot(N,B)
plot(N,C)
grid on;grid minor
title(['C_L VS Number of Panels for NACA 0015 at \alpha=5',char(176)])
xlabel('Number of Panel')
ylabel('C_L')
legend('Running C_L','95% of Accepted C_L','99% of Accepted C_L','99.9% of Accepted C_L')
C_l_percent = ((c_l - C_l)/c_l)*100;
figure
plot(N,C_l_percent)
grid on;grid minor
title('Percent Error in C_L VS Number of Panels')
xlabel('Numer of Panels')
ylabel('Error[%]')
%Compute flow over NACA 0015 at varying angles of attack using nominal
%number of panels 
%-5 degrees 
alph = [-5 0 5 10];
C_L = [];
for i = 1:length(alph)
    [c_l,x1,x2,x3,alpha,Cp,x] = Vortex_Panel_alt(0,0,15,1,400,alph(i));
    C_L(i) = c_l;
end
figure
plot(alph,C_L)
grid on;grid minor
title('C_L VS angle of attack NACA 0015')
xlabel('\alpha [Degrees]')
ylabel('C_L')
[c_l,x1,x2,x3,alpha,Cp1,x] = Vortex_Panel_alt(0,0,15,1,400,-5);
[c_l,x1,x2,x3,alpha,Cp2,x] = Vortex_Panel_alt(0,0,15,1,400,0);
[c_l,x1,x2,x3,alpha,Cp3,x] = Vortex_Panel_alt(0,0,15,1,400,5);
[c_l,x1,x2,x3,alpha,Cp4,x] = Vortex_Panel_alt(0,0,15,1,400,10);
figure
hold on 
plot(x(1:end - 1),Cp1)
plot(x(1:end - 1),Cp2)
plot(x(1:end - 1),Cp3)
plot(x(1:end - 1),Cp4)
grid on;grid minor
hold off
axis ij
title('Cp of NACA 0015 varying \alpha')
xlabel('X/C')
ylabel('Cp')
legend('\alpha=-5 degrees','\alpha=0 degrees','\alpha=5 degrees','\alpha=10 degrees')
%% Question 3 - plot c_l vs alpha for 4 requested airfoils, find lift slope and zero-lift aoa, and compare against TAT
% For NACA 0015
alph = [-5:10];
for i = 1:length(alph)
    [c_l,x1,x2,x3,alpha,Cp,x] = Vortex_Panel_alt(0,0,15,1,400,alph(i));
    C_L(i) = c_l;
end
figure
hold on
plot(alph,C_L)
grid on;grid minor
title('C_L VS angle of attack')
xlabel('\alpha [Degrees]')
ylabel('C_L')
alph = alph.*(pi/180);
P1 = polyfit(alph,C_L,1);
%For NACA 2415
alph = [-5:10];
for i = 1:length(alph)
    [c_l,x1,x2,x3,alpha,Cp,x] = Vortex_Panel_alt(2,4,15,1,400,alph(i));
    C_L(i) = c_l;
end
% figure
plot(alph,C_L)
grid on;grid minor
title('C_L VS angle of attack')
xlabel('\alpha [Degrees]')
ylabel('C_L')
alph = alph.*(pi/180);
P2 = polyfit(alph,C_L,1);
%For NACA 4415
alph = [-5:10];
for i = 1:length(alph)
    [c_l,x1,x2,x3,alpha,Cp,x] = Vortex_Panel_alt(4,4,15,1,400,alph(i));
    C_L(i) = c_l;
end
% figure
plot(alph,C_L)
grid on; grid minor
title('C_L VS angle of attack')
xlabel('\alpha [Degrees]')
ylabel('C_L')
alph = alph.*(pi/180);
P3 = polyfit(alph,C_L,1);
%For NACA 2424
alph = [-5:10];
for i = 1:length(alph)
    [c_l,x1,x2,x3,alpha,Cp,x] = Vortex_Panel_alt(2,4,24,1,400,alph(i));
    C_L(i) = c_l;
end
% figure
x_ax = linspace(0,0,length(alph));
y_ax = linspace(-1,2,length(alph));
hold on
plot(alph,C_L)
plot(alph,x_ax,'--k','LineWidth',2)
plot(x_ax,y_ax,'--k','LineWidth',2)
grid on; grid minor 
title('C_L VS angle of attack')
xlabel('\alpha [Degrees]')
ylabel('C_L')
legend('NACA 0015','NACA 2415','NACA 4415','NACA 2424')
alph = alph.*(pi/180);
P4 = polyfit(alph,C_L,1);
hold off
%Solving for alpha zero lift for each airfoil
%NACA 0015
alpha_0_L_0015 = -P1(2)/P1(1)*(180/pi);
%NACA 2415
alpha_0_L_2415 = -P2(2)/P2(1)*(180/pi);
%NACA 4415
alpha_0_L_4415 = -P3(2)/P3(1)*(180/pi);
%NACA2424
alpha_0_L_2424 = -P4(2)/P4(1)*(180/pi);

toc