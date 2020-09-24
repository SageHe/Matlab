%Cooling Sphere Fourier Project
%Author: Sage Herrin
%Created: 12/3/2018
%Modified: 12/3/2018

clear all; close all; clc
ice_bath_large_data = load('Ice_bath_large.txt');
ice_bath_small_data = load('Ice_bath_small.txt');
LN2_bath_small = load('LN2_bath_small.txt');
%Small stainless steel sphere
t = ice_bath_small_data(:,1);
alt_temp_data_small_bath = ice_bath_small_data(:,[2 3]);
alt_temp_data_small_bath = log(alt_temp_data_small_bath);

figure(1)
hold on
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,2))
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,3))
title('Original Small Stainless Steel Temperature Data')
xlabel('Time(S)')
ylabel('Temperature(C)')
legend('r=0','r=R/2')

figure(2)
hold on
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,2))
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,3))
title('Original Large Chrome Steel Temperature Data')
xlabel('Time(S)')
ylabel('Temperature(C)')
legend('r=0','r=R/2')

figure(3)
hold on
plot(t,alt_temp_data_small_bath(:,1))
plot(t,alt_temp_data_small_bath(:,2))
title('Log plot of small sphere data')
xlabel('Time(s)')
ylabel('Temp(C), Log Scale')
legend('r=0','r=R/2')
% plot(t,ice_bath_large_data(:,2))
% plot(t,ice_bath_large_data(:,3))

%Looking at linear sections of interest, outer thermocouple, time
%t=20-58.75 seconds
t_temp = ice_bath_small_data([81:236],1);
Temp_temp = log(ice_bath_small_data([81:236],2));
P1 = polyfit(t_temp,Temp_temp,1);
A1 = polyval(P1,t_temp);
figure(4)
hold on
plot(t_temp,Temp_temp)
plot(t_temp,A1)
title('Best fit plot small sphere')
xlabel('Time(s)')
ylabel('Temp(C), log scale')
legend('Data','Best Fit')
R_small = .0238; %radius of small sphere (stainless steel) in meters
K_SS = abs(P1(1)/((pi/R_small)^2)); %Thermal diffusivity stainless steel (SS) small sphere

%% Large chrome steel sphere
t = ice_bath_large_data(:,1);
alt_temp_data_large_bath = ice_bath_large_data(:,[2 3]);
alt_temp_data_large_bath = alt_temp_data_large_bath;
alt_temp_data_large_bath = log(alt_temp_data_large_bath);

figure(5)
hold on
plot(t,alt_temp_data_large_bath(:,1))
plot(t,alt_temp_data_large_bath(:,2))
title('Log plot of large sphere data')
xlabel('Time(s)')
ylabel('Temp(C), Log Scale')
legend('r=0','r=R/2')
t = ice_bath_large_data(:,1);
alt_temp_data_large_bath = ice_bath_large_data(:,[2 3]);
alt_temp_data_large_bath = log(alt_temp_data_large_bath);

%Looking at linear sections of interest, outer thermocouple, t = 18.75-75.25 seconds, 
t_temp = ice_bath_small_data([81:236],1);
Temp_temp = log(ice_bath_small_data([81:236],2));
P2 = polyfit(t_temp,Temp_temp,1);
A2 = polyval(P2,t_temp);
figure(6)
hold on
plot(t_temp,Temp_temp)
plot(t_temp,A2)
title('Best fit plot large sphere')
xlabel('Time(s)')
ylabel('Temp(C), log scale')
legend('Data','Best Fit')
R_large = .03125; %radius of small sphere (stainless steel) in meters
K_CS = abs(P2(1)/((pi/R_large)^2)); %Thermal diffusivity chrome steel (CS) large sphere

%Plotting a couple of modes of temp at r=0 and r = R/2 
%large sphere r = R/2
u = 0;
u_plot = 0;
t = 0:300;
n = 1:3;
k_large = 3.9843e-6;
R = .03125;%radius in meters
r = R/2;
for i = 1:length(t)
    u = 0;
    for j = 1:length(n)
        lambda = exp((-k_large*((n(j)*pi)/R)^2)*t(i));
%         u(j + 1) = u(j) + ((-44.5*R)/(pi*n(j))*cos(n(j)*pi)*lambda*sin((n(j)*pi*r)/R));
        u = u + ((-46.6*R)/(pi*n(j))*cos(n(j)*pi)*lambda*sin((n(j)*pi*r)/R));
    end
%     u = u(2:end);
    u = .2 + (1/r)*u;
    u_plot(i) = u;
end

figure(8)
hold on
plot(t,u_plot)
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,2))
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,3))
title('Numerical Approximation VS Data, Large Sphere r=R/2')
xlabel('Time(S)')
ylabel('Temperature(C)')
legend('Numerical model r=r/2','r=0','r=R/2')
%large sphere at r=0
u = 0;
u_plot = 0;
t = 0:300;
n = 1:30;
k_large = 3.3843e-6;
R = .03125;%radius in meters
r = 0;
for i = 1:length(t)
    u = 0;
    for j = 1:length(n)
        lambda = exp((-k_large*((n(j)*pi)/R)^2)*t(i));
%         u(j + 1) = u(j) + ((-44.5*R)/(pi*n(j))*cos(n(j)*pi)*lambda*sin((n(j)*pi*r)/R));
        u = u + (-46.6)/(n(j))*cos(n(j)*pi)*lambda;
    end
%     u = u(2:end);
    u = .2 + u;
    u_plot(i) = u;
end

figure(9)
hold on
plot(t,u_plot)
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,2))
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,3))
title('Numerical Approximation VS Data, Large Sphere r=0')
xlabel('Time(S)')
ylabel('Temperature(C)')
legend('Numerical Model r=0,k adjusted to 3.3843e-5','r=0','r=R/2')
%small sphere at r=R/2
u = 0;
u_plot = 0;
t = 0:300;
n = 1:3;
k_small = 2.811e-6; %2.311e-6
R = .0238125;%radius in meters
r = R/2;
for i = 1:length(t)
    u = 0;
    for j = 1:length(n)
        lambda = exp((-k_small*((n(j)*pi)/R)^2)*t(i));
%         u(j + 1) = u(j) + ((-44.5*R)/(pi*n(j))*cos(n(j)*pi)*lambda*sin((n(j)*pi*r)/R));
        u = u + ((-46.6*R)/(pi*n(j))*cos(n(j)*pi)*lambda*sin((n(j)*pi*r)/R));
    end
%     u = u(2:end);
    u = .2 + (1/r)*u;
    u_plot(i) = u;
end

figure(10)
hold on
plot(t,u_plot)
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,2))
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,3))
title('Numerial Model VS Data, Small Sphere R=r/2')
xlabel('Time(S)')
ylabel('Temperature(C)')
legend('Numerical Model r=r/2, k adjusted to 2.811e-5','r=0','r=R/2')
%small sphere r=0
u = 0;
u_plot = 0;
t = 0:300;
n = 1:3;
k_small = 2.311e-6;
R = .0238125;%radius in meters
r = 0;
for i = 1:length(t)
    u = 0;
    for j = 1:length(n)
        lambda = exp((-k_small*((n(j)*pi)/R)^2)*t(i));
%         u(j + 1) = u(j) + ((-44.5*R)/(pi*n(j))*cos(n(j)*pi)*lambda*sin((n(j)*pi*r)/R));
        u = u + (-46.6)/(n(j))*cos(n(j)*pi)*lambda;
    end
%     u = u(2:end);
    u = .2 + u;
    u_plot(i) = u;
end

figure(11)
hold on
plot(t,u_plot)
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,2))
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,3))
title('Numerical Model VS Data, Small Sphere R=0')
xlabel('Time(S)')
ylabel('Temperature(C)')
legend('Numerical Model r=0','r=0','r=R/2')
%Plot of two different regimes of LN2 bath
figure(7)
hold on
plot(LN2_bath_small(:,1),LN2_bath_small(:,2))
plot(LN2_bath_small(:,1),LN2_bath_small(:,3))
title('LN2 bath small sphere')
xlabel('Time(s)')
ylabel('Temp(C)')
legend('r=R/2','r=0')
%LN2 data split at t = 1000
t_temp = LN2_bath_small(:,1);
Temp_temp = log(abs(LN2_bath_small(:,[2:3])));
figure(12)
hold on
plot(t_temp,Temp_temp(:,1))
plot(t_temp,Temp_temp(:,2))
title('Log Plot of Small Sphere in LN2 Bath')
xlabel('Time(S)')
ylabel('Temperature(C),Log Scale')
legend('r=0','r=R/2')
%First regime        
t_temp = LN2_bath_small([365:1000],1);
Temp_temp = log(abs(LN2_bath_small([365:1000],2)));
P3 = polyfit(t_temp,Temp_temp,1);
A3 = polyval(P3,t_temp);
figure(13)
hold on
plot(t_temp,Temp_temp)
plot(t_temp,A3)
title('Best fit plot small sphere first region LN2')
xlabel('Time(s)')
ylabel('Temp(C), log scale')
legend('Data','Best Fit')
R_small = .0238125; %radius of small sphere (stainless steel) in meters
K_SS_LN2_R1 = abs(P3(1)/((pi/R_small)^2)); %Thermal diffusivity chrome steel (CS) large sphere
%Second Regime
t_temp = LN2_bath_small([1041:1101],1);
Temp_temp = log(abs(LN2_bath_small([1041:1101],2)));
P4 = polyfit(t_temp,Temp_temp,1);
A4 = polyval(P4,t_temp);
figure(14)
hold on
plot(t_temp,Temp_temp)
plot(t_temp,A4)
title('Best fit plot small sphere second region LN2')
xlabel('Time(s)')
ylabel('Temp(C), log scale')
legend('Data','Best Fit')
R_small = .0238125; %radius of small sphere (stainless steel) in meters
K_SS_LN2_R2 = abs(P4(1)/((pi/R_small)^2)); %Thermal diffusivity chrome steel (CS) large sphere

figure(15)
hold on
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,2))
plot(ice_bath_large_data(:,1),ice_bath_large_data(:,3))
time_line = linspace(-5,10.0,numel(ice_bath_large_data(:,1)));
x1_temp = 30.5*ones(1,length(ice_bath_large_data(:,1)));
plot(x1_temp,time_line,'--r')
x3_temp = 45.25*ones(1,length(ice_bath_large_data(:,1)));
time_line_2 = linspace(-5,10,numel(ice_bath_large_data(:,1)));
x4_temp = linspace(0,45.25,numel(ice_bath_large_data(:,1)));
temp_line_2 = 10.0*ones(1,numel(ice_bath_large_data(:,1)));
plot(x3_temp,time_line_2,'--g')
plot(x4_temp,temp_line_2,'--b')
title('Temperature of Large Sphere Over Time')
xlabel('Time(S)')
ylabel('Temperture(C)')
legend('r=0','r=R/2')

figure(16)
hold on
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,2))
plot(ice_bath_small_data(:,1),ice_bath_small_data(:,3))
x1_temp = 36.5*ones(1,length(ice_bath_small_data(:,1)));
time_line = linspace(-5,10,numel(ice_bath_small_data(:,1)));
plot(x1_temp,time_line,'--r')
x2_temp = 22*ones(1,length(ice_bath_small_data(:,1)));
time_line_2 = linspace(-5,10,numel(ice_bath_small_data(:,1)));
plot(x2_temp,time_line_2,'--g')
x3_temp = linspace(0,36.5,numel(ice_bath_small_data(:,1)));
temp_line =10.0*ones(1,numel(ice_bath_small_data(:,1)));
plot(x3_temp,temp_line,'--b')
title('Temperature of Small Sphere Over Time')
xlabel('Time(S)')
ylabel('Temperture(C)')
legend('r=0','r=R/2')
