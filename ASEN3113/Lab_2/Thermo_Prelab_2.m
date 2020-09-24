%%Prelab Question 3, determine distance x=0 to first thermocouple and H,
%%slope of steady state distribution
clear all;close all;clc
T_0 = 10; %Degress C, need to be in K?
T_x = [11.50,13.74,15.98,17.21,19.43,21.66,23.90,25.13];
x = linspace(0,.0889,numel(T_x));
P = polyfit(x,T_x,1);
y = P(1)*x + ((P(2)));
figure
hold on
scatter(x,T_x)
plot(x,y)
title('Given Data VS Linear Curve Fit')
xlabel('Thermocouple Position (m)')
ylabel('Temperature (C)')
legend('Given Data','Linear Curve Fit')
hold off
%Solve linear equation for when x = T_x, defined as x=0
x_0 = (10 - P(2))/(P(1));
x_0 = abs(x_0);%Meters
H = P(1); %Slope of steady state distribution
%%Prelab Question 5, plotting convergence for 0<=n<=10 for t = 1 and t =
%%1000
%H is same as in question 3
%Pull alpha for aluminium for table
alpha = 4.82e-5;
%L = distance from x=0 to heater, ~5 in
L = x_0 + .0889 + .0253;%distance to heater in m
%Location of last thermocouple
last_therm = x_0 + .0889;
%Plot convergence with variation in n at two t values
%t = 1
t = 1;
u_1 = 0;
for n = 1:10
    lambda = ((2*n - 1)*pi)/(2*L);
    b_n = ((8*H*L)*(-1)^n)/((2*n - 1)^2*(pi^2));
    current_u_1 = b_n*sin(lambda*last_therm)*exp((-lambda^2)*alpha*t);
    u_1(n+1) = u_1(n) + current_u_1;
end
u_1 = u_1(2:end);
u_1 = u_1 + T_0 + H*last_therm;
%t = 1000
t = 1000;
u_2 = 0;
for n = 1:10
    lambda = ((2*n - 1)*pi)/(2*L);
    b_n = ((8*H*L)*(-1)^n)/((2*n - 1)^2*(pi^2));
    current_u_2 = b_n*sin(lambda*last_therm)*exp((-lambda^2)*alpha*t);
    u_2(n + 1) = u_2(n) + current_u_2;
end
u_2 = u_2(2:end);
u_2 = u_2 + T_0 + H*last_therm;
n = [1:10];
figure
plot(n,u_1)
title('Temperature Convergence VS N Values at t=1')
xlabel('N Values')
ylabel('Temperature (C)')
figure
plot(n,u_2)
title('Temperature Convergence VS N Values at t=1000')
xlabel('N Values')
ylabel('Temperature (C)')
%%Question 6, plot temp over time with variation in k
n = 1;
u_a = 0;
for t = 1:1000
    lambda = ((2*n - 1)*pi)/(2*L);
    b_n = ((8*H*L)*(-1)^n)/((2*n - 1)^2*(pi^2));
    current_u = b_n*sin(lambda*last_therm)*exp((-lambda^2)*alpha*t);
    u_a(t + 1) = current_u;
end
u_a = u_a(2:end);
u_a = u_a + T_0 + H*last_therm;
u_b = 0;
alpha = .41e-5;
for t = 1:1000
    lambda = ((2*n - 1)*pi)/(2*L);
    b_n = ((8*H*L)*(-1)^n)/((2*n - 1)^2*(pi^2));
    current_u = b_n*sin(lambda*last_therm)*exp((-lambda^2)*alpha*t);
    u_b(t + 1) = current_u;
end
u_b = u_b(2:end);
u_b = u_b + T_0 + H*last_therm;
u_c = 0;
alpha = 3.56e-5;
for t = 1:1000
    lambda = ((2*n - 1)*pi)/(2*L);
    b_n = ((8*H*L)*(-1)^n)/((2*n - 1)^2*(pi^2));
    current_u = b_n*sin(lambda*last_therm)*exp((-lambda^2)*alpha*t);
    u_c(t + 1) = current_u;
end
u_c = u_c(2:end);
u_c = u_c + T_0 + H*last_therm;

t = 1:1000;
figure
hold on
plot(t,u_a)
plot(t,u_b)
plot(t,u_c)
title('Temperature VS Time for Variation in Diffusivity')
xlabel('Time (S)')
ylabel('Temperature (C)')
legend('Aluminum 7075-T651 \alpha = 4.82*10^-5 m^2/s','(Stainless Steel T-303 annealed \alpha = .41*10^-5 m^2/s','Brass C360 \alpha = 3.56*10^-5 m^2/s')
hold off