clear all
close all
clc
%Set the time span for ODE to run over 
tspan = [0 20];
%define given m paramter
m = 1;
%Define given k parameter
k = 1;
%Define w_n = sqrt(k/m)
w_n = sqrt(k/m);
%Define initial conditions of the system
x_0 = 4; % Initial displacement [m]
xdot_0 = 0; % Initial velocity [m/s]
inish_condish = [x_0 xdot_0];

%Call ODE45 to numerically solve the system
[t,y] = ode45('odefun',tspan,inish_condish);

%Plot to displacement and velocity of the mass over time 
figure
subplot(2,1,1)
hold on
plot(t(:,1),y(:,1))
plot(t,y(:,2))
xlabel('Time (s)')
ylabel('Displacement or Velocity')
title ('Displacement and Velocity VS Time (numerical solution)')
legend('Mass Position (m)','Mass Velocity (m/s)')
subplot(2,1,2)
hold on
plot(t,x_0*cos(w_n*t));
plot(t,x_0*-sin(w_n*t));
xlabel('Time (s)')
ylabel('Displacement or Velocity')
title('Displacement and Velocity VS Time (analytical solution')
xlabel('Time (s)')
legend('Mass Position (m)','Mass Velocity (m/s)')
