clear all; close all; clc
%Compute and plot trajectory of rocket and calculate final height.
tspan = [0,18];
y0 = [0 0 .167];

[t,y] = ode45(@odefun,tspan,y0);

hold on
plot(t,y);
title('Velocity, Mass, and Heght Vs Time')
xlabel('Time')
ylabel('Velocity, Height, and Mass')
legend('Heigh','Velocity','Mass')

fprintf('Max height reached by the rocket is %.4f meters \n',max(y(:,1)))

fileID = fopen('C:\Users\Owner\Documents\MATLAB\ASEN_2012\Homework_5_Problem_2.txt','w');
fprintf(fileID,'Max height reached by the rocket is %.4f meters \n',max(y(:,1)));
fclose(fileID);