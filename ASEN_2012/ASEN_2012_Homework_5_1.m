clear all;close all;clc
%Numerical Integration using Euler Explicit Method
x = 0:.1:2;
delta_x = .1;
y1 = [];
y1(1) = 1;
for i = 1:length(x) - 1
    y1(i + 1) = y1(i) + delta_x * ((x(i)*y1(i) + 3*x(i)));
end

%Numerical Integration using Runge Kutta method
y2 = [];
y2(1) = 1;
for j = 1:length(x) - 1
    k1 = (x(j)*y2(j) + 3*x(j));
    k2 = ((x(j) + (delta_x / 2))*(y2(j) + (.5*k1*delta_x)) + (3*(x(j) + (.5*delta_x))));
    k3 = ((x(j) + (delta_x / 2))*(y2(j) + (.5*k2*delta_x)) + (3*(x(j) + (.5*delta_x))));
    k4 = ((x(j) + delta_x)*(y2(j) + (k3 * delta_x)) + (3*(x(j) + delta_x)));
    y2(j + 1) = y2(j) + (delta_x / 6)*(k1 + 2*k2 + 2*k3 + k4);
end

analytical = 4 * exp((x.^2) / 2) - 3;

hold on
plot(x,y1)
plot(x,y2,'*')
plot(x,analytical)
title('Euler Vs Runge Kutta Vs Analytical Solutions')
legend('Euler','Runge Kutta','Analytical')
xlabel('X')
ylabel('Y')

fprintf('Euler Explicit solution yields y = %.6f @ x = 2 \n', y1(end))
fprintf('Runge Kutta solution yields y = %.6f @ x = 2 \n', y2(end))
fprintf('Analytical solution yields y = %.6f @ x = 2 \n', analytical(end))

fileID = fopen('C:\Users\Owner\Documents\MATLAB\ASEN_2012\Homework_5_Problem_1.txt','w');
fprintf(fileID,'Euler Explicit solution yields y = %.6f @ x = 2 \n', y1(end));
fprintf(fileID,'Runge Kutta solution yields y = %.6f @ x = 2 \n', y2(end));
fprintf(fileID,'Analytical solution yields y = %.6f @ x = 2 \n', analytical(end));
fclose(fileID);