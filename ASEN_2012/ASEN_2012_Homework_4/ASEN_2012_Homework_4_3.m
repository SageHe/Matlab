clear all;close all;clc
%Sage Herrin, ASEN 2012, Due 10/26/17, 2PM lecture
%The purpose of this problem was to estimate the area under a curve
%using the midpoint rule with a given number of subintervals on a given
%interval and to speculate on the largest contributed error from one of the
%rectangles.

% Estimate the area under a curve described by f(x) = 3x^2 - x^2 + 2 using
% the midpoint rule over the interval [-2,2] with a spacing of 0.5 on the 
% x-axis. Plot the rectangles used to calculate area with the midpoint rule
% on the same figure as the smooth function f(x) to produce smooth curve for
% your figure. Use built-in matlab function 'rectangle' to plot
% rectangles. Which rectangle causes greatest source of error in area
% estimation?

% Calculate area under the curve using midpoint rule with a step size of
%0.5.

start = -2;
finish = 2;
delta_x = 0.5;
x = (start + (delta_x/2)):delta_x:(finish - (delta_x/2));
y = 3.*x.^3 - x.^2 + 2;
area = sum(y * delta_x);

% To determine number of calculations/iterations, determine number of
% subintervals n by manipulating delta_x - (b - a)/n => n = (b -
% a)/delta_x => n = 8 for the given step size and interval of interest
% x_approx = (-2 + (delta_x / 2));
% y_approx = [];
% n = (2 - -2)/delta_x;
% for i = 1:n
%     y_approx = [y_approx y(x_approx)];
%     x_approx = x_approx + delta(x);
% end
x_real = -2:.001:2;
y_real = 3.*x_real.^3 - x_real.^2 + 2;

figure
hold on
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

plot(x_real,y_real)
rectangle('Position',[-2,y(1),.5,y(1)*-1]);
rectangle('Position',[-1.5,y(2),.5,y(2)*-1]);
rectangle('Position',[-1.0,0,.5,y(3)]);
rectangle('Position',[-0.5,0,.5,y(4)]);
rectangle('Position',[0,0,.5,y(5)]);
rectangle('Position',[0.5,0,.5,y(6)]);
rectangle('Position',[1,0,.5,y(7)]);
rectangle('Position',[1.5,0,0.5,y(8)]);
title('3X^3 - X^2 + 2 Midpoint Rule')
xlabel('X')
ylabel('Y')
legend('Y = 3X^3 - X^2 + 2')

fileID = fopen('C:\Users\Owner\Documents\MATLAB\ASEN_2012\Homework_4_Problem_3.txt','w');
fprintf(fileID,'The estimated area under the curve using the midpoint rule on [-2,2] with a .5 spacing was %.2f square units. \nBy speculation of the plotted rectangles and the precise curve of the given function, it appears that \nthe leftmost rectangle is the gretest source of error since it seems that that rectangle has \nthe most area lying outside the precise curve of the function, resulting in the largest amount of area \nthat is not actually under the curve. This extra area is error in the estimation \n\nThe precise graph of the curve was obtained by plotting the given function against x, where \nx had a small stepsize of .001.',area);






