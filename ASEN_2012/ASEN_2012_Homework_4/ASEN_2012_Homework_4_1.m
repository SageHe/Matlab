clear all;close all;clc
%Sage Herrin, ASEN 2012, Due 10/26/17, 2PM lecture
%Purpose of this problem is to estimate the area under a curve using both
%simpsons rule and the trapezoid rule and observe the accuracy in the two
%methods used with different step sizes.

% Problem 1, determine the distance traveled by object after 4 seconds with
% accelertion a(t) = 6*t m/s^2, initial velocity of 1 m/s at t = 0 s, using
% trapezoid rule and Simpson's rule. Use equally spaced intervals of 10 and
% 100. Compare results of trapezoid rule, Simpson's rule, and analytic
% solution. Comment on performance of two approximation methods with uses
% of different interval sizes.

% Since distance traveled is the area under the curve of the velocity
% function, first we need a velocity function. Given accleration a(t) =
% 6*t, integration of this acceleration function gives us v(t) = 1+3*t^2,
% given the object has an initial velocity of 1 m/s at time t=0. Using new
% found velocity function, we can use the trapezoid rule and Simpson's rule
% with various partiton numbers to approximate the distance traveled.
% clear all;close all;clc
% Trapezoid rule
% 10 equally spaced subintervals

t = linspace(0,4,11);
V = 1 + 3*t.^2;

delta_x = t(10) - t(9);
y = [];
for i = 2:numel(t)
    y = [y (V(i - 1) + V(i))];
end
y = sum(y);
y1 = (delta_x / 2) * y;

fprintf('Distance traveled using trapezoid rule and 10 evenly spaced subintervals = %.4f meters \n\n',y1);

%Trapezoid rule with 100 equally spaced subintervals

t = linspace(0,4,101);
V = 1 + 3*t.^2;

delta_x = t(10) - t(9);
y = [];
for i = 2:numel(t)
    y = [y (V(i - 1) + V(i))];
end
y = sum(y);
y2 = (delta_x / 2) * y;

fprintf('Distance traveled using trapezoid rule and 100 evenly spaced subintervals = %.4f meters \n\n',y2);

% Simpson's rule
% 10 equally spaced subintervals

t = linspace(0,4,11); %Creates 10 equally spaced subintervals
V = 1 + 3*t.^2; %Velocity equation

delta_x = t(10) - t(9); %Finds the delta x to be used in simpson's rule calculation
y = [];    
for j = 1:(numel(t) - 1)/2
    y = [y (V(2*j - 1) + 4*V(2*j) + V(2*j + 1))*(delta_x/3)];
end
y3 = sum(y);

fprintf('Distance traveled using Simpsons rule and 10 evenly spaced subintervals = %.4f meters \n\n',y3);

% Simpson's rule
% 100 equally spaced subintervals

t = linspace(0,4,11); %Creates 10 equally spaced subintervals
V = 1 + 3*t.^2; %Velocity equation

delta_x = t(10) - t(9); %Finds the delta x to be used in simpson's rule calculation
y = [];    
for j = 1:(numel(t) - 1)/2
    y = [y (V(2*j - 1) + 4*V(2*j) + V(2*j + 1))*(delta_x/3)];
end
y4 = sum(y);

fprintf('Distance traveled using Simpsons rule and 100 evenly spaced subintervals = %.4f meters \n\n',y4);
fprintf('The analytic solution when done by hand results in an answer of exactly 68 meters \n');
fileID = fopen('C:\Users\Owner\Documents\MATLAB\ASEN_2012\Homework_4_Problem_1.txt','w');
fprintf(fileID,'Distance traveled using trapezoid rule and 10 evenly spaced subintervals = %.4f meters \n\n',y1);
fprintf(fileID,'Distance traveled using trapezoid rule and 100 evenly spaced subintervals = %.4f meters \n\n',y1);
fprintf(fileID,'Distance traveled using Simpsons rule and 10 evenly spaced subintervals = %.4f meters \n\n',y3);
fprintf(fileID,'Distance traveled using Simpsons rule and 10 evenly spaced subintervals = %.4f meters \n\n',y3);
fprintf(fileID,'The trapezoid rule gave a good approximation to the distance traveled by the object by using both 10 and 100 evenly spaced subintervals, \nhowever Simpsons rule gave a better approximation using both 10 and 100 evenly spaced subintervals since the Simpson approximation was exact in this case.\nLooking at the trapezoid approximation, using a larger number of evenly spaced subintervals results in a better approximation of the area under the curve,\nwhich can be seen by examining the difference between the exact answer and the approximated answer using the trapezoid rule with 10 \nevenly spaced subintervals vs 100 evenly spaced subintervals.');
fclose(fileID);


            
