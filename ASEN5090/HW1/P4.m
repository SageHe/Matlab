clear all;close all;clc
%Write code to simulate motion of aircraft specified in problem (const.
%altitude and speed,traveling from -x0 to +x0. Compute and plot ideal
%zenith angle, range, and range rate as functions of horizontal position of
%aircraft, x.

%Given/knowns;
v0 = 50; %m/s
x0 = 250; %m
h0 = 100; %m
%Create vector of horizontal ranges
x = linspace(-x0,x0,1000);
%Calculating ideal zenith angle, range, and range rate based off of slide
%17 of lecture 3
for i = 1:numel(x)
    theta(i) = (pi/2) - atan(h0/abs(x(i)));
    R(i) = sqrt(x(i)^2 + h0^2);
    Rdot(i) = (x(i)/(sqrt(x(i)^2 + h0^2)))*v0;
end
theta = theta*(180/pi);
figure
plot(x,theta)
grid on 
grid minor
title('Ideal Zenith Angle Vs X')
xlabel('X (m)')
ylabel('Zenith Angle (degs)')
figure
plot(x,R)
grid on 
grid minor
title('Aircraft Range Vs X')
xlabel('X (m)')
ylabel('Range (m)')
figure
plot(x,Rdot)
grid on 
grid minor
title('Range Rate Vs X')
xlabel('X (m)')
ylabel('Range Rate (m/s)')