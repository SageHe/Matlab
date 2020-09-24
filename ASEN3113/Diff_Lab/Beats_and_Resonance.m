clear all;close all;clc

t = 0:.1:20;
y_h = cos(4*sqrt(3)*t);
y_p = sin(t);
y = cos(4*sqrt(3)*t) + sin(t);

subplot(2,1,1), plot(t,y_h,t,y_p)
subplot(2,1,2), plot(t,y)

y_p = t.*cos(4*sqrt(3)*t);
y = y_h + y_p;

figure
subplot(2,1,1), plot(t,y_h,t,y_p)
subplot(2,1,2), plot(t,y)
