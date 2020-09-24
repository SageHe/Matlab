%Roller coaster lab. Begins w/ initial height of 125m, starts with circ.
%hill, circ valley, loop, circ. valley, parab. hill, circ. valley, upward
%helix, banked turn@ constant height, circ. hill, straight line down, breaking section, end.

clear all;close all;clc

syms t
circle1_x = 0*t ;
circle1_y = cos(t) + 1;
circle1_z = sin(t) + 125;

fplot3(circle1_x,circle1_y,circle1_z,[0,pi]);

hold on

% circle2_x = 0*t;
% circle2_y = cos(t) + 1;
% circle2_z = sin(t);

xlabel('X')
ylabel('Y')
zlabel('Z')

% % fplot3(circle2_x,circle2_y,circle2_z,[0,pi/2]);

f = @(t) sqrt(((-sin(t)).^2) + (cos(t).^2));

length = integral(f,0,pi);