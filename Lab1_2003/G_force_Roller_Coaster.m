% Charles Puskar
% ASEN 2003
% Lab 1

clear all; close all; clc;

% Date Created: January 17, 2017
% Date Modified: January 24, 2017

% initial conditions
g = 9.807; % m/s^2
h0 = 125; % m
z_current = h0;

% defines symbolic parameter t
syms t

% defines velocity from current height
v = @(z) sqrt(2*g*(h0 - z));

%% PARAMETERIZING THE PATH

figure(1)

% circle 1
r1 = 12.5;

circle1_x = 0*t;
circle1_y = r1*sin(t);
circle1_z = r1*cos(t)+z_current-r1;

x_current = 0;
y_current = r1;
z_current = z_current - r1;

fplot3(circle1_x,circle1_y,circle1_z,[0,pi/2])

xlabel('x')
ylabel('y')
zlabel('z')

axis equal

hold on
%%
%circle 1 g force

i = 1
height0 = 125
% for t = 0:.1:(pi/2)
%     vel(i,:) = r1*-sin(t);
%     G(i,:) = -(vel(i,:)^2/(g*r1)) + cos(t)
%     z_current = z_current - r1;
%     i = i+1
time = 0
for height = h0:-.1:z_current
    if time <= pi/2
    vel = sqrt(2*g*(h0 - height));
    G1(i,:) = -(vel^2/(g*r1)) + cos(time)
    i = i+1
    time = time + .1
    end
end


% cicle 2
r2 = 12.5;

circle2_x = 0*t + x_current;
circle2_y = -r2*cos(t)+y_current+r2;
circle2_z = -r2*sin(t)+z_current;

y_current = y_current+r2;
z_current = z_current - r2;

fplot3(circle2_x,circle2_y,circle2_z,[0,pi/2])

%%
%circle2 g force
i = 1
time = pi/2
for height = (h0-r1):-.1:z_current
    if time <= (3*pi)/2
    vel = sqrt(2*g*(h0 - height));
    G2(i,:) = +(vel^2/(g*r1)) - cos(time)
    i = i+1
    time = time + .1
    end
end


%%


% loop 1
r3 = 10;

l1_x = 0*t;
l1_y = r3*t+r1+r2;
l1_z = z_current + 0*t;

fplot3(l1_x,l1_y,l1_z,[0,1])

circle3_x = 0*t;
circle3_y = r3*sin(t)+y_current+r3;
circle3_z = -r3*cos(t)+z_current+r3;

y_current = y_current+r3;

fplot3(circle3_x,circle3_y,circle3_z,[0,2*pi])

l2_x = 0*t;
l2_y = r3*t+y_current;
l2_z = z_current + 0*t;

y_current = y_current + r3;

fplot3(l2_x,l2_y,l2_z,[0,1])

% circle 3
r4 = 10;

circle4_x = 0*t;
circle4_y = r4*sin(t)+y_current;
circle4_z = -r4*cos(t)+z_current+r4;

theta1 = pi/4;

y_current = y_current + r4*sin(theta1);
z_current = z_current + r4*(1-cos(theta1));

fplot3(circle4_x,circle4_y,circle4_z,[0,theta1])

% parabola 1

parabola1_x = 0*t;
parabola1_y = v(z_current)*cos(theta1)*t+y_current;
parabola1_z = @(t) v(z_current)*sin(theta1)*t - 0.5*g*t.^2 + z_current;

find_max = @(t) diff(parabola1_z(t),t);
t_crit = double(solve(find_max(t)));

y_current = y_current+v(z_current)*cos(theta1)*2*t_crit;

fplot3(parabola1_x,parabola1_y,parabola1_z,[0,t_crit*2])

% circle 4
r5 = 10;

circle5_x = 0*t;
circle5_y = -r5*sin(t)+y_current+r5*sin(theta1);
circle5_z = -r5*cos(t)+z_current-r5*(1-cos(theta1))+r5;

y_current = y_current + r5*sin(theta1);
z_current = z_current - r5*(1-cos(theta1));

fplot3(circle5_x,circle5_y,circle5_z,[0,theta1])

% helix 1
r6 = 8;

l3_x = 0*t;
l3_y = r6*t+y_current;
l3_z = z_current + 0*t;

fplot3(l3_x,l3_y,l3_z,[0,1])

helix1_x = -r6*cos(t)+r6;
helix1_y = r6*sin(t) + y_current+r6;
helix1_z = z_current + 20*t/(4.5*pi);

x_current = x_current + r6;
y_current = y_current + 2*r6;
z_current = z_current + 20;

fplot3(helix1_x,helix1_y,helix1_z,[0,4.5*pi]);

l4_x = r6*t + r6;
l4_y = 0*t +y_current;
l4_z = z_current + 0*t;

x_current = x_current + r6;

fplot3(l4_x,l4_y,l4_z,[0,1])

% linear 1
r7 = 10;
theta2 = pi/6;

circle6_x = r7*sin(t)+x_current;
circle6_y = 0*t + y_current;
circle6_z = r7*cos(t)+z_current-r7;

x_current = x_current + r7*sin(theta2);
z_current = z_current - r7*(1-cos(theta2));

fplot3(circle6_x,circle6_y,circle6_z,[0,theta2])

drop1 = 40 - 2*r7*(1-cos(theta2));

l5_x = cot(theta2)*t + x_current;
l5_y = 0*t +y_current;
l5_z = -t+z_current;

x_current = x_current + cot(theta2)*drop1;
z_current = z_current - drop1;

fplot3(l5_x,l5_y,l5_z,[0,drop1])

r8 = 10;

circle7_x = -r8*sin(t)+x_current + r8*sin(theta2);
circle7_y = 0*t + y_current;
circle7_z = -r8*cos(t)+z_current +r8 - r8*(1-cos(theta2));

x_current = x_current + r7*sin(theta2);
z_current = z_current - r7*(1-cos(theta2));

fplot3(circle7_x,circle7_y,circle7_z,[0,theta2])

% banked turn
r9 = 50;

circle8_x = r9*sin(t)+x_current;
circle8_y = -r9*cos(t)+y_current-r9;
circle8_z = 0*t + z_current;

y_current = y_current - 2*r9;

fplot3(circle8_x,circle8_y,circle8_z,[0,pi])

% linear 2
r10 = 10;
theta3 = pi/6;

circle9_x = -r10*sin(t)+x_current;
circle9_y = 0*t + y_current;
circle9_z = r10*cos(t)+z_current-r10;

x_current = x_current - r10*sin(theta3);
z_current = z_current - r10*(1-cos(theta3));

fplot3(circle9_x,circle9_y,circle9_z,[0,theta3])

drop2 = 80 - 2*r10*(1-cos(theta3));

l6_x = -cot(theta3)*t + x_current;
l6_y = 0*t +y_current;
l6_z = -t+z_current;

x_current = x_current - cot(theta3)*drop2;
z_current = z_current - drop2;

fplot3(l6_x,l6_y,l6_z,[0,drop2])

r11 = 10;

circle10_x = r11*sin(t)+x_current - r11*sin(theta3);
circle10_y = 0*t + y_current;
circle10_z = -r11*cos(t)+z_current +r11 - r11*(1-cos(theta3));

x_current = x_current - r11*sin(theta3);
z_current = z_current - r11*(1-cos(theta3));

fplot3(circle10_x,circle10_y,circle10_z,[0,theta3])

% linear 3

l7_x = -t + x_current;
l7_y = 0*t + y_current;
l7_z = 0*t + z_current;

x_current = x_current - 100;

fplot3(l7_x,l7_y,l7_z,[0,100])