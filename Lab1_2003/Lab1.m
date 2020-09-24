% Charles Puskar
% Sage Herrin
% ASEN 2003
% Lab 1

clear all; close all; clc;

% Date Created: January 17, 2017
% Date Modified: January 26, 2017

% initial conditions
g = 9.807; % m/s^2
h0 = 125; % m
z_current = h0;

% defines symbolic parameter t
syms t z

% defines velocity from current height
v = @(z) sqrt(2*g*(h0 - z));
a = @(z) diff(v(z),z);

%% PARAMETERIZING THE PATH

figure(1)

%%
% circle 1
r1 = 40;
%plot
circle1_x = @(t) 0*t;
circle1_y = @(t) r1*sin(t);
circle1_z = @(t) r1*cos(t)+z_current-r1;
%Calculating G forces
G1_z = @(t) -(v(circle1_z(t)).^2)/(g*r1) + cos(t);
G1_y = @(t) sin(t);
G1_z(linspace(0,pi/4,20))
G1_y(linspace(0,pi/4,20))
%Re-initializing position
x_current = 0;
y_current = r1*sin(pi/4);
z_current = z_current - r1*(1-cos(pi/4));
t1 = linspace(0,pi/4,1000);
%finiding length of path
path_length = vpaintegral(sqrt((diff(circle1_y,t))^2 + ....
    (diff(circle1_z,t))^2),t1(1),t1(end));
%graphing
scatter3(circle1_x(t1),circle1_y(t1),circle1_z(t1),1,v(circle1_z(t1)))
hold on

%%
% circle 2
r2 = 40;
%plot
circle2_x = @(t) 0*t + x_current;
circle2_y = @(t) -r2*cos(t)+y_current+r2*sin(pi/4);
circle2_z = @(t) -r2*sin(t)+z_current+r2-r1*(1-cos(pi/4));
%Calculating G forces
G2_z = @(t) (v(circle2_z(t)).^2)/(g*r2) - sin(t);
G2_y = @(t) cos(t);
fprintf('G2')
G2_z(linspace(pi/4,pi/2,20))
G2_y(linspace(pi/4,pi/2,20))
%Re-initializing position
y_current = y_current+r2*sin(pi/4);
z_current = z_current - r2*(1-cos(pi/4));
t2 = linspace(pi/4,pi/2,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(circle2_x,t))^2 + ...
    (diff(circle2_y,t))^2 + (diff(circle2_z,t))^2),t2(1),t2(end));
%graphing
scatter3(circle2_x(t2),circle2_y(t2),circle2_z(t2),1,v(circle2_z(t2)))

%%
% loop 1
r3 = 10;

%Straight transition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1_x = @(t) 0*t+x_current;
l1_y = @(t) r3*t + y_current;
l1_z = @(t) z_current + 0*t;
%Calculating G forces
G3_z = @(t) 1 + 0*t;
G3_y = @(t) 0*t;
fprintf('G3')
G3_z(linspace(0,1,20))
G3_y(linspace(0,1,20))
t3 = linspace(0,1,1000);
%Finding length of path
path_length = path_length + vpaintegral(sqrt((diff(l1_x,t))^2 + ...
    (diff(l1_y,t))^2 + (diff(l1_z,t))^2),t3(1),t3(end));
%graphing
scatter3(l1_x(t3),l1_y(t3),l1_z(t3),1,v(l1_z(t3)))

%Actual part of loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%circle 3
circle3_x = @(t) -2*cos(t/2)+2+x_current;
circle3_y = @(t) r3*sin(t)+y_current+r3;
circle3_z = @(t) -r3*cos(t)+z_current+r3;
%Calculating G forces
G4_x_1 = @(t) diff(circle3_x(t),t);
G4_z = @(t) (v(circle3_z(t)).^2)/(g*r3) + cos(t);
G4_y = @(t) sin(t);
G4_x = @(x) vpa(subs(G4_x_1,t,x));

fprintf('G4')
G4_z(linspace(0,2*pi,20))
G4_y(linspace(0,2*pi,20))
G4_x(linspace(0,2*pi,20))
%Re-initializing position
y_current = y_current+r3;
x_current = 4;
t4 = linspace(0,2*pi,1000);
%Finding length of path
path_length = path_length + vpaintegral(sqrt((diff(circle3_x,t))^2 +...
    (diff(circle3_y,t))^2 + (diff(circle3_z,t))^2),t4(1),t4(end));
%graphing
scatter3(circle3_x(t4),circle3_y(t4),circle3_z(t4),1,v(circle3_z(t4)))

%Straight Transition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l2_x = @(t) 0*t+x_current;
l2_y = @(t) r3*t+y_current;
l2_z = @(t) z_current + 0*t;
%Calculating G forces
G5_z = @(t) 1 + 0*t;
G5_y = @(t) 0*t;
fprintf('G5')
G5_z(linspace(0,1,20))
G5_y(linspace(0,1,20))
%Re-initializing position
y_current = y_current + r3;
t5 = linspace(0,1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(l2_x,t))^2 + ...
    (diff(l2_y,t))^2 + (diff(l2_z,t))^2),t5(1),t5(end));
%graphing
scatter3(l2_x(t5),l2_y(t5),l2_z(t5),1,v(l2_z(t5)))

% circle 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r4 = 10;
%Circle transition to following parabola
%plot
circle4_x = @(t) 0*t+x_current;
circle4_y = @(t) r4*sin(t)+y_current;
circle4_z = @(t) -r4*cos(t)+z_current+r4;
theta1 = pi/4;
%Calculating G forces
G6_z = @(t) (v(circle4_z(t)).^2)/(g*r4) + cos(t);
G6_y = @(t) sin(t);
fprintf('G6')
G6_z(linspace(0,theta1,20))
G6_y(linspace(0,theta1,20))
%Re-initializing position
y_current = y_current + r4*sin(theta1);
z_current = z_current + r4*(1-cos(theta1));
t6 = linspace(0,theta1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(circle4_x,t))^2 + ...
    (diff(circle4_y,t))^2 + (diff(circle4_z,t))^2),t6(1),t6(end));
%graphing
scatter3(circle4_x(t6),circle4_y(t6),circle4_z(t6),1,v(circle4_z(t6)))

%%
% parabola 1
%plot
parabola1_x = @(t) 0*t+x_current;
parabola1_y = @(t) v(z_current)*cos(theta1)*t+y_current;
parabola1_z = @(t) v(z_current)*sin(theta1)*t - 0.5*g*t.^2 + z_current;
%Calcuating G forces for Parabola
find_max1 = @(t) diff(parabola1_z(t),t);
t_crit1 = double(solve(find_max1(t)));
%derivative equations used to find roh
deriv_z1 = @(t) diff(parabola1_z(t),t);
deriv_y1 = @(t) diff(parabola1_y(t),t);
deriv_1 = @(t) deriv_z1(t)./(deriv_y1(t).^2);
deriv2_1 = @(t) diff(deriv_1(t),t);

para_theta1 = @(t) atan(deriv_z1(t)./deriv_y1(t));
para_rho1 = @(t) ((1+(deriv_z1(t)./deriv_y1(t)).^2).^(3/2))./deriv2_1(t);
%Calculating G forces
G7_z_1 = @(t) -(v(parabola1_z(t)).^2)/(g*para_rho1(t)) - cos(para_theta1(t));
G7_y_1 = @(t) vpa(subs(a,z,parabola1_z(t))).*deriv_z1(t)./g + sin(para_theta1(t));

G7_z = @(x) vpa(subs(G7_z_1,t,x));
G7_y = @(x) vpa(subs(G7_y_1,t,x));

fprintf('G7')
G7_z(linspace(0,2*t_crit1,20))
G7_y(linspace(0,2*t_crit1,20))
%Re-initializing position
y_current = y_current+v(z_current)*cos(theta1)*2*t_crit1;
t7 = linspace(0,2*t_crit1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(parabola1_x,t))^2 + ...
    (diff(parabola1_y,t))^2 + (diff(parabola1_z,t))^2),t7(1),t7(end));
%graphing
scatter3(parabola1_x(t7),parabola1_y(t7),parabola1_z(t7),1,v(parabola1_z(t7)))

%circle transition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% circle 5
r5 = 10;
%plot
circle5_x = @(t) 0*t+x_current;
circle5_y = @(t) -r5*sin(theta1-t)+y_current+r5*sin(theta1);
circle5_z = @(t) -r5*cos(theta1-t)+z_current-r5*(1-cos(theta1))+r5;
%Calculating G forces
G8_z = @(t) (v(circle5_z(theta1-t)).^2)/(g*r5) + cos(theta1-t);
G8_y = @(t) sin(theta1-t);
fprintf('G8')
G8_z(linspace(0,theta1,20))
G8_y(linspace(0,theta1,20))
%Re-initializing position
y_current = y_current + r5*sin(theta1);
z_current = z_current - r5*(1-cos(theta1));
t8 = linspace(0,theta1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(circle5_x,t))^2 + ...
    (diff(circle5_y,t))^2 + (diff(circle5_z,t))^2),t8(1),t8(end));
%graphing
scatter3(circle5_x(t8),circle5_y(t8),circle5_z(t8),1,v(circle5_z(t8)))

%%
% banked turn
bank_angle = pi/4;
%calulating radius based on the angle
r6 = v(z_current).^2*cot(bank_angle)/g;

%Straight Transition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length3 = 20;
l3_theta = @(t) bank_angle*t;
l3_x = @(t) 0*t+x_current;
l3_y = @(t) length3*t+y_current;
l3_z = @(t) z_current + 0*t;
%Calculating G forces
G9_z = @(t) cos(l3_theta(t));
G9_y = @(t) 0*t;
G9_x = @(t) sin(l3_theta(t));
fprintf('G9')
G9_z(linspace(0,1,20))
G9_y(linspace(0,1,20))
%Re-initializing position
y_current = y_current + length3;
t9 = linspace(0,1,1000);
%finding lenght of path
path_length = path_length + vpaintegral(sqrt((diff(l3_x,t))^2 + ...
    (diff(l3_y,t))^2 + (diff(l3_z,t))^2),t9(1),t9(end));
%graphing
scatter3(l3_x(t9),l3_y(t9),l3_z(t9),1,v(l3_z(t9)))

%CircleTransition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%circle 6
circle8_x = @(t) r6*cos(t)+x_current+r6;
circle8_y = @(t) r6*sin(t)+y_current;
circle8_z = @(t) 0*t + z_current;
%Calculating G forces
G10_z = @(t) (v(circle8_z(t)).^2)*sin(bank_angle)/(g*r6) + cos(bank_angle);
G10_y = @(t) 0*t;
G10_x = @(t) (v(circle8_z(t)).^2)*cos(bank_angle)/(g*r6) - sin(bank_angle);
fprintf('G1')
G10_z(linspace(0,pi,20))
G10_y(linspace(0,pi,20))
G10_x(linspace(0,pi,20))
%Re-initializing position
x_current = x_current + 2*r6;
t10 = linspace(0,pi,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(circle8_x,t))^2 +...
    (diff(circle8_y,t))^2 + (diff(circle8_z,t))^2),t10(1),t10(end));
%graphing
scatter3(circle8_x(t10),circle8_y(t10),circle8_z(t10),1,v(circle8_z(t10)))

length4 = 20;
%Straight Transition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l4_theta = @(t) bank_angle*(1-t);
l4_x = @(t) 0*t+x_current;
l4_y = @(t) -length4*t+y_current;
l4_z = @(t) z_current + 0*t;
%calculating G
G11_z = @(t) cos(l4_theta(t));
G11_y = @(t) 0*t;
G11_x = @(t) sin(l4_theta(t));
fprintf('G1')
G11_z(linspace(0,1,20))
G11_y(linspace(0,1,20))
%Re-initializing position
y_current = y_current - length4;
t11 = linspace(0,1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(l4_x,t))^2 + ...
    (diff(l4_y,t))^2 + (diff(l4_z,t))^2),t11(1),t11(end));
%graphing
scatter3(l4_x(t11),l4_y(t11),l4_z(t11),1,v(l4_z(t11)))

%circle transition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% circle 9
r7 = 30;
theta2 = pi/8;
%circle 7
circle9_x = @(t) 0*t+x_current;
circle9_y = @(t) -r7*sin(t)+y_current;
circle9_z = @(t) r7*cos(t)+z_current-r7;
%calculating G forces
G12_z = @(t) -(v(circle9_z(t)).^2)/(g*r7) + cos(t);
G12_y = @(t) sin(t);
fprintf('G1')
G12_z(linspace(0,theta2,20))
G12_y(linspace(0,theta2,20))
%Re-initializing position
y_current = y_current - r7*sin(theta2);
z_current = z_current - r7*(1-cos(theta2));
t12 = linspace(0,theta2,1000);
%length of path
path_length = path_length + vpaintegral(sqrt((diff(circle9_x,t))^2 +...
    (diff(circle9_y,t))^2 + (diff(circle9_z,t))^2),t12(1),t12(end));
%graphing
scatter3(circle9_x(t12),circle9_y(t12),circle9_z(t12),1,v(circle9_z(t12)))

%%
% parabola 2
%plot
parabola2_x = @(t) 0*t+x_current;
parabola2_y = @(t) -v(z_current)*cos(theta2)*t+y_current;
parabola2_z = @(t) -(v(z_current)*sin(theta2)*t - 0.5*g*t.^2) + z_current;
%Calculatin G forces of Parabola
deriv_z2 = @(t) diff(parabola2_z(t),t);
deriv_y2 = @(t) diff(parabola2_y(t),t);
deriv_2 = @(t) deriv_z2(t)./(deriv_y2(t).^2);
deriv2_2 = @(t) diff(deriv_2(t),t);

find_max2 = @(t) diff(parabola2_z(t),t);
t_crit2 = double(solve(find_max2(t)));

para_theta2 = @(t) atan(deriv_z2(t)./deriv_y2(t));
para_rho2 = @(t) ((1+(deriv_z2(t)./deriv_y2(t)).^2).^(3/2))./deriv2_2(t);

G13_z_1 = @(t) (v(parabola2_z(t)).^2)/(g*para_rho2(t)) + cos(para_theta2(t));
G13_y_1 = @(t) vpa(subs(a,z,parabola2_z(t))).*deriv_z2(t)./g + sin(para_theta2(t));

G13_z = @(x) vpa(subs(G13_z_1,t,x));
G13_y = @(x) vpa(subs(G13_y_1,t,x));

fprintf('G13')
G13_z(linspace(0,2*t_crit2,20))
G13_y(linspace(0,2*t_crit2,20))
%Re-initializing position
y_current = y_current-v(z_current)*cos(theta2)*2*t_crit2;
t13 = linspace(0,2*t_crit2,1000);
%finidng length of path
path_length = path_length + vpaintegral(sqrt((diff(parabola2_x,t))^2 + ....
    (diff(parabola2_y,t))^2 + (diff(parabola2_z,t))^2),t13(1),t13(end));
%graphing
scatter3(parabola2_x(t13),parabola2_y(t13),parabola2_z(t13),1,v(parabola2_z(t13)))

%%
% parabola 3
theta3 = pi/3;
%plot
parabola3_x = @(t) 0*t+x_current;
parabola3_y = @(t) -v(z_current)*cos(theta2)*t+y_current;
parabola3_z = @(t) (v(z_current)*sin(theta2)*t - 0.5*g*t.^2) + z_current;
%Calculating G forces of Parabola
deriv_z3 = @(t) diff(parabola3_z(t),t);
deriv_y3 = @(t) diff(parabola3_y(t),t);
deriv_3 = @(t) deriv_z3(t)./(deriv_y3(t).^2);
deriv2_3 = @(t) diff(deriv_3(t),t);

find_max3 = @(t) diff(parabola3_z(t),t);
t_crit3 = double(solve(find_max3(t)));

t_total3 = double(solve(deriv_y3(t)-deriv_z3(t).*cot(theta3)));

para_theta3 = @(t) atan(deriv_z3(t)./deriv_y3(t));
para_rho3 = @(t) ((1+(deriv_z3(t)./deriv_y3(t)).^2).^(3/2))./deriv2_3(t);

G14_z_1 = @(t) -(v(parabola3_z(t)).^2)/(g*para_rho3(t)) - cos(para_theta3(t));
G14_y_1 = @(t) vpa(subs(a,z,parabola3_z(t))).*deriv_z3(t)./g - sin(para_theta3(t));

G14_z = @(x) vpa(subs(G14_z_1,t,x));
G14_y = @(x) vpa(subs(G14_y_1,t,x));

fprintf('G14')
G14_z(linspace(0,t_total3,20))
G14_y(linspace(0,t_total3,20))
%Re-initializing position
y_current = y_current-v(z_current)*cos(theta2)*t_total3;
z_current = parabola3_z(t_total3);

t14 = linspace(0,t_total3,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(parabola3_x,t))^2 + ...
    (diff(parabola3_y,t))^2 + (diff(parabola3_z,t))^2),t14(1),t14(end));
%graphing
scatter3(parabola3_x(t14),parabola3_y(t14),parabola3_z(t14),1,v(parabola3_z(t14)))

%%
% straight 1
r8 = 60;
drop1 = z_current-r8*(1-cos(theta3));
%plot
l5_x = @(t) 0*t+x_current;
l5_y = @(t) -drop1*cot(theta3)*t+y_current;
l5_z = @(t) z_current - drop1*t;
%Calculating G forces
G15_z = @(t) cos(theta3) + 0*t;
G15_y = @(t) sin(theta3) + 0*t;

fprintf('G15')
G15_z(linspace(0,1,20))
G15_y(linspace(0,1,20))
%Re-initializing position
y_current = y_current - drop1*cot(theta3);
z_current = z_current - drop1;

t15 = linspace(0,1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(l5_x,t))^2 + ...
    (diff(l5_y,t))^2 + (diff(l5_z,t))^2),t15(1),t15(end));
%graphing
scatter3(l5_x(t15),l5_y(t15),l5_z(t15),1,v(l5_z(t15)))

%circle trnasition???%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Circle 8
circle10_x = @(t) x_current+0*t;
circle10_y = @(t) r8*sin(theta3-t)+y_current-r8*sin(theta3);
circle10_z = @(t) -r8*cos(theta3-t)+z_current+r8-r8*(1-cos(theta3));
%Calculating G forces
G16_z = @(t) (v(circle10_z(t)).^2)/(g*r8) + cos(theta3-t);
G16_y = @(t) sin(theta3-t);

fprintf('G16')
G16_z(linspace(0,theta3,20))
G16_y(linspace(0,theta3,20))
%Re-initializing position
y_current = y_current - r8*sin(theta3);
z_current = z_current - r8*(1-cos(theta3));

t16 = linspace(0,theta3,1000);
%Finding length of path
path_length = path_length + vpaintegral(sqrt((diff(circle10_x,t))^2 + ...
    (diff(circle10_y,t))^2 + (diff(circle10_z,t))^2),t16(1),t16(end));
%graphing
scatter3(circle10_x(t16),circle10_y(t16),circle10_z(t16),1,v(circle10_z(t16)))

%%
% brake
decel = 2;
length6 = v(z_current)./(2*decel);
%plot
l6_x = @(t) 0*t+x_current;
l6_y = @(t) -length6*t+y_current;
l6_z = @(t) z_current + 0*t;
%Calculating G forces
G17_z = @(t) 1 + 0*t;
G17_y = @(t) -decel + 0*t;

fprintf('G17')
G17_z(linspace(0,1,20))
G17_y(linspace(0,1,20))
%Re-initializing position
y_current = y_current - length6;

t17 = linspace(0,1,1000);
%finding length of path
path_length = path_length + vpaintegral(sqrt((diff(l6_x,t))^2 + ...
    (diff(l6_y,t))^2 + (diff(l6_z,t))^2),t17(1),t17(end));
%graphing
scatter3(l6_x(t17),l6_y(t17),l6_z(t17),1,v(l6_z(t17))-decel.*t17)

%%
% plot formatting
c = colorbar;
ylabel(c,'Speed [m/s]')
colormap(jet)

title('Roller Coaster')
xlabel('x-position [m]')
ylabel('y-position [m]')
zlabel('z-position [m]')

az = -45;
el = 30;

view(az,el)

axis equal

fprintf('The Roller Coaster has a total length of %.2f meters \n',path_length)

