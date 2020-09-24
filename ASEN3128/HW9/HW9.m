%Sage Herrin, SID 106071909
%ASEN 3128, HW9, question 4
%Created:4/7/19
clear all;close all;clc
%Question 4a, calc. and validate values in SI units in table 6.7
%initialize relevant constants 
u0 = 235.9; %m/s
S = 511.0; %m^2
b = 59.64; %m
rho = .3045; %kg/m^3
%nondimensional derivative values from table 6.6
NDvals = [-.8771 -.2797 .1946;...
              0 -.3295 -.04073;...
              0 .304 -.2737];
Lat_dim_derivs = [.5*rho*u0*S*NDvals(1,1) .5*rho*u0*b*S*NDvals(1,2)...
                    .5*rho*u0*b*S*NDvals(1,3);...
               .25*rho*u0*b*S*NDvals(2,1) .25*rho*u0*b^2*S*NDvals(2,2)...
               .25*rho*u0*b^2*S*NDvals(2,3);
               .25*rho*u0*b*S*NDvals(3,1) .25*rho*u0*b^2*S*NDvals(3,2)...
               .25*rho*u0*b^2*S*NDvals(3,3);];
               
%4bi, change in roll moment due to roll rate change delta_p = 0.05 rad/s
delta_L1 = Lat_dim_derivs(2,2)*0.05;
%4bii, change in yaw moment due to sudden yaw rate change delta_r = -0.05
%rad/s
delta_N1 = Lat_dim_derivs(3,3)*-0.05;
%4biii, change in roll moment due to sudden yaw rate change delta_r =
%0.01rad/s
delta_L2 = Lat_dim_derivs(3,2)*0.01;
%4biv, change in yaw moment due to sudden roll rate change of delta_p =
%-0.7rad/s
delta_N2 = Lat_dim_derivs(2,3)*-0.7;
%4bv, change in side force due to sudden roll rate change of delta_p = 0.15
%rad/s and side velocity delta_y = 2.04 m/s
delta_Y = Lat_dim_derivs(1,1)*2.04 + Lat_dim_derivs(2,1)*.15;
%rbvi, change in yaw moment due to simultaneous sudden changes in side
%velocity delta_y = -1.3 rad/s, roll rate delta_p = 0.5 rad/s, and yaw rate
%delta_r = 0.37 rad/s
delta_N3 = Lat_dim_derivs(1,3)*-1.3 + Lat_dim_derivs(2,3)*0.5 + Lat_dim_derivs(3,3)*0.37;



