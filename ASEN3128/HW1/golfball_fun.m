%Author:Sage Herrin, SID 106071909
%ASEN 3128 Homework 1
%Created:1/25/19
%Modified:1/27/19
function [dydt] = golfball_fun(t,y)
global g cd m D N_0 E_0 D_0 vN_0 vE_0 vD_0 rho area q_0 V_E0 R_E W 
vn = y(1); %assign velocity and position to columns in y
ve = y(2);
vd = y(3);
N = y(4);
E = y(5);
D = y(6);
% wn = y(7);
% we = y(8);
% wd = y(9);
% v_rel = [(y(1) - y(7)) (y(2) - y(8)) (y(3) - y(9))];
v_inertial = [y(1) y(2) y(3)]; %Define intertial velocity 
v_rel = v_inertial - W; %Calculate relative wind 
vrel_mag = norm(v_rel); %Determine magnitude of relative wind 
vrel_unit = v_rel/vrel_mag; %Form unit vector from relative wind 
f_drag = -.5*rho*vrel_mag^2*cd*area*vrel_unit; %Calculate magnitude of drage from relative wind and given ball geometry and parameters

dvndt = (f_drag(1))/m; %Split up drag N,E, and D components 
dvedt = (f_drag(2))/m;
dvddt = (f_drag(3) + m*g)/m;
 
dndt = vn; %Split up velocity N,E, and D components
dedt = ve;
dddt = vd;

% dwndt = 0;
% dwedt = 0;
% dwddt = 0;

dydt(1) = dvndt; %Assigns dydt columns 
dydt(2) = dvedt;
dydt(3) = dvddt;
dydt(4) = vn;
dydt(5) = ve;
dydt(6) = vd;

dydt = dydt'; %transpose dydt