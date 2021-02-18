clear all;close all;clc
%% Deriving the EOMs of the satellite under influence of drag and J2 -- using drag eqn from TSB pg 67
% assumes that station positions are given in ECI frame to match the frame that the S/C state is expressed in 
syms x y z xdot ydot zdot xs1 ys1 zs1 xs2 ys2 zs2 xs3 ys3 zs3 mu J2 Re CD We A m_sc
r = sqrt(x^2 + y^2 + z^2);
%density calculations for drag equation
rho0 = 3.614e-4;
r0 = 700 + Re;
H = 88.6670;
rho = rho0*exp(-(r - r0)/H);
%calculate relative velocity between sat and atm.
atm_vel = cross([0 0 We],[x y z]);
V_rvec = [xdot ydot zdot] - atm_vel;
Vr = norm(V_rvec);

R = [x y z];
Rs1 = [xs1 ys1 zs1];
Rs2 = [xs2 ys2 zs2];
Rs3 = [xs3 ys3 zs3];

Vs1 = cross([0 0 We],Rs1);
Vs2 = cross([0 0 We],Rs2);
Vs3 = cross([0 0 We],Rs3);

V = [xdot ydot zdot];

F_drag = -.5*rho*((CD*A)/m_sc)*Vr*V_rvec;

U = (mu/r)*(1 - ((J2*Re^2)/(2*r^4))*(3*z^2 - r^2)); 

rdd = simplify(jacobian(U,[x y z]));

rdd = rdd + F_drag;

state = [R V mu J2 CD Rs1 Rs2 Rs3];

f = [V rdd 0 0 0 Vs1 Vs2 Vs3];

Amat = simplify(jacobian(f,state));
